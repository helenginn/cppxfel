/*
 * Indexer.cpp
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#include "Indexer.h"
#include "Vector.h"
#include <cmath>
#include "gaussianfit.h"
#include "Holder.h"
#include "parameters.h"
#include <map>
#include <iostream>
#include "misc.h"
#include <fstream>
#include "Panel.h"
#include "FileParser.h"

#define BIG_BANDWIDTH 0.015
#define DISTANCE_TOLERANCE 0.02
#define WAVELENGTH_TOLERANCE 0.0001

#define HROT_TOLERANCE 0.001
#define KROT_TOLERANCE 0.001
#define ANGLE_TOLERANCE 0.00005
#define SPOT_DISTANCE_TOLERANCE 5

double Indexer::intensityThreshold;
bool Indexer::absoluteIntensity;

Indexer::Indexer(Image *newImage, MatrixPtr matrix)
{
    spaceGroup = ccp4spg_load_by_standard_num(197);
    initialStep = 0.2;
    testWavelength = 0;
    testDistance = 0;
    testBandwidth = BIG_BANDWIDTH;
    testSpotSize = 0.0002;
    maxResolution = 1.4;
    search = 5;
    searchSize = 8;
    reference = NULL;
    hRot = 0;
    kRot = 0;
    lastTotal = 0;
    lastStdev = 0;
    expectedSpots = 30;
    intensityThreshold = INTENSITY_THRESHOLD;
    refinement = RefinementTypeOrientationMatrixEarly;
    image = newImage;
    absoluteIntensity = FileParser::getKey("ABSOLUTE_INTENSITY", false);
    this->matrix = matrix;
    
    //	calculateNearbyMillers(true);
}

double Indexer::getDetectorDistance()
{
    return image->getDetectorDistance();
}

double Indexer::getWavelength()
{
    return image->getWavelength();
}

void Indexer::dropMillers()
{
    nearbyMillers.clear();
    vector<MillerPtr>().swap(nearbyMillers);
}

bool Indexer::millerReachesThreshold(MillerPtr miller)
{
    double iSigI = miller->getRawIntensity() / miller->getCountingSigma();
    
    if (absoluteIntensity)
        return (miller->getRawIntensity() > intensityThreshold);
        
    return (iSigI > intensityThreshold);
}

void Indexer::getWavelengthHistogram(vector<double> &wavelengths,
                                     vector<int> &frequencies, LogLevel level)
{
    wavelengths.clear();
    frequencies.clear();
    
    double interval = testBandwidth / 8;
    
    double wavelength = image->getWavelength();
    vector<double> totals;
    
    ostringstream logged;
    logged << "Wavelength histogram for " << this->image->getFilename() << std::endl;
    
    for (double i = wavelength * (1 - testBandwidth);
         i < wavelength * (1 + testBandwidth); i += interval)
    {
        wavelengths.push_back(i);
        int frequency = 0;
        int total = 0;
        
        for (int j = 0; j < millers.size(); j++)
        {
            double ewald = millers[j]->getWavelength(hRot, kRot);
            
            if (ewald < i || ewald > i + interval)
                continue;
            
            bool strong = millerReachesThreshold(millers[j]);
            
            total++;
            if (strong)
                frequency++;
        }
        
        logged << i << "\t";
        
        for (int i=0; i < frequency; i++)
        {
            logged << ".";
        }
        
        logged << std::endl;
        
        frequencies.push_back(frequency);
        totals.push_back(total);
    }
    
    int strong = getTotalReflections();
    
    logged << std::endl;
    logged << "Total strong reflections: " << strong << std::endl;

    
    Logger::mainLogger->addStream(&logged, level);
}

void Indexer::calculateNearbyMillers(bool rough)
{
    MatrixPtr matrix = getMatrix();
    double wavelength = image->getWavelength();
    
    double minBandwidth = wavelength * (1 - testBandwidth * 2);
    double maxBandwidth = wavelength * (1 + testBandwidth * 2);
    
    double sphereThickness = FileParser::getKey("SPHERE_THICKNESS", 0.01);
    
    double minSphere = 1 / wavelength * (1 - sphereThickness);
    double maxSphere = 1 / wavelength * (1 + sphereThickness);
    
    nearbyMillers.clear();
    
    int maxMillers[3];
    
    MatrixPtr newMatrix = matrix->copy();
    
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    
    if (!(hRad == 0 && kRad == 0))
        newMatrix->rotate(hRad, kRad, 0);
    
    for (int i = 0; i < 3; i++)
    {
        vec testHKL = new_vector(i == 0, i == 1, i == 2);
    
        newMatrix->multiplyVector(&testHKL);
        
        double maxD = 1 / maxResolution;
        double hklLength = length_of_vector(testHKL);
        
        maxMillers[i] = abs(maxD / hklLength);
    }
    
    logged << "Integrating to maximum Miller indices: (" << maxMillers[0] << ", " << maxMillers[1] << ", " << maxMillers[2] << ")" << std::endl;
    
    int overRes = 0;
    
    for (int h = -maxMillers[0]; h < maxMillers[0]; h++)
    {
        for (int k = -maxMillers[1]; k < maxMillers[1]; k++)
        {
            for (int l = -maxMillers[2]; l < maxMillers[2]; l++)
            {
                if (ccp4spg_is_sysabs(spaceGroup, h, k, l))
                    continue;
                
                MillerPtr newMiller = MillerPtr(new Miller(NULL, h, k, l));
                
                newMiller->setSelf(newMiller);
                newMiller->setImageAndIndexer(image, this);
                
                if (rough == false)
                {
                    nearbyMillers.push_back(newMiller);
                    continue;
                }

                vec hkl = new_vector(h, k, l);
                newMatrix->multiplyVector(&hkl);
                
                if (hkl.l > 0)
                    continue;
                
                vec beam = new_vector(0, 0, -1 / wavelength);
                vec beamToMiller = vector_between_vectors(beam, hkl);
                
                double sphereRadius = length_of_vector(beamToMiller);
                
                if (sphereRadius < minSphere || sphereRadius > maxSphere)
                {
                    double res = 1 / length_of_vector(hkl);
                    
                    if (res > 7)
                    {
                        nearbyMillers.push_back(newMiller);
                        continue;
                    }
                    
                    double ewaldSphere = getEwaldSphereNoMatrix(hkl);
                    
                    if (ewaldSphere < minBandwidth || ewaldSphere > maxBandwidth)
                    {
                        continue;
                    }
                    
                    continue;
                }
                
                double res = length_of_vector(hkl);
                
                if (res > 1 / maxResolution)
                {
                    overRes++;
                    continue;
                }

                if (h == 0 && k == 0 && l == 0)
                    continue;
                
               nearbyMillers.push_back(newMiller);
            }
        }
    }
    
    logged << "Rejected " << overRes << " due to being over resolution edge of " << maxResolution << " Å." << std::endl;
    
    sendLog(LogLevelDetailed);
}

void Indexer::checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox)
{
    MatrixPtr matrix = getMatrix();
    double wavelength = image->getWavelength();
    double maxD = 1 / maxResolution;
    if (maxResolution == 0)
        maxD = FLT_MAX;
    
    if (testDistance != 0)
        image->setDetectorDistance(testDistance);
    
    if (testWavelength != 0)
    {
        image->setWavelength(testWavelength);
        wavelength = testWavelength;
    }
    
    millers.clear();
    
    double averageEwald = 0;
    
    logged << "Testing " << nearbyMillers.size() << " reflections close to the Ewald sphere." << std::endl;
    
    int cutResolution = 0;
    int partialityTooLow = 0;
    int unacceptableIntensity = 0;
    
    for (int i = 0; i < nearbyMillers.size(); i++)
    {
        MillerPtr miller = nearbyMillers[i];
        
        vec hkl = new_vector(miller->getH(), miller->getK(), miller->getL());
        matrix->multiplyVector(&hkl);
        
        double d = length_of_vector(hkl);
        if (d > maxD)
        {
            cutResolution++;
            continue;
        }
        
        miller->setMatrix(matrix);
        miller->recalculatePartiality(hRot, kRot, 0.03, testSpotSize,
                                      wavelength, bandwidth, 1.5);
        
        if (i == 0)
            logged << "Calculated partiality with parameters: hRot " << hRot << ", kRot " << kRot << ", spot size " << testSpotSize << ", wavelength " << wavelength << ", bandwidth " << bandwidth << std::endl;
        
        if (miller->getPartiality() <= 0)
        {
            partialityTooLow++;
            continue;
        }
        
        miller->setPartialityModel(PartialityModelScaled);
        miller->getWavelength(hRot, kRot);
        
        if (complexShoebox)
        {
            double initialBandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
            double initialRlpSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
            double initialMosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);
            
            miller->makeComplexShoebox(wavelength, initialBandwidth, initialMosaicity, initialRlpSize);
        }
        
        miller->integrateIntensity(hRot, kRot);
        
        double rawIntensity = miller->getRawIntensity();
        
        if (rawIntensity != rawIntensity || (int) rawIntensity == 0)
        {
            unacceptableIntensity++;
            continue;
        }
        
        averageEwald += miller->getWavelength();
        
        this->millers.push_back(miller);
    }
    
    logged << "Beyond resolution cutoff: " << cutResolution << std::endl;
    logged << "Partiality equal to 0: " << partialityTooLow << std::endl;
    logged << "Image pixels were masked/flagged: " << unacceptableIntensity << std::endl;
    logged << "Reflections accepted: " << millers.size() << std::endl;
    
    sendLog(LogLevelDetailed);
}

void Indexer::minimizeParameter(double *meanStep, double *param)
{
    double param_trials[3];
    double param_scores[3];
    
    int j = 0;
    double param_min_score = FLT_MAX;
    int param_min_num = 1;
    
    double bestParam = *param;
    
    ostringstream logged;
    logged << "Scores for " << image->getFilename() << ": ";
    
    for (double i = bestParam - *meanStep; j < 3; i += *meanStep)
    {
        *param = i;
        this->checkAllMillers(maxResolution, testBandwidth);
        param_scores[j] = score();
        logged << param_scores[j] << ", ";
        param_trials[j] = i;
        j++;
    }
    
    param_min_score = param_scores[1];
    
    for (int i = 0; i < 3; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }
    
    *param = param_trials[param_min_num];
    this->checkAllMillers(maxResolution, testBandwidth);
    
    logged << "chosen no. " << param_min_num << std::endl;
    
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
    
    if (param_min_num == 1)
        *meanStep /= 2;
}

void Indexer::minimizeTwoParameters(double *meanStep1, double *meanStep2,
                                    double *param1, double *param2)
{
    double param_trials1[9];
    double param_trials2[9];
    double param_scores[9];
    
    int j = 0;
    double param_min_score = FLT_MAX;
    int param_min_num = 4;
    
    double bestParam1 = *param1;
    double bestParam2 = *param2;
    
    for (double i = bestParam1 - *meanStep1; j < 3; i += *meanStep1)
    {
        int l = 0;
        
        for (double k = bestParam2 - *meanStep2; l < 3; k += *meanStep2)
        {
            *param1 = i;
            *param2 = k;
            this->checkAllMillers(maxResolution, testBandwidth);
            param_scores[j * 3 + l] = score();
            param_trials1[j * 3 + l] = i;
            param_trials2[j * 3 + l] = k;
            l++;
        }
        j++;
    }
    
    param_min_score = param_scores[4];
    
    for (int i = 0; i < 9; i++)
    {
        if ((param_scores[i] < param_min_score))
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }
        
    }
    
    *param1 = param_trials1[param_min_num];
    *param2 = param_trials2[param_min_num];
    
    if (param_min_num == 4)
    {
        *meanStep1 /= 2;
        *meanStep2 /= 2;
    }
}

int Indexer::getTotalReflections(double threshold)
{
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (millers[i]->getRawIntensity() > threshold)
            count++;
    }
    
    return count;
}

int Indexer::getTotalReflections()
{
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (millerReachesThreshold(millers[i]))
            count++;
    }
    
    return count;
}

int Indexer::getTotalReflectionsWithinBandwidth()
{
    int count = 0;
    double minBandwidth = image->getWavelength() * (1 - testBandwidth);
    double maxBandwidth = image->getWavelength() * (1 + testBandwidth);
    
    for (int i = 0; i < millers.size(); i++)
    {
        double wavelength = millers[i]->getWavelength();
        
        if (wavelength < minBandwidth || wavelength > maxBandwidth)
            continue;
        
        if (millerReachesThreshold(millers[i]))
            count++;
    }
    
    return count;
}

double Indexer::getTotalIntegratedSignal()
{
    double totalIntensity = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (!isfinite(millers[i]->getRawIntensity()))
            continue;
        
        if (millers[i]->getRawIntensity() <= 0)
            continue;
        
        totalIntensity += millers[i]->getRawIntensity();
    }
    
    return totalIntensity;
}

void Indexer::findSpots()
{
    int tolerance = 60;
    
    //	image->printBox(1180, 660, 20);
    
    for (int i = 0; i < image->getXDim(); i += tolerance)
    {
        for (int j = 0; j < image->getYDim(); j += tolerance)
        {
            Spot *spot = new Spot(image);
            
            double maxLift = 0;
            double maxX = 0;
            double maxY = 0;
            
            for (int tolX = 0; tolX < tolerance; tolX++)
            {
                for (int tolY = 0; tolY < tolerance / 2; tolY++)
                {
                    double x = i + tolX;
                    double y = j + tolY;
                    
                    double lift = spot->maximumLift(image, x, y);
                    
                    if (lift > maxLift)
                    {
                        maxLift = lift;
                        maxX = x;
                        maxY = y;
                    }
                }
            }
            
            if (maxLift > 0)
            {
                spot->setXY(maxX, maxY);
                spots.push_back(spot);
                
                image->addSpotCover(maxX - 30, maxY - 30, maxX + 30, maxY + 30);
                //	image->printBox(maxX, maxY, 8);
            }
        }
    }
    
    Spot::sortSpots(&spots);
    
    std::string name = "spots-" + image->getFilename();
    int lastindex = (int)name.find_last_of(".");
    std::string rootName = name.substr(0, lastindex);
    std::string datName = rootName + ".dat";
    writeDatFromSpots(datName);
    
    ostringstream logged;
    logged << "Found " << spots.size() << " spots" << std::endl;
    Logger::mainLogger->addStream(&logged, LogLevelNormal);
}

void Indexer::duplicateSpots(std::vector<Image *> images)
{
    std::map<std::vector<int>, int> frequencies =
    std::map<std::vector<int>, int>();
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->indexerCount(); j++)
        {
            IndexerPtr indexer = images[i]->getIndexer(j);
            
            std::vector<Spot *> spots = indexer->getSpots();
            
            std::cout << "Image: " << i << " Spot count: " << spots.size()
            << std::endl;
            
            for (int j = 0; j < spots.size(); j++)
            {
                std::vector<int> coord = std::vector<int>();
                coord.push_back(spots[j]->x);
                coord.push_back(spots[j]->y);
                
                if (frequencies.count(coord) == 0)
                {
                    frequencies[coord] = 1;
                }
                else
                {
                    frequencies[coord]++;
                }
            }
        }
    }
    
    int threshold = 0.2 * images.size();
    int spotsRemoved = 0;
    
    if (images.size() < 4)
        return;
    
    for (map<std::vector<int>, int>::iterator it = frequencies.begin();
         it != frequencies.end(); ++it)
    {
        std::cout << it->first[0] << "\t" << it->first[1] << "\t"
        << frequencies[it->first] << std::endl;
        
        if (frequencies[it->first] >= threshold)
        {
            int startX = it->first[0] - 2;
            int endX = it->first[0] + 2;
            int startY = it->first[1] - 2;
            int endY = it->first[1] + 2;
            
            Image::applyMaskToImages(images, startX, startY, endX, endY);
            spotsRemoved++;
        }
    }
    
    std::cout << "Spots removed: " << spotsRemoved << std::endl;
}

void Indexer::scatterSpots(std::vector<Image *> images)
{
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->indexerCount(); j++)
        {
            IndexerPtr indexer = images[i]->getIndexer(j);
            
            std::vector<Spot *> spots = indexer->getSpots();
            
            for (int j = 0; j < spots.size(); j++)
            {
                std::cout << spots[j]->scatteringAngle(images[i]) << std::endl;
            }
        }
    }
}

void Indexer::writeDatFromSpots(std::string filename)
{
    std::ofstream dat;
    dat.open(filename);
    int count = 0;
    
    for (int i = 0; i < expectedSpots && i < spots.size(); i++)
    {
        //	if (spots[i]->isAcceptable(image))
        //	{
        dat << "0\t0\t0\t1\t1\t1\t" << spots[i]->x << "\t" << spots[i]->y
        << "\t1" << std::endl;
        
        std::cout << "Spot lift " << spots[i]->weight() << std::endl;
        
        count++;
        //	}
    }
    
    std::cout << count << " spots remain on image." << std::endl;
    
    dat.close();
}


int Indexer::identicalSpotsAndMillers()
{
    int count = 0;
    
    for (int i = 0; i < expectedSpots && i < spots.size(); i++)
    {
        for (int j = 0; j < millers.size(); j++)
        {
            int x1 = spots[i]->x;
            int x2 = millers[j]->getLastX();
            
            int y1 = spots[i]->y;
            int y2 = millers[j]->getLastY();
            
            int diffX = abs(x2 - x1);
            int diffY = abs(y2 - y1);
            
            if (diffX < 4 && diffY < 4)
                count++;
        }
    }
    
    return count;
}

double Indexer::medianIntensity()
{
    vector<double> intensities = vector<double>();
    
    for (int i = 0; i < millers.size(); i++)
    {
        intensities.push_back(millers[i]->getRawIntensity());
    }
    
    std::sort(intensities.begin(), intensities.end());
    
    int mid = (int)intensities.size() / 2;
    
    if (intensities.size() % 2 == 0)
    {
        return intensities[mid];
    }
    else
    {
        return (intensities[mid] + intensities[mid + 1]) / 2;
    }
    
}

double Indexer::score()
{
    if (refinement == RefinementTypeDetectorWavelength)
        return 0 - getTotalReflections();
    
    if (refinement == RefinementTypeOrientationMatrixEarly)
    {
        vector<double> wavelengths;
        vector<int> frequencies;
        
        getWavelengthHistogram(wavelengths, frequencies);
        
        double mean = 0;
        double stdev = 0;
        histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
        
        double total = getTotalReflectionsWithinBandwidth();
        //	double totalIntensity = getTotalIntegratedSignal();
        
        double totalWeight = 19;
        double stdevWeight = 15;
        
        double totalChange = total / lastTotal - 1;
        double totalStdev = stdev / lastStdev - 1;
        
        double score = (totalStdev * stdevWeight - totalChange * totalWeight);
        
     //   return 0 - total;

        return score;
    }
    
    if (refinement == RefinementTypeOrientationMatrixSpots)
    {
        double score = 0;
        checkAllMillers(maxResolution, testBandwidth);
        
        for (int j = 0; j < expectedSpots && j < spots.size(); j++)
        {
            for (int i = 0; i < millers.size(); i++)
            {
                double millerAngle = millers[i]->scatteringAngle(image);
                double spotAngle = spots[j]->scatteringAngle(image);
                
                if (abs(spotAngle - millerAngle) < ANGLE_TOLERANCE)
                {
                    spots[j]->setParentImage(image);
                    score += spots[j]->weight();
                    j++;
                    i = 0;
                }
            }
        }
        
        //	score /= millers.size();
        
        return 0 - score;
    }
    
    if (refinement == RefinementTypeOrientationMatrixExactSpots)
    {
        double score = 0;
        checkAllMillers(maxResolution, testBandwidth);
        
        for (int j = 0; j < expectedSpots && j < spots.size(); j++)
        {
            for (int i = 0; i < millers.size(); i++)
            {
                double millerX = millers[i]->getLastX();
                double millerY = millers[i]->getLastY();
                
                double spotX = spots[j]->x;
                double spotY = spots[j]->y;
                
                double distance = sqrt(pow(spotY - millerY, 2) + pow(spotX - millerX, 2));
                
                if (distance < SPOT_DISTANCE_TOLERANCE)
                {
                    spots[j]->setParentImage(image);
                    score += spots[j]->weight();
                    j++;
                    i = 0;
                }
            }
        }
        
        return 0 - score;
    }
    
    if (refinement == RefinementTypeOrientationMatrixMedian)
    {
        return 0 - identicalSpotsAndMillers();
    }
    
    if (refinement == RefinementTypeOrientationMatrixRough)
    {
        std::vector<double> wavelengths;
        
        for (int i = 0; i < millers.size(); i++)
        {
            if (millerReachesThreshold(millers[i]))
            {
                wavelengths.push_back(millers[i]->getWavelength());
            }
        }
        
        double stdev = standard_deviation(&wavelengths);
        double num = wavelengths.size();
        
        double newScore = stdev / num;
        
        return newScore;
    }
    
    if (refinement == RefinementTypeOrientationMatrixLate)
    {
        vector<double> wavelengths;
        vector<int> frequencies;
        
        getWavelengthHistogram(wavelengths, frequencies);
        
        double leastSquares = least_squares_gaussian_fit(&wavelengths,
                                                        &frequencies);
        
        return leastSquares;
    }
    
    return 0;
}

void Indexer::refineDetectorAndWavelength(MtzManager *reference)
{
    int oldSearch = getSearchSize();
    setSearchSize(0);
    this->reference = reference;
    this->calculateNearbyMillers(true);
    checkAllMillers(maxResolution, testBandwidth);
    refinement = RefinementTypeDetectorWavelength;
    
    testDistance = image->getDetectorDistance();
    testWavelength = image->getWavelength();
    double oldDistance = testDistance;
    double oldWavelength = testWavelength;
    
    int count = 0;
    double distStep = 0.2;
    double waveStep = 0.01;
    
    bool refinedDist = true;
    bool refinedWave = false;
    
    vector<double> wavelengths;
    vector<int> frequencies;
    
    double mean = 0;
    double stdev = 0;
    getWavelengthHistogram(wavelengths, frequencies);
    histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
    
    double oldScore = score();
    double newScore = 0;
    
    while (!(refinedDist && refinedWave) && count < 50)
    {
        if (!refinedDist)
            this->minimizeParameter(&distStep, &testDistance);
        
        if (!refinedWave)
            this->minimizeParameter(&waveStep, &testWavelength);
        
     //   this->minimizeTwoParameters(&distStep, &waveStep, &testDistance, &testWavelength);
        
        newScore = score();
        
        sendLog(LogLevelNormal);

        if (distStep < DISTANCE_TOLERANCE)
            refinedDist = true;
        
        if (waveStep < WAVELENGTH_TOLERANCE)
            refinedWave = true;
    }
    
    getWavelengthHistogram(wavelengths, frequencies);
    histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
    
 //   testWavelength = mean;
 //   this->image->setWavelength(mean);
    
    logged << "Distance refined for " << image->getFilename() << ":\t" << oldDistance << "\t" << oldWavelength << "\t" << testDistance << "\t" << testWavelength << "\t" << oldScore << "\t" << newScore << std::endl;
    
    sendLog(LogLevelNormal);
    
    image->setDetectorDistance(testDistance);
    
    setSearchSize(oldSearch);
}

// in degrees, returns degree rotation
double Indexer::refineRoundBeamAxis(double start, double end, double wedge,
                                   bool allSolutions)
{
    double startRad = start * M_PI / 180;
    vector<double> scores = vector<double>();
    vector<double> wedges = vector<double>();
    vector<double> hRads = vector<double>();
    vector<double> kRads = vector<double>();
    vector<MatrixPtr> matrices = vector<MatrixPtr>();
    
    refinement = RefinementTypeOrientationMatrixExactSpots;
    
    MatrixPtr copyMatrix = this->getMatrix()->copy();
    
    int solutionCount = allSolutions ? (int)solutions.size() : 1;
    
    for (int j = 0; j < solutionCount; j++)
    {
        double hRad = allSolutions ? solutions[j][0] : 0;
        double kRad = allSolutions ? solutions[j][1] : 0;
        
        for (double i = start; i < end; i += wedge)
        {
            double radians = i * M_PI / 180;
            
            checkAllMillers(maxResolution, testBandwidth);
            this->getMatrix()->rotate(hRad, kRad, radians);
            
            double total = score();
            hRads.push_back(hRad);
            kRads.push_back(kRad);
            scores.push_back(total);
            wedges.push_back(radians);
            matrices.push_back(this->getMatrix()->copy());
            logged << hRad << "\t" << kRad << "\t" << i << "\t" << total
            << std::endl;
            sendLog(LogLevelNormal);
            
            this->setMatrix(copyMatrix);
        }

    }
    
    double bestScore = 0;
    double bestWedge = 0;
    double bestHRad = 0;
    double bestKRad = 0;
    MatrixPtr bestMatrix = MatrixPtr();
    
    for (int i = 0; i < scores.size(); i++)
    {
        if (scores[i] < bestScore)
        {
            bestScore = scores[i];
            bestWedge = wedges[i];
            bestHRad = hRads[i];
            bestKRad = kRads[i];
            bestMatrix = matrices[i];
        }
    }
    
    std::cout << "Rotating to best score for " << bestWedge
    << "º rotation for " << bestHRad << ", " << bestKRad
    << std::endl;
    this->setMatrix(bestMatrix);
    std::cout << "Score: " << score() << std::endl;
    
    sendLog(LogLevelNormal);

    
    return bestWedge;
}

void Indexer::refineRoundBeamAxis()
{
    refineRoundBeamAxis(0, 360, 10, true);
    double betterWedge = refineRoundBeamAxis(-6, 6, 1, false);
    
    double bestRadians = betterWedge * M_PI / 180;
    checkAllMillers(maxResolution, testBandwidth);
    
    double total = getTotalIntegratedSignal();
    total /= millers.size();
    
    std::cout << "Signal is: " << total << std::endl;
    
    this->getMatrix()->printDescription();
}

void Indexer::matchMatrixToSpots()
{
    matchMatrixToSpots(RefinementTypeOrientationMatrixSpots);
}

bool compareScore(pair<vector<double>, double> a, pair<vector<double>, double> b)
{
    return (a.second < b.second);
}

void Indexer::matchMatrixToSpots(RefinementType refinement)
{
    double wedge = 45; // degrees
    //	map<vector<double>, double> scores = map<vector<double>, double>();
    vector<pair<vector<double>, double> > scores;
    
    MatrixPtr copyMatrix = this->getMatrix()->copy();
    
    this->refinement = refinement;
    
    calculateNearbyMillers(false);
    
    for (double hr = 0; hr < 360; hr += wedge)
    {
        double hRad = hr * M_PI / 180;
        
        std::cout << hr / 360 * 100 << "% ..." << std::endl;
        
        for (double kr = 0; kr < 360; kr += wedge)
        {
            double kRad = kr * M_PI / 180;
            
            this->getMatrix()->rotate(hRad, kRad, 0);
            checkAllMillers(maxResolution, testBandwidth);
            
            double theScore = score();
            
            vector<double> rotation = vector<double>();
            rotation.push_back(hRad);
            rotation.push_back(kRad);
            
            pair<vector<double>, double> pair = std::make_pair(rotation,
                                                             theScore);
            scores.push_back(pair);
            
            std::cout << hRad << "\t" << kRad << "\t" << theScore << std::endl;
        }
        
        this->setMatrix(copyMatrix);
    }
    
    std::cout << "100%" << std::endl;
    
    vector<double> bestNum = vector<double>();
    
    std::sort(scores.begin(), scores.end(), compareScore);
    
    std::cout << std::endl << "Solutions to try: " << std::endl;
    
    for (int i = 0; i < 3; i++)
    {
        std::cout << scores[i].first[0] << "\t" << scores[i].first[1] << "\t"
        << scores[i].second << std::endl;
        
        solutions.push_back(scores[i].first);
    }
    
    std::cout << std::endl;
}

void Indexer::refineOrientationMatrix()
{
    int orientationScore = FileParser::getKey("ORIENTATION_SCORE", 0);
    RefinementType refinementType = (RefinementType)orientationScore;
    
    this->refineOrientationMatrix(refinementType);
}

void Indexer::refineOrientationMatrix(RefinementType refinementType)
{
    refinement = refinementType;
    this->calculateNearbyMillers(true);
    
    testDistance = image->getDetectorDistance();
    testWavelength = image->getWavelength();
    
    int count = 0;
    double hRotStep = initialStep;
    double kRotStep = initialStep;
    
    bool refinedH = false;
    bool refinedK = false;
    
    vector<double> wavelengths;
    vector<int> frequencies;
    
    checkAllMillers(maxResolution, testBandwidth);
    sendLog(LogLevelDetailed);
    getWavelengthHistogram(wavelengths, frequencies);
    
    double mean = 0;
    double stdev = 0;
    double theScore = 0;
    
    histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
    
    lastStdev = stdev;
    lastTotal = getTotalReflections();
    
    image->setWavelength(mean);
    
    bool recalculated = false;
    
    while (!(refinedH && refinedK) && count < 50)
    {
      //  if (count % 2 == 1)
      //      this->calculateNearbyMillers(true);
        
        this->minimizeTwoParameters(&hRotStep, &kRotStep, &hRot, &kRot);
        
        checkAllMillers(maxResolution, testBandwidth);
        
        getWavelengthHistogram(wavelengths, frequencies);
        histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
        lastStdev = stdev;
        lastTotal = getTotalReflections();
        
        double newScore = score();
        lastScore = newScore;
        
        logged << hRot << "\t" << kRot << "\t" << newScore << "\t" << hRotStep << std::endl;
        
        if (refinement == RefinementTypeOrientationMatrixLate)
        {
            image->setWavelength(mean);
        }
        
        if (hRotStep < 0.1 && kRotStep < 0.1 && !recalculated)
        {
            recalculated = true;
            this->calculateNearbyMillers(true);
    //      refinement = RefinementTypeOrientationMatrixLate;
        }
   
        if (hRotStep < HROT_TOLERANCE)
            refinedH = true;
        
        if (kRotStep < KROT_TOLERANCE)
            refinedK = true;
        
        count++;
    }
    
    logged << "Rotation result:\t" << image->getFilename() << "\t" << hRot
    << "\t" << kRot << "\t" << getTotalReflections() << "\t" << getLastScore() << std::endl;
    
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    
    getMatrix()->rotate(hRad, kRad, 0);
    
    bestHRot = hRot;
    bestKRot = kRot;
    
    hRot = 0;
    kRot = 0;
    
    checkAllMillers(maxResolution, testBandwidth);
    getWavelengthHistogram(wavelengths, frequencies, LogLevelNormal);
    gaussian_fit(wavelengths, frequencies, (int)wavelengths.size(), &mean, &stdev,
                 &theScore, true);
    
    sendLog(LogLevelNormal);
}

MtzPtr Indexer::newMtz(int index)
{
    calculateNearbyMillers(true);
    setSearch(searchSize);
    checkAllMillers(maxResolution, testBandwidth);
    vector<double> wavelengths;
    vector<int> frequencies;
    
    double mean = 0;
    double stdev = 0;
    double theScore = 0;
    
    getWavelengthHistogram(wavelengths, frequencies, LogLevelNormal);
    gaussian_fit(wavelengths, frequencies, (int)wavelengths.size(), &mean, &stdev,
                 &theScore, true);
    
    bool complexShoebox = FileParser::getKey("COMPLEX_SHOEBOX", false);
    
    checkAllMillers(maxResolution, testBandwidth, complexShoebox);
    
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    
    MatrixPtr newMat = getMatrix()->copy();
    
    newMat->rotate(hRad, kRad, 0);
    
    MtzPtr mtz = MtzPtr(new MtzManager());
    mtz->setWavelength(mean);
    mtz->setFilename(image->filenameRoot() + "_" + i_to_str(index) + ".mtz");
    mtz->setSpaceGroup(spaceGroup->spg_num);
    mtz->setUnitCell(unitCell);
    
    mtz->setMatrix(newMat);
    double distance = image->getDetectorDistance();
    mtz->setDetectorDistance(distance);
    
    char *hallSymbol = ccp4spg_symbol_Hall(spaceGroup);
    
    space_group _spaceGroup = space_group(hallSymbol);
    space_group_type spgType = space_group_type(_spaceGroup);
    asu asymmetricUnit = asu(spgType);
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        miller->incrementOverlapMask();
        miller->setMtzParent(&*mtz);
        
        int index = Holder::indexForReflection(miller->getH(), miller->getK(), miller->getL(),
                                               mtz->getLowGroup(), false);
        
        Holder *found = NULL;
        mtz->findHolderWithId(index, &found);
        
        Panel::addMillerToPanelArray(miller);
        
        if (found != NULL)
        {
            found->addMiller(miller);
            miller->setParent(found);
        }
        else
        {
            Holder *holder = new Holder();
            holder->setSpaceGroup(spaceGroup, spgType, asymmetricUnit);
            holder->addMiller(miller);
            holder->calculateResolution(&*mtz);
            miller->setParent(holder);
            mtz->addHolder(holder);
            mtz->sortLastHolder();
        }
    }
    
    this->sendLog(LogLevelDetailed);
    
    double cutoff = FileParser::getKey("SIGMA_RESOLUTION_CUTOFF", SIGMA_RESOLUTION_CUTOFF);
    
    if (cutoff != 0)
        mtz->cutToResolutionWithSigma(cutoff);
    
    string imgFilename = "img-" + image->filenameRoot() + "_" + i_to_str(index) + ".mtz";
    mtz->writeToFile(imgFilename, false, true);
    mtz->writeToDat();

    nearbyMillers.clear();

    return mtz;
}

Indexer::~Indexer()
{
    nearbyMillers.clear();
    vector<MillerPtr>().swap(nearbyMillers);
    
    millers.clear();
    vector<MillerPtr>().swap(millers);
    
    if (spaceGroup != NULL)
        ccp4spg_free(&spaceGroup);
}

void Indexer::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}


