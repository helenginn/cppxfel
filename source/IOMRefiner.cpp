/*
 * IOMRefiner.cpp
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#include "IOMRefiner.h"
#include "Vector.h"
#include <cmath>
#include "Reflection.h"
#include "parameters.h"
#include <map>
#include <iostream>
#include "misc.h"
#include <fstream>
#include "FileParser.h"
#include "Miller.h"
#include "Hdf5Crystal.h"
#include "RefinementStrategy.h"
#include "RefinementStepSearch.h"
#include "Detector.h"

#define BIG_BANDWIDTH 0.015
#define DISTANCE_TOLERANCE 0.01
#define WAVELENGTH_TOLERANCE 0.0001

#define ANGLE_TOLERANCE 0.0001
#define SPOT_DISTANCE_TOLERANCE 5

bool IOMRefiner::lowIntensityPenalty = false;

IOMRefiner::IOMRefiner(ImagePtr newImage, MatrixPtr matrix)
{
    image = newImage;
    int spgNum = FileParser::getKey("SPACE_GROUP", -1);

    if (spgNum > 0)
    {
        Reflection::setSpaceGroup(spgNum);
        spaceGroup = ccp4spg_load_by_standard_num(spgNum);
    }
    
    initialStep = FileParser::getKey("INITIAL_ORIENTATION_STEP", INITIAL_ORIENTATION_STEP);
    
    testWavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    testBandwidth = FileParser::getKey("OVER_PRED_BANDWIDTH",
                                       OVER_PRED_BANDWIDTH) / 2;
    testSpotSize = FileParser::getKey("OVER_PRED_RLP_SIZE",
                                      OVER_PRED_SPOT_SIZE);
    maxResolution = FileParser::getKey(
                                       "MAX_INTEGRATED_RESOLUTION", MAX_INTEGRATED_RESOLUTION);;
    minResolution = FileParser::getKey(
                                       "MIN_INTEGRATED_RESOLUTION", 0.0);
    
    if (!FileParser::hasKey("MIN_INTEGRATED_RESOLUTION") || !FileParser::hasKey("MAX_INTEGRATED_RESOLUTION"))
    {
        if (Detector::isActive())
        {
            double min, max;
            Detector::getMaster()->resolutionLimits(&min, &max, getImage()->getWavelength());
            
            if (!FileParser::hasKey("MIN_INTEGRATED_RESOLUTION"))
            {
                minResolution = 1 / min;
                if (min == 0) minResolution = 0;
            }
            
            if (!FileParser::hasKey("MAX_INTEGRATED_RESOLUTION"))
            {
                maxResolution = 1 / max;
                if (max == 0) maxResolution = 0;
            }
            
            logged << "Setting minimum resolution to " << minResolution << " Å and max resolution to " << maxResolution << " Å." << std::endl;
            sendLog(LogLevelDetailed);
        }
    }
    
    searchSize = FileParser::getKey("METROLOGY_SEARCH_SIZE",
                                    METROLOGY_SEARCH_SIZE);
    
    reference = NULL;
    lastMtz = MtzPtr();
    roughCalculation = FileParser::getKey("ROUGH_CALCULATION", true);
    hRot = 0;
    kRot = 0;
    lRot = 0;
    bestHRot = 0;
    bestKRot = 0;
    bestLRot = 0;
    lastTotal = 0;
    lastStdev = 0;
    expectedSpots = FileParser::getKey("EXPECTED_SPOTS", 30);
    lowIntensityPenalty = FileParser::getKey("LOW_INTENSITY_PENALTY", false);
    this->matrix = matrix;
    unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    orientationTolerance = FileParser::getKey("INDEXING_ORIENTATION_TOLERANCE", INDEXING_ORIENTATION_TOLERANCE);
    needsReintegrating = true;
    recalculateMillerPositions = false;
    
    complexUnitCell = false;
}

double IOMRefiner::getDetectorDistance()
{
    return getImage()->getDetectorDistance();
}

double IOMRefiner::getWavelength()
{
    return getImage()->getWavelength();
}

void IOMRefiner::setComplexMatrix()
{
    complexUnitCell = true;
    bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);
    
    if (unitCell.size() == 0)
    {
        vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
        
        if (unitCell.size() == 0)
        {
            fixUnitCell = false;
        }
    }
    
    if (fixUnitCell)
    {
        getImage()->setUnitCell(unitCell);
        setUnitCell(unitCell);
        
        matrix->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
    }
}

void IOMRefiner::dropMillers()
{
    nearbyMillers.clear();
    vector<MillerPtr>().swap(nearbyMillers);
}

void IOMRefiner::getWavelengthHistogram(vector<double> &wavelengths,
                                        vector<int> &frequencies, LogLevel level, int whichAxis)
{
    bool shouldPrint = (Logger::getPriorityLevel() >= level);
    
    wavelengths.clear();
    frequencies.clear();
    
    if (shouldPrint)
    {
        logged << "Wavelength histogram for " << this->getImage()->getFilename() << std::endl;
        sendLog(level);
    }
    
    double wavelength = getImage()->getWavelength();
    vector<double> totals;
    
    vector<double> wavelengthRange = FileParser::getKey("WAVELENGTH_RANGE", vector<double>(2, 0));
    
    double spread = testBandwidth * 2;
    double interval = (wavelength * spread * 2) / 20;
    
    double minLength = wavelength * (1 - spread);
    double maxLength = wavelength * (1 + spread);
    
    if (wavelengthRange[0] != 0 || wavelengthRange[1] != 0)
    {
        minLength = wavelengthRange[0];
        maxLength = wavelengthRange[1];
        interval = (maxLength - minLength) / 20;
    }
    
    for (double i = minLength; i < maxLength; i += interval)
    {
        wavelengths.push_back(i);
        double frequency = 0;
        double total = 0;
        
        for (int j = 0; j < millers.size(); j++)
        {
            double ewald = millers[j]->getWavelength();
            
            if (ewald < i || ewald > i + interval)
                continue;
            
            bool strong = millers[j]->reachesThreshold();
            
            double weight = 1;
            
            total += weight;
            if (strong)
            {
                frequency += weight;
            }
            else if (lowIntensityPenalty)
            {
                frequency -= weight / 4;
            }
            
            if (frequency < 0)
            {
                frequency = 0;
            }
        }
        
        if (shouldPrint)
        {
            logged << std::setprecision(4) << i << "\t";
            
            for (int i=0; i < frequency; i++)
            {
                logged << ".";
            }
        }
        
        if (shouldPrint)
        {
            logged << std::endl;
        }
        
        frequencies.push_back(frequency);
        totals.push_back(total);
    }
    
    int strong = getTotalReflections();
    
    if (shouldPrint)
    {
        logged << std::endl;
        logged << "Total strong reflections: " << strong << std::endl;
        
        Logger::mainLogger->addStream(&logged, level);
    }
}

void IOMRefiner::calculateNearbyMillers(bool rough)
{
    double wavelength = getImage()->getWavelength();
    
    double minBandwidth = wavelength * (1 - testBandwidth * 2);
    double maxBandwidth = wavelength * (1 + testBandwidth * 2);
    
    double sphereThickness = FileParser::getKey("SPHERE_THICKNESS", 0.01);
    
    double minSphere = 1 / wavelength * (1 - sphereThickness);
    double maxSphere = 1 / wavelength * (1 + sphereThickness);
    
    double minSphereSquared = pow(minSphere, 2);
    double maxSphereSquared = pow(maxSphere, 2);
    
    nearbyMillers.clear();
    
    int maxMillers[3];
    
    Miller::rotateMatrixHKL(hRot, kRot, lRot, matrix, &lastRotatedMatrix);
    
    lastRotatedMatrix->maxMillers(maxMillers, maxResolution);
    
    logged << "Integrating to maximum Miller indices: (" << maxMillers[0] << ", " << maxMillers[1] << ", " << maxMillers[2] << ")" << std::endl;
    
    int overRes = 0;
    int underRes = 0;
    double maxDistSquared = pow(1 / maxResolution, 2);
    double minResolutionSquared = pow(1 / minResolution, 2);
    
    for (int h = -maxMillers[0]; h < maxMillers[0]; h++)
    {
        for (int k = -maxMillers[1]; k < maxMillers[1]; k++)
        {
            for (int l = -maxMillers[2]; l < maxMillers[2]; l++)
            {
                if (ccp4spg_is_sysabs(spaceGroup, h, k, l))
                    continue;
                
                vec hkl = new_vector(h, k, l);
                lastRotatedMatrix->multiplyVector(&hkl);
                
                if (hkl.l > 0.01)
                    continue;
                
                vec beam = new_vector(0, 0, -1 / wavelength);
                vec beamToMiller = vector_between_vectors(beam, hkl);
                
                double sphereRadiusSquared = length_of_vector_squared(beamToMiller);
                
                if (sphereRadiusSquared < minSphereSquared || sphereRadiusSquared > maxSphereSquared)
                {
                    double ewaldSphere = getEwaldSphereNoMatrix(hkl);
                    
                    if (ewaldSphere < minBandwidth || ewaldSphere > maxBandwidth)
                    {
                        continue;
                    }
                }
                
                double res = length_of_vector_squared(hkl);
                
                if (res > maxDistSquared)
                {
                    overRes++;
                    continue;
                }
                
                if (minResolution > 0 && res < minResolutionSquared)
                {
                    underRes++;
                    continue;
                }
                
                if (h == 0 && k == 0 && l == 0)
                    continue;
                
                MillerPtr newMiller = MillerPtr(new Miller(NULL, h, k, l));
                newMiller->setImageAndIOMRefiner(getImage(), shared_from_this());
                newMiller->setMatrix(lastRotatedMatrix);
                nearbyMillers.push_back(newMiller);
            }
        }
    }
    
    logged << "Rejected " << overRes << " due to being over resolution edge of " << maxResolution << " Å." << std::endl;
    
    if (minResolution > 0)
        logged << "Rejected " << underRes << " due to being under resolution edge of " << minResolution << " Å." << std::endl;
    
    sendLog(LogLevelDetailed);
    needsReintegrating = true;
}

/* merge with others... */
void IOMRefiner::lockUnitCellDimensions()
{
    double spgNum = spaceGroup->spg_num;
    if (spgNum >= 75 && spgNum <= 194)
    {
        unitCell[1] = unitCell[0];
    }
    if (spgNum >= 195)
    {
        unitCell[1] = unitCell[0];
        unitCell[2] = unitCell[0];
    }
}

void IOMRefiner::checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox, bool perfectCalculation)
{
    MatrixPtr matrix = getMatrix();
    
    if (complexUnitCell)
    {
        lockUnitCellDimensions();
        matrix->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
    }
    
    millers.clear();
    
    double wavelength = getImage()->getWavelength();
    
    double maxD = 1 / maxResolution;
    if (maxResolution == 0)
        maxD = FLT_MAX;
    
    double averageEwald = 0;
    
    logged << "Testing " << nearbyMillers.size()
    << " reflections close to the Ewald sphere with wavelength "
    << wavelength << std::endl;
    
    int cutResolution = 0;
    int partialityTooLow = 0;
    int unacceptableIntensity = 0;
    
    logged << "Rotating by " << hRot << ", " << kRot << ", " << lRot << std::endl;
    sendLog(LogLevelDetailed);
    
    Miller::rotateMatrixHKL(hRot, kRot, lRot, matrix, &lastRotatedMatrix);
    
    std::vector<MillerPtr> *chosenMillerArray = &nearbyMillers;
    
    if (!perfectCalculation && roughCalculation && !needsReintegrating && roughMillers.size())
    {
        chosenMillerArray = &roughMillers;
    }
    else if (needsReintegrating)
    {
        bandwidth *= 2;
        roughMillers.clear();
        std::vector<MillerPtr>().swap(roughMillers);
    }
    
    logged << "Checking " << chosenMillerArray->size() << " reflections." << std::endl;
    sendLog(LogLevelDetailed);
    
    for (int i = 0; i < chosenMillerArray->size(); i++)
    {
        MillerPtr miller = (*chosenMillerArray)[i];
        miller->setMatrix(lastRotatedMatrix);
        
        vec hkl = new_vector(miller->getH(), miller->getK(), miller->getL());
        matrix->multiplyVector(&hkl);
        
        double d = length_of_vector(hkl);
        if (d > maxD)
        {
            cutResolution++;
            continue;
        }
        
        miller->crossesBeamRoughly(lastRotatedMatrix, 0.0, testSpotSize, wavelength, bandwidth);
        
        if (i == 0)
        {
            logged << "Calculated partiality with parameters: hRot " << hRot << ", kRot " << kRot << ", spot size " << testSpotSize << ", wavelength " << wavelength << ", bandwidth " << bandwidth << std::endl;
            sendLog(LogLevelDebug);
        }
        
        if (miller->getPartiality() <= 0.05)
        {
            partialityTooLow++;
            continue;
        }
        
        miller->getWavelength(lastRotatedMatrix);
        
        if (complexShoebox)
        {
            double initialBandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
            double initialRlpSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
            double initialMosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);
            
            miller->makeComplexShoebox(wavelength, initialBandwidth, initialMosaicity, initialRlpSize);
        }
        
        double rawIntensity = miller->getRawIntensity();
        
        bool intensityMissing = (rawIntensity != rawIntensity);
        
        if (!roughCalculation || intensityMissing || needsReintegrating || complexShoebox)
        {
            miller->integrateIntensity(lastRotatedMatrix);
        }

        if (recalculateMillerPositions)
        {
            int x, y;
            miller->positionOnDetector(lastRotatedMatrix, &x, &y, false);
        }
        
        rawIntensity = miller->getRawIntensity();
        
        if (rawIntensity != rawIntensity)
        {
            unacceptableIntensity++;
            chosenMillerArray->erase(chosenMillerArray->begin() + i);
            i--;
            continue;
        }
        
        averageEwald += miller->getWavelength();
        
        this->millers.push_back(miller);
        
        if (needsReintegrating && miller->reachesThreshold())
        {
            roughMillers.push_back(miller);
        }
    }
    
    logged << "Using wavelength " << wavelength << " Å and distance " << getImage()->getDetectorDistance() << " mm" << std::endl;
    logged << "Beyond resolution cutoff: " << cutResolution << std::endl;
    logged << "Partiality equal to 0: " << partialityTooLow << std::endl;
    logged << "Image pixels were masked/flagged: " << unacceptableIntensity << std::endl;
    logged << "Reflections accepted: " << millers.size() << std::endl;
    
    sendLog(LogLevelDetailed);
    
    needsReintegrating = false;
}

double IOMRefiner::getReflectionWavelengthStdev()
{
    std::vector<double> values;
    
    bool unbalanced = FileParser::getKey("UNBALANCED_REFLECTIONS", false);
    double mean = getWavelength();
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (millers[i]->reachesThreshold() && millerWithinBandwidth(millers[i]))
        {
            values.push_back(millers[i]->getWavelength());
        }
    }
    
    if (!unbalanced)
    {
        mean = weighted_mean(&values);
    }
    
    double stdev = standard_deviation(&values, NULL, mean);
    
    return stdev;
}

int IOMRefiner::getTotalReflections()
{
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (millers[i]->reachesThreshold())
            count++;
    }
    
    return count;
}

bool IOMRefiner::millerWithinBandwidth(MillerPtr miller)
{
    double minBandwidth = getImage()->getWavelength() * (1 - testBandwidth * 2);
    double maxBandwidth = getImage()->getWavelength() * (1 + testBandwidth * 2);
    
    if (FileParser::hasKey("WAVELENGTH_RANGE"))
    {
        std::vector<double> range = FileParser::getKey("WAVELENGTH_RANGE", std::vector<double>());
        
        if (range.size() >= 2)
        {
            minBandwidth = range[0];
            maxBandwidth = range[1];
        }
    }
    
    double wavelength = miller->getWavelength();
    
    return (wavelength > minBandwidth && wavelength < maxBandwidth);
}

int IOMRefiner::getTotalReflectionsWithinBandwidth()
{
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (!millerWithinBandwidth(millers[i]))
            continue;
        
        if (millers[i]->reachesThreshold())
            count++;
    }
    
    return count;
}

void IOMRefiner::duplicateSpots(vector<ImagePtr> images)
{
    std::map<vector<int>, int> frequencies =
    std::map<vector<int>, int>();
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->IOMRefinerCount(); j++)
        {
            IOMRefinerPtr IOMRefiner = images[i]->getIOMRefiner(j);
            
            vector<Spot *> spots = IOMRefiner->getSpots();
            
            std::cout << "Image: " << i << " Spot count: " << spots.size()
            << std::endl;
            
            for (int j = 0; j < spots.size(); j++)
            {
                vector<int> coord = vector<int>();
                coord.push_back(spots[j]->getX());
                coord.push_back(spots[j]->getY());
                
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
    
    for (std::map<vector<int>, int>::iterator it = frequencies.begin();
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

void IOMRefiner::recalculateMillers(bool thorough)
{
    this->checkAllMillers(maxResolution, testBandwidth, false, true);
}

double IOMRefiner::hkScoreFinalise(void *object)
{
    IOMRefiner *refiner = static_cast<IOMRefiner *>(object);
    refiner->recalculateMillers();
    return refiner->hkScore(false, true);
}

double IOMRefiner::hkScoreWrapper(void *object)
{
    IOMRefiner *refiner = static_cast<IOMRefiner *>(object);
    refiner->recalculateMillers();
    return refiner->hkScore(false);
}

double IOMRefiner::hkScoreStdevWrapper(void *object)
{
    IOMRefiner *refiner = static_cast<IOMRefiner *>(object);
    refiner->recalculateMillers();
    return refiner->hkScore(true);
}

double IOMRefiner::lScoreWrapper(void *object)
{
    IOMRefiner *refiner = static_cast<IOMRefiner *>(object);
    refiner->recalculateMillerPositions = true;
    refiner->recalculateMillers();
    return refiner->lScore(true);
}

double IOMRefiner::hkScore(bool reportStdev, bool finalise)
{
    vector<double> wavelengths;
    vector<double> throw1; vector<int> throw2;
    
    if (Logger::getPriorityLevel() >= LogLevelDetailed)
    {
        getWavelengthHistogram(throw1, throw2, LogLevelDetailed);
    }
    
    bool unbalanced = FileParser::getKey("UNBALANCED_REFLECTIONS", false);
    double mean = getWavelength();
    
    if (!unbalanced)
    {
        mean = weighted_mean(&wavelengths);
    }
    
    double stdev = getReflectionWavelengthStdev();
    double refTotal = getTotalReflections();
    
    if (finalise)
    {
        lastTotal = refTotal;
        lastStdev = stdev;
    }

    double totalWeight = 20;
    double stdevWeight = 15;
    
    double totalChange = refTotal / lastTotal - 1;
    double totalStdev = stdev / lastStdev - 1;
    
    double score = (totalStdev * stdevWeight - totalChange * totalWeight);
    
    logged << "Score is " << score << " from " << wavelengths.size() << " reflections" << std::endl;
    sendLog(LogLevelDetailed);
    
    return score;
}

double IOMRefiner::lScore(bool silent)
{
    double averageShift = 0;
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (millers[i]->reachesThreshold())
        {
            std::pair<double, double> shift = millers[i]->getShift();
            double shiftDistance = sqrt(pow(shift.first, 2) + pow(shift.second, 2));
            
            averageShift += shiftDistance;
            count++;
        }
    }
    
    averageShift /= count;
    
    logged << "Average shift: " << averageShift << std::endl;
    sendLog(LogLevelDetailed);
    
    return averageShift;
}

double IOMRefiner::getRot(int rotNum)
{
    switch (rotNum)
    {
        case 0:
            return hRot;
        case 1:
            return kRot;
        case 2:
            return lRot;
        default:
            return 0;
            break;
    }
}

void IOMRefiner::refineOrientationMatrix()
{
    this->calculateNearbyMillers();
    
    vector<double> wavelengths;
    vector<int> frequencies;
    
    checkAllMillers(maxResolution, testBandwidth);
    logged << "Wavelength histogram before refinement for " << getImage()->getFilename() << std::endl;
    sendLog(LogLevelDetailed);
    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);
    
    lastStdev = getReflectionWavelengthStdev();
    lastTotal = getTotalReflections();
    
    int oldSearchSize = searchSize;
    int bigSize = FileParser::getKey("METROLOGY_SEARCH_SIZE_BIG", 6);
    
    // FIXME: read support for unit cell dimensions
    
    recalculateMillerPositions = true;
    searchSize = bigSize;
    
    RefinementStrategyPtr lStrategy = RefinementStrategy::userChosenStrategy();
    lStrategy->setVerbose(Logger::getPriorityLevel() >= LogLevelDetailed);
    lStrategy->setEvaluationFunction(lScoreWrapper, this);
    lStrategy->setJobName("Refining in-detector-plane angle for " + getImage()->getFilename());
    lStrategy->addParameter(this, getLRot, setLRot, initialStep, orientationTolerance, "lRot");
    lStrategy->refine();

    searchSize = oldSearchSize;
    recalculateMillerPositions = false;

    if (initialStep >= 1.5)
    {
        hkScoreFinalise(this);
        
        RefinementStepSearchPtr hkStrategy = RefinementStepSearchPtr(new RefinementStepSearch());
        hkStrategy->setVerbose(Logger::getPriorityLevel() >= LogLevelDetailed);
        hkStrategy->setEvaluationFunction(hkScoreWrapper, this);
        hkStrategy->setAfterCycleFunction(hkScoreFinalise, this);
        hkStrategy->setJobName("Refining angles for " + getImage()->getFilename());
        hkStrategy->addParameter(this, getHRot, setHRot, initialStep, orientationTolerance, "hRot");
        hkStrategy->addCoupledParameter(this, getKRot, setKRot, initialStep, orientationTolerance, "kRot");
        hkStrategy->refine();
    }
    
    this->calculateNearbyMillers();

 //   if (false)
    {
        hkScoreFinalise(this);
        RefinementStepSearchPtr hkStrategy = RefinementStepSearchPtr(new RefinementStepSearch());
        hkStrategy->setVerbose(Logger::getPriorityLevel() >= LogLevelDetailed);
        hkStrategy->setEvaluationFunction(hkScoreStdevWrapper, this);
        hkStrategy->setAfterCycleFunction(hkScoreFinalise, this);
        hkStrategy->setJobName("Refining angles for " + getImage()->getFilename());
        hkStrategy->addParameter(this, getHRot, setHRot, std::max(0.5, initialStep), orientationTolerance, "hRot");
        hkStrategy->addCoupledParameter(this, getKRot, setKRot, std::max(0.5, initialStep), orientationTolerance, "kRot");
        hkStrategy->refine();
    }
    
    bestHRot = hRot;
    bestKRot = kRot;
    bestLRot = lRot;
    
    //needsReintegrating = true;
    //this->calculateNearbyMillers();
    checkAllMillers(maxResolution, testBandwidth);
    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);
    
    matrix = lastRotatedMatrix;
    
    hRot = 0;
    kRot = 0;
    lRot = 0;
}

/*
 void IOMRefiner::refineOrientationMatrix(RefinementType refinementType)
 {
 refinement = refinementType;
 this->calculateNearbyMillers(true);
 
 testWavelength = getImage()->getWavelength();
 
 vector<double> wavelengths;
 vector<int> frequencies;
 
 checkAllMillers(maxResolution, testBandwidth);
 Logger::mainLogger->addString("Wavelength histogram before refinement", LogLevelDetailed);
 sendLog(LogLevelDetailed);
 getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);
 
 double mean = 0;
 double stdev = 0;
 double theScore = 0;
 
 histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
 
 lastStdev = stdev;
 lastTotal = getTotalReflections();
 
 getImage()->setWavelength(mean);
 
 bool recalculated = false;
 
 for (int i = 0; i < 1; i++)
 {
 double hRotStep = initialStep;
 double kRotStep = initialStep;
 double lRotStep = initialStep;
 
 double alphaStep = FileParser::getKey("STEP_UNIT_CELL_ALPHA", 0.5);
 double betaStep = FileParser::getKey("STEP_UNIT_CELL_BETA", 0.5);
 double gammaStep = FileParser::getKey("STEP_UNIT_CELL_GAMMA", 0.5);
 
 bool refinedH = false;
 bool refinedK = false;
 bool refinedL = !FileParser::getKey("REFINE_IN_PLANE_OF_DETECTOR", true);
 bool refinedAlpha = !FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
 bool refinedBeta = !FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
 bool refinedGamma = !FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
 int bigSize = FileParser::getKey("METROLOGY_SEARCH_SIZE_BIG", 6);
 
 int count = 0;
 
 while (!(refinedH && refinedK && refinedL) && count < 20)
 {
 if (!refinedL)
 {
 int oldSearchSize = searchSize;
 
 if (searchSize < bigSize)
 {
 searchSize = bigSize;
 }
 
 recalculateMillerPositions = true;
 refinement = RefinementTypeRefineLAxis;
 
 this->minimizeParameter(&lRotStep, &lRot);
 refinement = refinementType;
 recalculateMillerPositions = false;
 
 searchSize = oldSearchSize;
 }
 if (!refinedH && !refinedK)
 {
 this->minimizeTwoParameters(&hRotStep, &kRotStep, &hRot, &kRot);
 }
 
 //     refinementType = RefinementTypeOrientationMatrixPanelStdev;
 
 if (!refinedAlpha)
 {
 minimizeParameter(&alphaStep, &unitCell[3]);
 }
 
 if (!refinedBeta)
 {
 minimizeParameter(&betaStep, &unitCell[4]);
 }
 
 if (!refinedGamma)
 {
 minimizeParameter(&gammaStep, &unitCell[5]);
 }
 
 //     refinementType = RefinementTypeOrientationMatrixEarly;
 
 checkAllMillers(maxResolution, testBandwidth, false, false);
 
 getWavelengthHistogram(wavelengths, frequencies);
 histogram_gaussian(&wavelengths, &frequencies, mean, stdev);
 lastStdev = stdev;
 lastTotal = getTotalReflections();
 
 double newScore = score();
 lastScore = newScore;
 
 logged << getRot(0) << "\t" << getRot(1) << "\t" << getRot(2) << "\t" << newScore << std::endl;
 sendLog(LogLevelDetailed);
 
 if (hRotStep < 0.25 && kRotStep < 0.25)
 {
 if (!recalculated)
 {
 recalculated = true;
 this->calculateNearbyMillers(true);
 }
 }
 
 if (hRotStep < orientationTolerance)
 refinedH = true;
 if (kRotStep < orientationTolerance)
 refinedK = true;
 
 if (lRotStep < orientationTolerance)
 refinedL = true;
 
 if (alphaStep < 0.01)
 refinedAlpha = true;
 if (betaStep < 0.01)
 refinedBeta = true;
 if (gammaStep < 0.01)
 refinedGamma = true;
 
 count++;
 }
 
 count = 0;
 }
 
 lastScore = score();
 
 logged << "Current wavelength: " << testWavelength << " Å." << std::endl;
 logged << "Rotation result:\t" << getImage()->getFilename() << "\t" << hRot
 << "\t" << kRot << "\t" << getTotalReflections() << "\t" << getLastScore() << std::endl;
 
 double hRad = getRot(0) * M_PI / 180;
 double kRad = getRot(1) * M_PI / 180;
 double lRad = getRot(2) * M_PI / 180;
 
 vector<double> originalUnitCell = FileParser::getKey("UNIT_CELL", vector<double>());
 double *lengths = new double[3];
 getMatrix()->unitCellLengths(&lengths);
 
 double aRatio = originalUnitCell[0] / lengths[0];
 double bRatio = originalUnitCell[1] / lengths[1];
 double cRatio = originalUnitCell[2] / lengths[2];
 
 double aveRatio = (aRatio + bRatio + cRatio) / 3;
 
 lengths[0] *= aveRatio;
 lengths[1] *= aveRatio;
 lengths[2] *= aveRatio;
 
 delete [] lengths;
 
 getMatrix()->rotate(hRad, kRad, lRad);
 
 bestHRot = getRot(0);
 bestKRot = getRot(1);
 bestLRot = getRot(2);
 
 hRot = 0;
 kRot = 0;
 lRot = 0;
 
 needsReintegrating = true;
 checkAllMillers(maxResolution, testBandwidth);
 getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);
 
 sendLog(LogLevelNormal);
 }
 */

bool IOMRefiner::isGoodSolution()
{
    if (getImage()->getClass() == ImageClassDIALS)
    {
        logged << "Solution refinement has not been enabled for DIALS images yet - assuming good solution." << std::endl;
        sendLog();
        
        return true;
    }
    
    bool good = false;
    double goodSolutionStdev = FileParser::getKey("GOOD_SOLUTION_ST_DEV", 0.04);
    double badSolutionStdev = FileParser::getKey("BAD_SOLUTION_ST_DEV", 0.066);
    double goodSolutionSumRatio = FileParser::getKey("GOOD_SOLUTION_SUM_RATIO", 6.5);
    int goodSolutionHighestPeak = FileParser::getKey("GOOD_SOLUTION_HIGHEST_PEAK", 17);
    int badSolutionHighestPeak = FileParser::getKey("BAD_SOLUTION_HIGHEST_PEAK", 10);
    int minimumReflections = FileParser::getKey("MINIMUM_REFLECTION_CUTOFF", 30);
    std::ostringstream details;
    
    logged << "(" << getImage()->getFilename() << ") Standard deviation: " << lastStdev << std::endl;
    sendLog(LogLevelNormal);
    
    vector<double> wavelengths;
    vector<int> frequencies;
    
    calculateOnce();
    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed);
    
    double totalMean = 0;
    double totalStdev = 0;
    int maxFrequency = 0;
    
    histogram_gaussian(&wavelengths, &frequencies, totalMean, totalStdev);
    
    for (int i = 0; i < frequencies.size(); i++)
    {
        if (frequencies[i] > maxFrequency)
        {
            maxFrequency = frequencies[i];
        }
    }
    
    double diffSquared = 0;
    double worstSquared = 0;
    double maxHeight = super_gaussian(totalMean, totalMean, goodSolutionStdev * 0.8, 2) ;
    
    for (int i = 0; i < frequencies.size(); i++)
    {
        double length = wavelengths[i];
        double expected = super_gaussian(length, totalMean, goodSolutionStdev * 0.8, 2) * (double)maxFrequency / maxHeight;
        
        //    logged << "expected: " << expected << ", frequency: " << frequencies[i] << std::endl;
        
        diffSquared += pow(expected - frequencies[i], 2);
        worstSquared += pow(expected, 2);
    }
    
    diffSquared = sqrt(diffSquared);
    worstSquared = sqrt(worstSquared);
    
    std::sort(frequencies.begin(), frequencies.end(), std::greater<int>());
    
    double highSum = 0;
    std::vector<double> lowOnly;
    int cutoff = 3;
    
    for (int i = 0; i < frequencies.size(); i++)
    {
        if (i < cutoff)
        {
            highSum += frequencies[i];
            continue;
        }
        
        lowOnly.push_back(frequencies[i]);
    }
    
    highSum /= cutoff;
    
    double stdevLow = standard_deviation(&lowOnly, NULL, 0);
    
    logged << "Stdev low: " << stdevLow << ", highAverage " << highSum << std::endl;
    sendLog();
    
    if (lastStdev < goodSolutionStdev)
    {
        good = true;
        details << "(" << getImage()->getFilename() << ") Standard deviation is sufficiently low (" << lastStdev << " vs " << goodSolutionStdev << ")" << std::endl;
    }
    
    if (highSum > stdevLow * goodSolutionSumRatio)
    {
        good = true;
        details << "(" << getImage()->getFilename() << ") Sum ratio is sufficiently high (" << highSum << " vs " << stdevLow << ")" << std::endl;
    }
    
    /*  if (highSum <= 5)
     {
     details << "(" << getImage()->getFilename() << ") However, high sum not high enough (" << highSum << ")" << std::endl;
     good = false;
     }
     */
    if (frequencies[0] > goodSolutionHighestPeak)
    {
        details << "(" << getImage()->getFilename() << ") Highest peak is high enough (" << frequencies[0] << " vs " << goodSolutionHighestPeak << ")" << std::endl;
        good = true;
    }
    
    if (frequencies[0] < badSolutionHighestPeak)
    {
        details << "(" << getImage()->getFilename() << ") Highest peak is too low (" << frequencies[0] << " vs " << badSolutionHighestPeak << ")" << std::endl;
        good = false;
    }
    
    //    if (lastScore < 3)
    //       good = true;
    
    if (getTotalReflections() < minimumReflections)
    {
        details << "(" << getImage()->getFilename() << ") However, not enough reflections (" << getTotalReflections() << " vs " << minimumReflections << ")" << std::endl;
        good = false;
    }
    
    if (lastStdev > badSolutionStdev)
    {
        good = false;
        details << "(" << getImage()->getFilename() << ") However, standard deviation too high (" << lastStdev << " vs " << badSolutionStdev << ")" << std::endl;
    }
    
    Logger::mainLogger->addStream(&details, LogLevelNormal);
    
    return good;
}

void IOMRefiner::calculateOnce()
{
    needsReintegrating = true;
    calculateNearbyMillers();
    checkAllMillers(maxResolution, testBandwidth);
    
    vector<double> wavelengths;
    vector<int> frequencies;
    
    getWavelengthHistogram(wavelengths, frequencies, LogLevelDetailed, 0);
}

void IOMRefiner::showHistogram(bool silent)
{
    vector<double> wavelengths;
    vector<int> frequencies;
    
    getWavelengthHistogram(wavelengths, frequencies, silent ? LogLevelDebug : LogLevelNormal, 0);
}

MtzPtr IOMRefiner::newMtz(int index, bool silent)
{
    calculateNearbyMillers();
    checkAllMillers(maxResolution, testBandwidth);
    vector<double> wavelengths;
    vector<int> frequencies;
    
    std::ostringstream logged;
    logged << "Wavelength histogram for " << this->getImage()->getFilename() << std::endl;
    sendLog(silent ? LogLevelDebug : LogLevelNormal);
    getWavelengthHistogram(wavelengths, frequencies, silent ? LogLevelDebug : LogLevelNormal, 0);
    
    needsReintegrating = true;
    bool complexShoebox = FileParser::getKey("COMPLEX_SHOEBOX", false);
    
    checkAllMillers(maxResolution, testBandwidth, complexShoebox);
    
    MatrixPtr newMat = lastRotatedMatrix->copy();
    
    MtzPtr mtz;
    
    std::string crystalFileName = getImage()->filenameRoot() + "_" + i_to_str(index) + ".mtz";
    
    if (getImage()->getClass() == ImageClassHdf5)
    {
        Hdf5CrystalPtr crystal = Hdf5CrystalPtr(new Hdf5Crystal(crystalFileName));
        
        mtz = boost::static_pointer_cast<MtzManager>(crystal);
    }
    else
    {
        mtz = MtzPtr(new MtzManager());
        mtz->setFilename(crystalFileName);
    }
    mtz->setWavelength(0);
    mtz->setImage(image);
    mtz->setSpaceGroup(spaceGroup->spg_num);
    mtz->setUnitCell(unitCell);
    
    mtz->setMatrix(newMat);
    double distance = getImage()->getDetectorDistance();
    mtz->setDetectorDistance(distance);
    
    char *hallSymbol = ccp4spg_symbol_Hall(spaceGroup);
    
    space_group _spaceGroup = space_group(hallSymbol);
    space_group_type spgType = space_group_type(_spaceGroup);
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        
        miller->incrementOverlapMask();
        miller->setMtzParent(&*mtz);
        
        int index = Reflection::indexForReflection(miller->getH(), miller->getK(), miller->getL(),
                                                   mtz->getLowGroup(), false);
        
        ReflectionPtr found = ReflectionPtr();
        mtz->findReflectionWithId(index, &found);
        
        if (found != NULL)
        {
            found->addMiller(miller);
            miller->setParent(found);
        }
        else
        {
            ReflectionPtr reflection = ReflectionPtr(new Reflection());
            reflection->setSpaceGroup(spaceGroup->spg_num);
            reflection->addMiller(miller);
            reflection->setUnitCellDouble(&unitCell[0]);
            reflection->calculateResolution(&*mtz);
            miller->setParent(reflection);
            mtz->addReflection(reflection);
            mtz->sortLastReflection();
        }
    }
    
    this->sendLog(LogLevelDetailed);
    
    double cutoff = FileParser::getKey("SIGMA_RESOLUTION_CUTOFF", SIGMA_RESOLUTION_CUTOFF);
    
    if (cutoff != 0)
        mtz->cutToResolutionWithSigma(cutoff);
    
    nearbyMillers.clear();
    lastMtz = mtz;
    
    return mtz;
}

IOMRefiner::~IOMRefiner()
{
    //    std::cout << "Deallocating IOMRefiner." << std::endl;
    
    lastMtz = MtzPtr();
    
    nearbyMillers.clear();
    vector<MillerPtr>().swap(nearbyMillers);
    
    millers.clear();
    vector<MillerPtr>().swap(millers);
    
    // FIXME: work out when this should and should not be freed
    //   if (spaceGroup != NULL)
    //       ccp4spg_free(&spaceGroup);
}

void IOMRefiner::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

std::string IOMRefiner::refinementSummaryHeader()
{
    return "Filename\tRefl num\tScore\tHoriz rot\tVert rot\tBeam rot\tStdev\tWavelength\tDistance\tUnitCellA\tUnitCellB\tUnitCellC\tAlpha\tBeta\tGamma";
}

std::string IOMRefiner::refinementSummary()
{
    std::string filename = getImage()->getFilename();
    int totalReflections = getTotalReflections();
    double lastScore = getLastScore();
    double wavelength = getWavelength();
    double distance = getDetectorDistance();
    MatrixPtr matrix = getMatrix();
    double *lengths = new double[3];
    
    for (int i = 0; i < 3; i++)
    {
        lengths[i] = 0;
    }
    
    matrix->unitCellLengths(&lengths);
    
    std::ostringstream summary;
    
    summary << filename << "\t" << totalReflections << "\t" << lastScore << "\t"
    << bestHRot << "\t" << bestKRot << "\t" << bestLRot << "\t" << lastStdev << "\t" << wavelength << "\t" << distance << "\t" << lengths[0] << "\t" << lengths[1] << "\t" << lengths[2] << "\t" << unitCell[3] << "\t" << unitCell[4] << "\t" << unitCell[5];
    
    delete [] lengths;
    
    return summary.str();
}

