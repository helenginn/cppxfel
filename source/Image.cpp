/*
 * Image.cpp
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#include "parameters.h"
#include "Image.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "FileReader.h"
#include "Logger.h"
#include "FileParser.h"
#include "Panel.h"
#include "misc.h"
#include "Shoebox.h"
#include "Spot.h"
#include "Vector.h"
#include "IndexingSolution.h"
#include "StatisticsManager.h"
#include "CSV.h"
#include "polyfit.hpp"
#include "SolventMask.h"
#include "PNGFile.h"
#include "Miller.h"
#include "SpotFinderQuick.h"
#include "SpotFinderCorrelation.h"
#include "Detector.h"

double Image::globalDetectorDistance = 0;
double Image::globalBeamX = INT_MAX;
double Image::globalBeamY = INT_MAX;

Image::Image(std::string filename, double wavelength,
             double distance)
{
    vector<double> dims = FileParser::getKey("DETECTOR_SIZE", vector<double>());
    
    xDim = 1765;
    yDim = 1765;
    
    beamX = INT_MAX;
    beamY = INT_MAX;
    
    spotsFile = "";
    
    if (dims.size())
    {
        xDim = dims[0];
        yDim = dims[1];
    }
    
    minimumSolutionNetworkCount = FileParser::getKey("MINIMUM_SOLUTION_NETWORK_COUNT", 20);
    indexingFailureCount = 0;
    data = vector<int>();
    shortData = vector<short>();
    mmPerPixel = FileParser::getKey("MM_PER_PIXEL", MM_PER_PIXEL);
    
    metrologyMoveThreshold = FileParser::getKey("METROLOGY_MOVE_THRESHOLD", 1.0);
    
    shouldMaskValue = FileParser::hasKey("IMAGE_MASKED_VALUE");
    shouldMaskUnderValue = FileParser::hasKey("IMAGE_IGNORE_UNDER_VALUE");
    
    maskedValue = 0;
    maskedUnderValue = 0;
    useShortData = false;
    
    if (shouldMaskValue)
        maskedValue = FileParser::getKey("IMAGE_MASKED_VALUE", 0);
    
    if (shouldMaskUnderValue)
        maskedUnderValue = FileParser::getKey("IMAGE_MASKED_UNDER_VALUE", 0);
    
    detectorGain = FileParser::getKey("DETECTOR_GAIN", 1.0);
    
    if (globalBeamX == INT_MAX || globalBeamY == INT_MAX)
    {
        vector<double> beam = FileParser::getKey("BEAM_CENTRE", vector<double>());
        
        if (beam.size() == 0)
        {
            beam.push_back(0);
            beam.push_back(0);
        }
        
        globalBeamX = beam[0];
        globalBeamY = beam[1];
    }
    
    _hasSeeded = false;
    
    pinPoint = true;
    int tempShoebox[7][7] =
    {
        { 1, 1, 1, 1, 1, 1, 1 },
        { 1, 0, 0, 0, 0, 0, 1 },
        { 1, 0, 2, 2, 2, 0, 1 },
        { 1, 0, 2, 2, 2, 0, 1 },
        { 1, 0, 2, 2, 2, 0, 1 },
        { 1, 0, 0, 0, 0, 0, 1 },
        { 1, 1, 1, 1, 1, 1, 1 } };
    
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            shoebox[i][j] = tempShoebox[i][j];
    
    this->filename = filename;
    this->wavelength = wavelength;
    detectorDistance = distance;
    this->fitBackgroundAsPlane = FileParser::getKey("FIT_BACKGROUND_AS_PLANE", false);
    
    loadedSpots = false;
    
    pixelCountCutoff = FileParser::getKey("PIXEL_COUNT_CUTOFF", 0);
}

std::string Image::filenameRoot()
{
    vector<std::string> components = FileReader::split(filename, '.');
    
    std::string root = "";
    
    for (int i = 0; i < components.size() - 1; i++)
    {
        root += components[i] + ".";
    }
    
    return root.substr(0, root.size() - 1);
}

void Image::setUpIOMRefiner(MatrixPtr matrix)
{
    IOMRefinerPtr indexer = IOMRefinerPtr(new IOMRefiner(shared_from_this(), matrix));
    
    if (matrix->isComplex())
        indexer->setComplexMatrix();
    
    indexers.push_back(indexer);
}

void Image::setUpIOMRefiner(MatrixPtr unitcell, MatrixPtr rotation)
{
    IOMRefinerPtr indexer = IOMRefinerPtr(new IOMRefiner(shared_from_this(), MatrixPtr(new Matrix())));
    
    indexers.push_back(indexer);
}

Image::~Image()
{
//    std::cout << "Deallocating image." << std::endl;
    
    data.clear();
    vector<int>().swap(data);
    
    overlapMask.clear();
    vector<unsigned char>().swap(overlapMask);
    
    /*
    spots.clear();
    vector<SpotPtr>().swap(spots);
    
    spotVectors.clear();
    vector<SpotVectorPtr>().swap(spotVectors);*/
}

void Image::addMask(int startX, int startY, int endX, int endY)
{
    vector<int> mask = vector<int>();
    
    mask.push_back(startX);
    mask.push_back(startY);
    mask.push_back(endX);
    mask.push_back(endY);
    
    masks.push_back(mask);
}

void Image::addSpotCover(int startX, int startY, int endX, int endY)
{
    vector<int> mask = vector<int>();
    
    mask.push_back(startX);
    mask.push_back(startY);
    mask.push_back(endX);
    mask.push_back(endY);
    
    spotCovers.push_back(mask);
}

void Image::applyMaskToImages(vector<ImagePtr> images, int startX,
                              int startY, int endX, int endY)
{
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->addMask(startX, startY, endX, endY);
    }
}

bool Image::isLoaded()
{
    
    return (data.size() > 0 || shortData.size() > 0);
}

void Image::setImageData(vector<int> newData)
{
    data.resize(newData.size());
    
    memcpy(&data[0], &newData[0], newData.size() * sizeof(int));
    
    logged << "Image with wavelength " << this->wavelength << ", distance " << getDetectorDistance() << std::endl;
    sendLog();
}

void Image::newImage()
{
    data = std::vector<int>(xDim * yDim, 0);
}

void Image::loadImage()
{
    if (isLoaded())
    {
        return;
    }
    
    std::streampos size;
    vector<char> memblock;
    
    std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
    if (file.is_open())
    {
        size = file.tellg();
        memblock = vector<char>();
        
        char c = file.get();
        
        while (file.good())
        {
            memblock.push_back(c);
            c = file.get();
        }
        
        int bitsPerPixel = FileParser::getKey("BITS_PER_PIXEL", 32);
        useShortData = (bitsPerPixel == 16);
        
        if (!useShortData)
            data.resize(memblock.size() / sizeof(int));
        else
            shortData.resize(memblock.size() / sizeof(short));
        
        overlapMask = vector<unsigned char>(memblock.size(), 0);
        
        logged << "Image size: " << memblock.size() << " for image: "
        << filename << std::endl;
        sendLog();
        
        if (!useShortData)
            memcpy(&data[0], &memblock[0], memblock.size());
        
        if (useShortData)
        {
            memcpy(&shortData[0], &memblock[0], memblock.size());
        }
    }
    else
    {
        Logger::mainLogger->addString("Unable to open file " + getFilename());
    }
    
    bool fromTiffs = FileParser::getKey("FROM_TIFFS", false);
    
    if (fromTiffs)
    {
        shortData.erase(shortData.begin(), shortData.begin() + 4);
    }
}

void Image::dropImage()
{
    data.clear();
    vector<int>().swap(data);
    
    shortData.clear();
    vector<short>().swap(shortData);
    
    overlapMask.clear();
    vector<unsigned char>().swap(overlapMask);
    
    for (int i = 0; i < IOMRefinerCount(); i++)
        getIOMRefiner(i)->dropMillers();
}

bool Image::coveredBySpot(int x, int y)
{
    for (int i = 0; i < spotCovers.size(); i++)
    {
        int startX = spotCovers[i][0];
        int startY = spotCovers[i][1];
        int endX = spotCovers[i][2];
        int endY = spotCovers[i][3];
        
        if (x >= startX && x <= endX && y >= startY && y <= endY)
            return true;
    }
    
    return false;
}

void Image::addValueAt(int x, int y, int addedValue)
{
    if (!isLoaded())
    {
        newImage();
    }

    if (x < 0 || y < 0)
        return;
    
    if (x > xDim || y > yDim)
        return;
    
    int position = y * xDim + x;
    
    data[position] = std::max(data[position], addedValue);
}

int Image::rawValueAt(int x, int y)
{
    if (!isLoaded())
    {
        loadImage();
    }
    
    if (x < 0 || y < 0)
        return 0;
    
    if (x > xDim || y > yDim)
        return 0;
    
    for (int i = 0; i < masks.size(); i++)
    {
        int startX = masks[i][0];
        int startY = masks[i][1];
        int endX = masks[i][2];
        int endY = masks[i][3];
        
        if (x >= startX && x <= endX && y >= startY && y <= endY)
            return 0;
    }
    
    int position = y * xDim + x;
    
    if (!useShortData)
    {
        if (position < 0 || position >= data.size())
            return 0;
    }
    else
    {
        if (position < 0 || position >= shortData.size())
            return 0;
    }
    
    return (useShortData ? shortData[position] : data[position]);
}

int Image::valueAt(int x, int y)
{
    double rawValue = rawValueAt(x, y);
    
    PanelPtr panel = Panel::panelForSpotCoord(std::make_pair(x, y));
    double panelGain = detectorGain;
    
    if (panel)
        panelGain *= panel->getGainScale();
    
    return rawValue * panelGain;
}

void Image::focusOnAverageMax(int *x, int *y, int tolerance1, int tolerance2, bool even)
{
    int maxValue = 0;
    int newX = *x;
    int newY = *y;
    std::string bestPixels;
    std::vector<double> values;
    
    for (int j = *y - tolerance1; j <= *y + tolerance1 + 1; j++)
    {
        for (int i = *x - tolerance1; i <= *x + tolerance1 + 1; i++)
        {
            double newValue = valueAt(i, j);
            
            if (newValue > maxValue)
            {
                newX = i;
                newY = j;
                maxValue = newValue;
            }
        }
    }
    
    *x = newX;
    *y = newY;
}

void Image::focusOnMaximum(int *x, int *y, int tolerance, double shiftX, double shiftY)
{
    int newX = *x;
    int newY = *y;
    int midValue = valueAt(newX, newY);
    
    for (int i = *x - tolerance; i <= *x + tolerance; i++)
    {
        for (int j = *y - tolerance; j <= *y + tolerance; j++)
        {
            if (valueAt(i, j) > midValue)
            {
                newX = i;
                newY = j;
                midValue = valueAt(i, j);
            }
        }
    }
    
  //  if (tolerance > 0 && accepted(newX, newY))
  //      printBox(newX, newY, tolerance);
    
    *x = newX;
    *y = newY;
}

int Image::shoeboxLength()
{
    double totalSize = 7;
    
    return totalSize;
}

Mask Image::flagAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y)
{
    Mask flag = MaskNeither;
    
    double value = (*shoebox)[x][y];
    
    if (value == -1)
        flag = MaskBackground;
    
    if (value > 0 && value <= 1)
        flag = MaskForeground;
    
    return flag;
}

double Image::weightAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y)
{
    double value = (*shoebox)[x][y];
    
    if (value == 0)
        return 0;
    
    return 1;
}

void Image::makeMaximumFromImages(std::vector<ImagePtr> images)
{
    newImage();
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < yDim; j++)
        {
            for (int k = 0; k < xDim; k++)
            {
                addValueAt(k, j, images[i]->valueAt(k, j));
            }
        }
        
        images[i]->dropImage();
    }
    
    loadedSpots = true;
    
    drawSpotsOnPNG();
    dumpImage();
}

void Image::dumpImage()
{
    std::ofstream imgStream;
    imgStream.open(getFilename().c_str(), std::ios::binary);
    
    long int size = useShortData ? shortData.size() : data.size();
    size *= useShortData ? sizeof(short) : sizeof(int);
    char *start = useShortData ? (char *)&shortData[0] : (char *)&data[0];
    
    imgStream.write(start, size);
    
    imgStream.close();
}

// @TODO
bool Image::checkShoebox(ShoeboxPtr shoebox, int x, int y)
{
    int slowSide = 0;
    int fastSide = 0;
    
    shoebox->sideLengths(&slowSide, &fastSide);
    
    int centreX = 0;
    int centreY = 0;
    
    shoebox->centre(&centreX, &centreY);
    
    int zeroCount = 0;
    int count = 0;
    
    for (int i = 0; i < slowSide; i++)
    {
        int panelPixelX = (i - centreX) + x;
        
        for (int j = 0; j < fastSide; j++)
        {
            int panelPixelY = (j - centreY) + y;
            
            if (!accepted(panelPixelX, panelPixelY))
            {
                return false;
            }
            
            count++;
            if (valueAt(panelPixelX, panelPixelY) == 0)
                zeroCount++;
        }
    }
    
    return true;
}

std::pair<double, double> Image::reciprocalCoordinatesToPixels(vec hkl, double myWavelength)
{
    if (myWavelength == 0)
    {
        myWavelength = wavelength;
    }
    double x_mm = (-hkl.h * getDetectorDistance() / (1 / myWavelength + hkl.l));
    double y_mm = (hkl.k * getDetectorDistance() / (1 / myWavelength + hkl.l));
    
    double x_coord = getBeamX() - x_mm / mmPerPixel;
    double y_coord = getBeamY() - y_mm / mmPerPixel;
    
    return std::make_pair(x_coord, y_coord);
}

vec Image::pixelsToReciprocalCoordinates(double xPix, double yPix)
{
    double mmX = xPix * getMmPerPixel();
    double mmY = yPix * getMmPerPixel();
    
    return millimetresToReciprocalCoordinates(mmX, mmY);
}

vec Image::millimetresToReciprocalCoordinates(double xmm, double ymm)
{
    double mmBeamX = getBeamX() * getMmPerPixel();
    double mmBeamY = getBeamY() * getMmPerPixel();
    
    vec crystalVec = new_vector(mmBeamX, mmBeamY, 0 - getDetectorDistance());
    vec spotVec = new_vector(xmm, ymm, 0);
    vec reciprocalCrystalVec = new_vector(0, 0, 0 - 1 / getWavelength());
    
    vec crystalToSpot = vector_between_vectors(crystalVec, spotVec);
    scale_vector_to_distance(&crystalToSpot, 1 / getWavelength());
    add_vector_to_vector(&reciprocalCrystalVec, crystalToSpot);
    
    reciprocalCrystalVec.k *= -1;
    
    return reciprocalCrystalVec;
}

void Image::printBox(int x, int y, int tolerance)
{
    std::ostringstream logged;
    
    logged << "Print box at (" << x << ", " << y << "), radius " << tolerance << std::endl;
    
    for (int i = x - tolerance; i <= x + tolerance; i++)
    {
        for (int j = y - tolerance; j <= y + tolerance; j++)
        {
            logged << valueAt(i, j) << "\t";
        }
        logged << std::endl;
    }
    
    logged << std::endl;
    
    Logger::mainLogger->addStream(&logged, LogLevelNormal);
    
}

double Image::integrateFitBackgroundPlane(int x, int y, ShoeboxPtr shoebox, float *error)
{
    int centreX = 0;
    int centreY = 0;
    
    shoebox->centre(&centreX, &centreY);
    
    int slowSide = 0;
    int fastSide = 0;
    
    shoebox->sideLengths(&slowSide, &fastSide);
    
    std::vector<double> xxs, xys, xs, yys, ys, xzs, yzs, zs, allXs, allYs, allZs;
    
    for (int i = 0; i < slowSide; i++)
    {
        int panelPixelX = (i - centreX) + x;
        
        for (int j = 0; j < fastSide; j++)
        {
            int panelPixelY = (j - centreY) + y;
            
            Mask flag = flagAtShoeboxIndex(shoebox, i, j);
            
            if (!accepted(panelPixelX, panelPixelY))
            {
                return std::nan(" ");
            }
            
            if (flag == MaskForeground || flag == MaskNeither)
                continue;
            
            double newX = panelPixelX;
            double newY = panelPixelY;
            double newZ = valueAt(panelPixelX, panelPixelY);
            
            allXs.push_back(newX);
            allYs.push_back(newY);
            allZs.push_back(newZ);
        }
    }
    
    double meanZ = weighted_mean(&allZs);
    double stdevZ = standard_deviation(&allZs);
    int rejected = 0;
    
    for (int i = 0; i < allZs.size(); i++)
    {
        double newZ = allZs[i];
        double diffZ = fabs(newZ - meanZ);
        
        if (diffZ > stdevZ * 2.2)
        {
            allXs.erase(allXs.begin() + i);
            allYs.erase(allYs.begin() + i);
            allZs.erase(allZs.begin() + i);
            i--;
            rejected++;
        }
    }
    
    logged << "Rejected background pixels: " << rejected << std::endl;
    sendLog(LogLevelDebug);
    
    for (int i = 0; i < allZs.size(); i++)
    {
        double newX = allXs[i];
        double newY = allYs[i];
        double newZ = allZs[i];
        
        xxs.push_back(newX * newX);
        yys.push_back(newY * newY);
        xys.push_back(newX * newY);
        xs.push_back(newX);
        ys.push_back(newY);
        xzs.push_back(newX * newZ);
        yzs.push_back(newY * newZ);
        zs.push_back(newZ);
    }
    
    double xxSum = sum(xxs);
    double yySum = sum(yys);
    double xySum = sum(xys);
    double xSum = sum(xs);
    double ySum = sum(ys);
    double zSum = sum(zs);
    double xzSum = sum(xzs);
    double yzSum = sum(yzs);
    
    double aveBackgroundPixel = zSum / zs.size();
    
    MatrixPtr matrix = MatrixPtr(new Matrix());
    
    matrix->components[0] = xxSum;
    matrix->components[1] = xySum;
    matrix->components[2] = xSum;
    
    matrix->components[4] = xySum;
    matrix->components[5] = yySum;
    matrix->components[6] = ySum;
    
    matrix->components[8] = xSum;
    matrix->components[9] = ySum;
    matrix->components[10] = xs.size();
    
    vec b = new_vector(xzSum, yzSum, zSum);
    
    MatrixPtr inverse = matrix->inverse3DMatrix();
    
    if (inverse->isIdentity())
    {
        return nan(" ");
    }
    
    inverse->multiplyVector(&b);
    
    // plane is now pa + qb + rc (components of b)
    
    double p = b.h;
    double q = b.k;
    double r = b.l;
    
    double backgroundInSignal = 0;
    
 //   double lowestPoint = FLT_MAX;
 //   double highestPoint = -FLT_MAX;
    std::vector<double> corners;
    /*
     corners.push_back(p * (startX) + q * (startY) + r);
     corners.push_back(p * (startX + slowSide) + q * (startY) + r);
     corners.push_back(p * (startX) + q * (startY + fastSide) + r);
     corners.push_back(p * (startX + slowSide) + q * (startY + fastSide) + r);
     
     for (int i = 0; i < corners.size(); i++)
     {
     if (corners[i] < lowestPoint)
     lowestPoint = corners[i];
     if (corners[i] > highestPoint)
     highestPoint = corners[i];
     }
     
     double topHeight = (highestPoint - lowestPoint);
     
     double bottomVolume = lowestPoint * slowSide * fastSide;
     double topVolume = 0.5 * topHeight * slowSide * fastSide;
     
     backgroundInSignal = topVolume + bottomVolume;
     */
    double foreground = 0;
    int num = 0;
    
    logged << "Foreground pixels: ";
    
    for (int i = 0; i < slowSide; i++)
    {
        int panelPixelX = (i - centreX) + x;
        
        for (int j = 0; j < fastSide; j++)
        {
            int panelPixelY = (j - centreY) + y;
            
            Mask flag = flagAtShoeboxIndex(shoebox, i, j);
            
            if (flag == MaskNeither)
                continue;
            
            else if (flag == MaskForeground)
            {
                double weight = weightAtShoeboxIndex(shoebox, i, j);
                double total = valueAt(panelPixelX, panelPixelY);
                
                logged << "F:" << total << ", ";
                
                foreground += total * weight;
                num++;
                
                double backTotal = (p * panelPixelX + q * panelPixelY + r);
                logged << "B:" << backTotal << ", ";
                
                backgroundInSignal += backTotal * weight;
            }
            
            
        }
    }
    
    
    logged << "Background: " << backgroundInSignal << std::endl;
    logged << "Foreground: " << foreground << std::endl;
    
    double signalOnly = foreground - backgroundInSignal;
    *error = sqrt(foreground);
    
    logged << std::endl;
    sendLog(LogLevelDebug);
    
    return signalOnly;
}

double Image::integrateSimpleSummation(int x, int y, ShoeboxPtr shoebox, float *error)
{
    int centreX = 0;
    int centreY = 0;
    
    shoebox->centre(&centreX, &centreY);
    
    int slowSide = 0;
    int fastSide = 0;
    
    shoebox->sideLengths(&slowSide, &fastSide);
    
    int foreground = 0;
    int foreNum = 0;
    
    int background = 0;
    int backNum = 0;
    
    for (int i = 0; i < slowSide; i++)
    {
        int panelPixelX = (i - centreX) + x;
        
        for (int j = 0; j < fastSide; j++)
        {
            int panelPixelY = (j - centreY) + y;
            
            double value = valueAt(panelPixelX, panelPixelY);
            
            Mask flag = flagAtShoeboxIndex(shoebox, i, j);
            
            if (!accepted(panelPixelX, panelPixelY))
            {
                return std::nan(" ");
            }
            
            if (flag == MaskForeground)
            {
                double weight = weightAtShoeboxIndex(shoebox, i, j);
                foreNum += weight;
                foreground += value * weight;
                
                if (value > pixelCountCutoff && pixelCountCutoff > 0)
                {
                    return std::nan(" ");
                }
            }
            else if (flag == MaskBackground)
            {
                backNum++;
                background += value;
            }
        }
    }
  
    double totalSigmaForBackground = sqrt(background);
    double averageSigmaForBackground = totalSigmaForBackground / (double)backNum;
    
    double aveBackground = (double) background / (double) backNum;
    double backgroundInForeground = aveBackground * (double) foreNum;
    
    double totalSigmaForForeground = sqrt(foreground);
    double backgroundSigmaInForeground = averageSigmaForBackground * (double)foreNum;
    *error = backgroundSigmaInForeground + totalSigmaForForeground;
    
    double intensity = (foreground - backgroundInForeground);
    
    return intensity;
}

double Image::integrateWithShoebox(int x, int y, ShoeboxPtr shoebox, float *error)
{
    if (!fitBackgroundAsPlane)
    {
        double intensity = integrateSimpleSummation(x, y, shoebox, error);
        
        return intensity;
    }
    else
    {
        return integrateFitBackgroundPlane(x, y, shoebox, error);
    }
}

double Image::intensityAt(int x, int y, ShoeboxPtr shoebox, float *error, int tolerance)
{
    int x1 = x;
    int y1 = y;
    
    if (tolerance > 0)
    {
        if (pinPoint)
            focusOnMaximum(&x1, &y1, tolerance);
        else
            focusOnAverageMax(&x1, &y1, tolerance, 1, shoebox->isEven());
    }
    
    double integral = integrateWithShoebox(x1, y1, shoebox, error);
    
    return integral;
}

bool Image::accepted(int x, int y)
{
    if (Detector::isActive())
    {
        DetectorPtr detector = Detector::getMaster()->findDetectorPanelForSpotCoord(x, y);
        
        if (!detector)
            return false;
    }
    else
    {
        PanelPtr panel = Panel::panelForSpotCoord(std::make_pair(x, y));
        
        if (!panel)
        {
            return false;
        }
    }
    
    double value = rawValueAt(x, y);
    
    if (shouldMaskValue)
    {
        if (value == maskedValue)
        {
            return false;
        }
    }
    
    if (shouldMaskUnderValue)
    {
        if (value <= maskedUnderValue)
        {
            return false;
        }
    }
    
    if (value == -100000)
        return false;
    
    return true;
}

void Image::index()
{
    if (indexers.size() == 0)
    {
        logged << "No orientation matrices, cannot index/integrate." << std::endl;
        sendLog();
        return;
    }
    
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->checkAllMillers(indexers[i]->getMaxResolution(),
                                     indexers[i]->getTestBandwidth());
    }
}

void Image::refineOrientations()
{
    if (indexers.size() == 0)
    {
        logged << "No orientation matrices, cannot index/integrate." << std::endl;
        sendLog();
        return;
    }
    
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->refineOrientationMatrix();
    }
}

vector<MtzPtr> Image::currentMtzs()
{
    vector<MtzPtr> mtzs;
    int count = 0;
    bool rejecting = !FileParser::getKey("REJECT_OVERLAPPING_REFLECTIONS", true);
    
    for (int i = 0; i < IOMRefinerCount(); i++)
    {
        MtzPtr newMtz = indexers[i]->newMtz(i);
        newMtz->removeStrongSpots(&spots, true);
        logged << "After integrating crystal " << i << " on " << filenameRoot() << ", spot count now " << spots.size() << std::endl;
        sendLog();

        std::string imgFilename = "img-" + filenameRoot() + "_" + i_to_str(i) + ".mtz";
        newMtz->writeToFile(imgFilename, true);
        newMtz->writeToDat();
        
        if (rejecting)
            count += newMtz->rejectOverlaps();
        mtzs.push_back(newMtz);
    }
    
    logged << "Generated " << mtzs.size() << " mtzs with " << count << " rejected reflections." << std::endl;
    sendLog();
    
    // fix me
    
    return mtzs;
}

int Image::throwAwayIntegratedSpots(std::vector<MtzPtr> mtzs)
{
    int thrown = 0;
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        thrown += mtzs[i]->removeStrongSpots(&spots);
    }
    
    return thrown;
}

void Image::setSpaceGroup(CCP4SPG *spg)
{
    for (int i = 0; i < IOMRefinerCount(); i++)
    {
        indexers[i]->setSpaceGroup(spg);
    }
}

void Image::setMaxResolution(double res)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setMaxResolution(res);
    }
}

void Image::setSearchSize(int searchSize)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setSearchSize(searchSize);
    }
}

void Image::setUnitCell(vector<double> dims)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setUnitCell(dims);
    }
}

void Image::setInitialStep(double step)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setInitialStep(step);
    }
}

void Image::setTestSpotSize(double spotSize)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setTestSpotSize(spotSize);
    }
}

void Image::setOrientationTolerance(double newTolerance)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setOrientationTolerance(newTolerance);
    }
}

void Image::setTestBandwidth(double bandwidth)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setTestBandwidth(bandwidth);
    }
}

void Image::incrementOverlapMask(int x, int y)
{
    int position = y * yDim + x;
    
    if (position < 0 || position >= overlapMask.size())
        return;
    
    overlapMask[position]++;
}

void Image::incrementOverlapMask(int x, int y, ShoeboxPtr shoebox)
{
    int slowSide = 0;
    int fastSide = 0;
    
    shoebox->sideLengths(&slowSide, &fastSide);
    int slowTol = (slowSide - 1) / 2;
    int fastTol = (fastSide - 1) / 2;
    
    int startX = x - slowTol;
    int startY = y - fastTol;
    for (int i = startX; i < startX + slowSide; i++)
    {
        for (int j = startY; j < startY + fastSide; j++)
        {
            incrementOverlapMask(i, j);
        }
    }
}

unsigned char Image::overlapAt(int x, int y)
{
    int position = y * yDim + x;
    
    if (position < 0 || position >= overlapMask.size())
        return 0;
    
    return overlapMask[position];
}

unsigned char Image::maximumOverlapMask(int x, int y, ShoeboxPtr shoebox)
{
    unsigned char max = 0;
    
    int slowSide = 0;
    int fastSide = 0;
    
    shoebox->sideLengths(&slowSide, &fastSide);
    int slowTol = (slowSide - 1) / 2;
    int fastTol = (fastSide - 1) / 2;
    
    int startX = x - slowTol;
    int startY = y - fastTol;
    for (int i = startX; i < startX + slowSide; i++)
    {
        for (int j = startY; j < startY + fastSide; j++)
        {
            if (overlapAt(x, y) > max)
                max = overlapAt(x, y);
        }
    }
    
    return max;
}

bool Image::checkUnitCell(double trueA, double trueB, double trueC, double tolerance)
{
    for (int i = 0; i < IOMRefinerCount(); i++)
    {
        bool isOk = true;
        double *cellDims = new double[3];
        getIOMRefiner(i)->getMatrix()->unitCellLengths(&cellDims);
        
        if (cellDims[0] < trueA - tolerance || cellDims[0] > trueA + tolerance)
            isOk = false;
        if (cellDims[1] < trueB - tolerance || cellDims[1] > trueB + tolerance)
            isOk = false;
        if (cellDims[2] < trueC - tolerance || cellDims[2] > trueC + tolerance)
            isOk = false;
        
        if (!isOk)
        {
            indexers.erase(indexers.begin() + i);
            i--;
        }
    }
    
    return IOMRefinerCount() > 0;
}

void Image::updateAllSpots()
{
    for (int i = 0; i < spots.size(); i++)
    {
        spots[i]->setUpdate();
    }
    
    for (int i = 0; i < spotVectors.size(); i++)
    {
        spotVectors[i]->setUpdate();
    }
}

void Image::addSpotIfNotMasked(SpotPtr newSpot)
{
    bool masked = SolventMask::isMasked(newSpot);
    
    if (!masked)
        spots.push_back(newSpot);
}

void Image::findSpots()
{
    int algorithm = FileParser::getKey("SPOT_FINDING_ALGORITHM", 0);
    std::vector<SpotPtr> tempSpots;
    loadImage();
    SpotFinderPtr spotFinder;
    
    if (algorithm == 1)
    {
        spotFinder = SpotFinderPtr(new SpotFinderQuick(shared_from_this()));
    }
    else if (algorithm == 0)
    {
        spotFinder = SpotFinderPtr(new SpotFinderCorrelation(shared_from_this()));
    }
    else
    {
        return;
    }
    
    tempSpots = spotFinder->findSpots();
    
    for (int i = 0; i < tempSpots.size(); i++)
    {
        addSpotIfNotMasked(tempSpots[i]);
    }
    
    loadedSpots = true;
    dropImage();

    logged << "(" << getBasename() << ") found " << spotCount() << " spots." << std::endl;
    sendLog();
    
    std::string basename = getBasename();
    Spot::writeDatFromSpots(basename + "_spots.csv", spots);
    writeSpotsList("_" + basename + "_strong.list");
    
    bool alwaysFilterSpots = FileParser::getKey("ALWAYS_FILTER_SPOTS", false);

    if (alwaysFilterSpots)
    {
        compileDistancesFromSpots(0, 0, alwaysFilterSpots);
    }
    
    return;
    
}

void Image::processSpotList()
{
    std::string spotContents;
    bool forceSpotFinding = FileParser::getKey("FORCE_SPOT_FINDING", false);
    
    if (!spotsFile.length() && !forceSpotFinding)
    {
        spotsFile = "_" + getBasename() + "_strong.list";
    }
    
    if (spotsFile == "find" || forceSpotFinding)
    {
        logged << "Finding spots using cppxfel" << std::endl;
        sendLog();
        findSpots();
    }
    else if (!FileReader::exists(spotsFile))
    {
        logged << "Cannot find spot file " << spotsFile << std::endl;
        loadImage();
        findSpots();
        sendLog();
    }
    else
    {
        bool spotsAreReciprocalCoordinates = FileParser::getKey("SPOTS_ARE_RECIPROCAL_COORDINATES", false);
        
        try
        {
            spotContents = FileReader::get_file_contents(spotsFile.c_str());
        }
        catch(int e)
        {
            logged << "Error reading spot file " << filename << std::endl;
            sendLog();
            
            return;
        }
        
        spots.clear();
        
        vector<std::string> spotLines = FileReader::split(spotContents, '\n');
        
        double x = getBeamX();
        double y = getBeamY();
        vec beamXY = new_vector(getBeamX(), getBeamY(), 1);
        vec newXY = new_vector(x, y, 1);
        vec xyVec = vector_between_vectors(beamXY, newXY);
        MatrixPtr rotateMat = MatrixPtr(new Matrix());
        rotateMat->rotate(0, 0, M_PI);
        rotateMat->multiplyVector(&xyVec);
        
        SpotPtr newSpot = SpotPtr(new Spot(shared_from_this()));
        newSpot->setXY(getBeamX() - xyVec.h, getBeamY() - xyVec.k);
        
        addSpotIfNotMasked(newSpot);
        double tooCloseDistance = 0;
        bool rejectCloseSpots = FileParser::getKey("REJECT_CLOSE_SPOTS", false);
        if (rejectCloseSpots)
        {
            tooCloseDistance = IndexingSolution::getMinDistance() * 0.7;
        }
        
        for (int i = 0; i < spotLines.size(); i++)
        {
            std::string line = spotLines[i];
            vector<std::string> components = FileReader::split(line, '\t');
            
            if (components.size() < 2)
                continue;
            
            if (spotsAreReciprocalCoordinates && components.size() < 3)
                continue;
            
            double x = 0; double y = 0;
            
            if (spotsAreReciprocalCoordinates)
            {
                double k = atof(components[0].c_str());
                double h = atof(components[1].c_str());
                double l = atof(components[2].c_str());
                
                vec hkl = new_vector(h, k, l);
                
                Coord coord = reciprocalCoordinatesToPixels(hkl);
                
                x = coord.first;
                y = coord.second;
                
                logged << "Coord:\t" << h << "\t" << k << "\t" << l << "\t" << x << "\t" << y << std::endl;
                sendLog();
            }
            else
            {
                x = atof(components[0].c_str());
                y = atof(components[1].c_str());
            }
            
            vec beamXY = new_vector(getBeamX(), getBeamY(), 1);
            vec newXY = new_vector(x, y, 1);
            vec xyVec = vector_between_vectors(beamXY, newXY);
            MatrixPtr rotateMat = MatrixPtr(new Matrix());
            rotateMat->rotate(0, 0, M_PI);
            rotateMat->multiplyVector(&xyVec);
            
            SpotPtr newSpot = SpotPtr(new Spot(shared_from_this()));
            newSpot->setXY(getBeamX() - xyVec.h, getBeamY() - xyVec.k);
            bool add = true;
            
            vec myVec = newSpot->estimatedVector();
            
            for (int j = 0; j < spots.size(); j++)
            {
                SpotPtr testSpot = spots[j];
                
                if (rejectCloseSpots && tooCloseDistance > 0)
                {
                    vec testVec = testSpot->estimatedVector();
                    vec copyVec = copy_vector(myVec);
                    take_vector_away_from_vector(testVec, &copyVec);
                    
                    double distance = length_of_vector(copyVec);
                    if (distance < tooCloseDistance)
                    {
                        add = false;
                    }
                }
            }
            
            for (int j = 0; j < spots.size(); j++)
            {
                if (newSpot->isSameAs(spots[j]))
                {
                    add = false;
                }
            }
            
            logged << "SPOT\t" << i << "\t" << newSpot->getRawXY().first << "\t"
            << newSpot->getRawXY().second
            << "\t" << newSpot->getX() << "\t" << newSpot->getY() << std::endl;
            sendLog(LogLevelDebug);
            
            if (!Detector::isActive() && Panel::panelForSpot(&*newSpot) == PanelPtr())
            {
                add = false;
            }
            
            if (add)
            {
                addSpotIfNotMasked(newSpot);
            }
        }
        
        logged << "Loaded " << spots.size() << " spots from list " << spotsFile << std::endl;
        sendLog(LogLevelNormal);
    }
    
    double fractionToExclude = FileParser::getKey("EXCLUDE_WEAKEST_SPOT_FRACTION", 0.0);
    
    if (fractionToExclude > 0)
    {
        excludeWeakestSpots(fractionToExclude);
    }
    
    bool recentreSpots = FileParser::getKey("FORCE_RECENTRE_SPOTS", false);
    
    if (recentreSpots)
    {
        for (int i = 0; i < spotCount(); i++)
        {
            spots[i]->recentreInWindow();
        }
        
        writeSpotsList("_" + getBasename() + "_strong.list");
    }
    
    loadedSpots = true;
}


void Image::rotatedSpotPositions(MatrixPtr rotationMatrix, std::vector<vec> *spotPositions, std::vector<std::string> *spotElements)
{
    double maxRes = 1.0;
    
    for (int i = 0; i < spots.size(); i++)
    {
        vec spotPos = spots[i]->estimatedVector();
        bool goodSpot = (spots[i]->successfulLineCount() > 0);
        
        if (length_of_vector(spotPos) > 1 / maxRes)
            continue;
        
        rotationMatrix->multiplyVector(&spotPos);
        spotPositions->push_back(spotPos);
        
        std::string element = goodSpot ? "N" : "O";
        spotElements->push_back(element);
    }
}

bool Image::acceptableSpotCount()
{
    int maxSpots = FileParser::getKey("REJECT_OVER_SPOT_COUNT", 4000);
    int minSpots = FileParser::getKey("REJECT_UNDER_SPOT_COUNT", 0);

    return (spotCount() > minSpots) && (spotCount() < maxSpots);
    
}

void Image::compileDistancesFromSpots(double maxReciprocalDistance, double tooCloseDistance, bool filter)
{
    if (!loadedSpots)
    {
        processSpotList();
    }
    
    bool rejectCloseSpots = FileParser::getKey("REJECT_CLOSE_SPOTS", false);
    double minResolution = FileParser::getKey("INDEXING_MIN_RESOLUTION", 0.0);
    double maxResolution = FileParser::getKey("MAX_INTEGRATED_RESOLUTION", 0.0);
    
    if (rejectCloseSpots && tooCloseDistance == 0)
    {
        tooCloseDistance = IndexingSolution::getMinDistance() * 0.7;
    }
    
    int maxSpots = FileParser::getKey("REJECT_OVER_SPOT_COUNT", 4000);
    int minSpots = FileParser::getKey("REJECT_UNDER_SPOT_COUNT", 0);
    
    if (maxSpots > 0)
    {
        if (spotCount() > maxSpots)
        {
            logged << "(" << getFilename() << ") Aborting image due to too many spots." << std::endl;
            sendLog();
            spotVectors.clear();
            std::vector<SpotVectorPtr>().swap(spotVectors);
            return;
        }
    }
    
    if (spotCount() < minSpots)
    {
        logged << "(" << getFilename() << ") Aborting image due to too few spots." << std::endl;
        sendLog();
        
        if (IOMRefinerCount() == 0)
        {
            sendLog();
        }
        spotVectors.clear();
        std::vector<SpotVectorPtr>().swap(spotVectors);
        return;
    }
    
    
    if (maxReciprocalDistance == 0)
    {
        maxReciprocalDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    }
    
    spotVectors.clear();
    std::vector<SpotVectorPtr>().swap(spotVectors);
    
    if (spots.size() == 0)
        return;
    
    for (int i = 0; i < spots.size() - 1; i++)
    {
        if (spots[i]->isRejected() && rejectCloseSpots)
        {
            i++;
            continue;
        }
        
        vec spotPos1 = spots[i]->estimatedVector();
        
        if (minResolution != 0)
        {
            if (length_of_vector(spotPos1) < 1 / minResolution)
                continue;
        }
        
        if (maxResolution > 0)
        {
            if (length_of_vector(spotPos1) > 1 / maxResolution)
                continue;
        }
        
        for (int j = i + 1; j < spots.size(); j++)
        {
            if (spots[j]->isRejected())
            {
                j++;
                continue;
            }
            
            vec spotPos2 = spots[j]->estimatedVector();
            
            bool close = within_vicinity(spotPos1, spotPos2, maxReciprocalDistance);
            
            bool tooClose = within_vicinity(spotPos1, spotPos2, tooCloseDistance);
            
            if (tooClose && rejectCloseSpots)
            {
                spots[j]->setRejected();
                spots[i]->setRejected();
                i++;
            }
            
            if (close)
            {
                SpotVectorPtr newVec = SpotVectorPtr(new SpotVector(spots[i], spots[j]));
                
                double distance = newVec->distance();
                
                if (distance == 0)
                    continue;
                
                if (distance > maxReciprocalDistance)
                    continue;
                
                logged << "vec\t" << spots[i]->getX() << "\t" << spots[i]->getY() << "\t" << spots[j]->getX() << "\t" << spots[j]->getY() << "\t0\t0\t0\t" << distance << std::endl;
                sendLog(LogLevelDebug);
                
                spotVectors.push_back(newVec);
            }
        }
        
        sendLog(LogLevelDetailed);
    }
    
    
    if (filter)
    {
        filterSpotVectors();
    }
    
    bool scramble = FileParser::getKey("SCRAMBLE_SPOT_VECTORS", true);
    
    if (scramble)
    {
        //   std::sort(spotVectors.begin(), spotVectors.end(), SpotVector::isGreaterThan);
        
        std::random_shuffle(spotVectors.begin(), spotVectors.end());
    }
    
    logged << "(" << filename << ") " << spotCount() << " spots produced " << spotVectorCount() << " spot vectors." << std::endl;
    sendLog();
}

bool solutionBetterThanSolution(IndexingSolutionPtr one, IndexingSolutionPtr two)
{
    return (one->spotVectorCount() > two->spotVectorCount());
}

void Image::filterSpotVectors()
{
    int spotsPerLattice = FileParser::getKey("SPOTS_PER_LATTICE", 100);
    
    std::vector<double> scoresOnly;
    std::map<SpotVectorPtr, double> spotVectorMap;
    double reciprocalTolerance = FileParser::getKey("INITIAL_RLP_SIZE", 0.0015);
    
    int totalSpots = spotCount();
    double expectedLatticesFraction = (double)totalSpots / (double)spotsPerLattice;
    int goodHits = round(expectedLatticesFraction);
    int maxVectors = 12000;
    
    double goodFraction = proportion(goodHits);
    
    logged << "From " << totalSpots << " spots there is an estimated " << expectedLatticesFraction << " lattices" << std::endl;
    logged << "Fraction of spot vectors which should be good: " << goodFraction << std::endl;
    
    sendLog();
    
    for (int i = 0; i < spotVectorCount() && i < maxVectors; i++)
    {
        SpotVectorPtr spotVec1 = spotVector(i);
        double score = 0;
        
        for (int j = 0; j < spotVectorCount() && j < maxVectors; j++)
        {
            if (j == i)
                continue;
            
            SpotVectorPtr spotVec2 = spotVector(j);
            
            if (spotVec1->isCloseToSpotVector(spotVec2, reciprocalTolerance))
            {
                double interDistance = spotVec1->similarityToSpotVector(spotVec2);
                
                if (interDistance < reciprocalTolerance)
                    score += reciprocalTolerance - interDistance;
            }
        }
        
        spotVectorMap[spotVec1] = score;
        scoresOnly.push_back(score);
    }
    
    std::sort(scoresOnly.begin(), scoresOnly.end(), std::greater<int>());
    
    int vectorsToKeep = int((double)goodFraction * (double)scoresOnly.size());
    
    logged << "Keeping vectors: " << vectorsToKeep << " out of total: " << scoresOnly.size() << std::endl;
    sendLog();
    
    if (vectorsToKeep == scoresOnly.size())
        vectorsToKeep--;
    
    if (scoresOnly.size() == 0)
        return;
    
    if (vectorsToKeep == 0)
        return;
    
    double threshold = scoresOnly[vectorsToKeep];
    int deleted = 0;
    
    logged << "Threshold: " << threshold << std::endl;
    
    for (int i = 0; i < spotVectorCount(); i++)
    {
        if (spotVectorMap.count(spotVector(i)) == 0 || spotVectorMap[spotVector(i)] <= threshold)
        {
            spotVectors.erase(spotVectors.begin() + i);
            i--;
            deleted++;
        }
    }
    
    logged << "Deleted: " << deleted << std::endl;
    
    sendLog();
}

bool Image::checkIndexingSolutionDuplicates(MatrixPtr newSolution, bool excludeLast)
{
    for (int i = 0; i < IOMRefinerCount() - excludeLast; i++)
    {
        MatrixPtr oldSolution = getIOMRefiner(i)->getMatrix();
        
        bool similar = IndexingSolution::matrixSimilarToMatrix(newSolution, oldSolution, true);
        
        if (similar)
            return true;
    }
    
    return false;
}

IndexingSolutionStatus Image::tryIndexingSolution(IndexingSolutionPtr solutionPtr, bool modify, int *spotsRemoved, IOMRefinerPtr *aRefiner)
{
    logged << "(" << filename << ") Trying solution from " << solutionPtr->spotVectorCount() << " vectors." << std::endl;
    sendLog(LogLevelNormal);
    
    if (modify)
    {
        solutionPtr->setModMatrix(Reflection::getCustomAmbiguity());
    }
    else
    {
        solutionPtr->setModMatrix(MatrixPtr());
    }
    
    MatrixPtr solutionMatrix = solutionPtr->createSolution();
    bool similar = checkIndexingSolutionDuplicates(solutionMatrix);
    
    if (similar)
    {
        logged << "Indexing solution too similar to previous solution. Continuing..." << std::endl;
        sendLog(LogLevelNormal);
        solutionPtr->removeSpotVectors(&spotVectors);
        
        return IndexingSolutionTrialDuplicate;
    }
    
    bool refineOrientations = FileParser::getKey("REFINE_ORIENTATIONS", true);
    int minimumSpotsExplained = FileParser::getKey("MINIMUM_SPOTS_EXPLAINED", 20);
    
    logged << solutionPtr->printNetwork();
    logged << solutionPtr->getNetworkPDB();
    
    sendLog(LogLevelDetailed);
    
    setUpIOMRefiner(solutionMatrix);
    int lastRefiner = IOMRefinerCount() - 1;
    IOMRefinerPtr refiner = getIOMRefiner(lastRefiner);
    *aRefiner = refiner;

    if (refineOrientations)
    {
        refiner->refineOrientationMatrix();
        bool similar = checkIndexingSolutionDuplicates(refiner->getMatrix(), true);
        
        if (similar)
        {
            removeRefiner(lastRefiner);
            
            logged << "Indexing solution too similar to previous solution after refinement. Continuing..." << std::endl;
            sendLog(LogLevelNormal);
            
            return IndexingSolutionTrialDuplicate;
        }
    }
    else
    {
        refiner->calculateOnce();
    }
    
    bool successfulImage = refiner->isGoodSolution();
    
    MtzPtr mtz = refiner->newMtz(lastRefiner, true);
    *spotsRemoved = mtz->removeStrongSpots(&spots, false);
    int minSpotsForImage = std::min(minimumSpotsExplained, (int)(spotCount() * 0.8));
    
    if (*spotsRemoved < minSpotsForImage)
    {
        logged << "(" << getFilename() << ") However, does not explain enough spots (" << *spotsRemoved << " vs  " << minSpotsForImage << ")" << std::endl;
        sendLog();
        successfulImage = false;
    }
    else
    {
        logged << "(" << getFilename() << ") Enough spots are explained (" << *spotsRemoved << " vs  " << minimumSpotsExplained << ")" << std::endl;
        sendLog();
        
    }
    
    removeRefiner(lastRefiner);
    
    return successfulImage ? IndexingSolutionTrialSuccess : IndexingSolutionTrialFailure;
}

IndexingSolutionStatus Image::tryIndexingSolution(IndexingSolutionPtr solutionPtr)
{
    bool acceptAllSolutions = FileParser::getKey("ACCEPT_ALL_SOLUTIONS", false);
    IOMRefinerPtr normRefiner, modRefiner;
    int spotsRemoved = 0;
    int modSpotsRemoved = 0;
    IndexingSolutionStatus status = tryIndexingSolution(solutionPtr, false, &spotsRemoved, &normRefiner);
    IndexingSolutionStatus modStatus = IndexingSolutionTrialFailure;
    
    if (Reflection::getCustomAmbiguity())
    {
        modStatus = tryIndexingSolution(solutionPtr, true, &modSpotsRemoved, &modRefiner);
    }
    
    bool modBetter = false;

    if ((status != IndexingSolutionTrialSuccess && modStatus != IndexingSolutionTrialSuccess) && !acceptAllSolutions)
    {
        badSolutions.push_back(solutionPtr);
        indexingFailureCount++;
        return status;
    }
    else if (status != IndexingSolutionTrialSuccess && modStatus == IndexingSolutionTrialSuccess)
    {
        modBetter = true;
    }
    else if (status == IndexingSolutionTrialSuccess && modStatus == IndexingSolutionTrialSuccess)
    {
        if (modSpotsRemoved > spotsRemoved)
        {
            modBetter = true;
        }
    }
    
    if (modBetter)
    {
        solutionPtr->setModMatrix(Reflection::getCustomAmbiguity());
    }
    else
    {
        solutionPtr->setModMatrix(MatrixPtr());
    }
    
    logged << "(" << getFilename() << ") " << (modBetter ? "custom ambiguity" : "expected ambiguity") << " is better." << std::endl;
    sendLog();

    IOMRefinerPtr myRefiner = modBetter ? modRefiner : normRefiner;

    int lastRefiner = IOMRefinerCount();
    addIOMRefiner(myRefiner);
    MtzPtr mtz = myRefiner->newMtz(lastRefiner, true);
    
    logged << "Successful crystal for " << getFilename() << std::endl;
    Logger::log(logged); logged.str("");
    goodSolutions.push_back(solutionPtr);
    int spotCountBefore = (int)spots.size();
    myRefiner->showHistogram(false);
    
    mtz->removeStrongSpots(&spots);
    compileDistancesFromSpots();
    IndexingSolution::calculateSimilarStandardVectorsForImageVectors(spotVectors);
    
    int spotCountAfter = (int)spots.size();
    
    logged << "Removed spots; from " << spotCountBefore << " to " << spotCountAfter << "." << std::endl;
    
    Logger::log(logged); logged.str("");
    
    addMtz(mtz);
    std::string imgFilename = "img-" + mtz->getFilename();
    mtz->writeToFile(imgFilename, true);
    mtz->writeToDat("img-");
    
    return IndexingSolutionTrialSuccess;
}

IndexingSolutionStatus Image::extendIndexingSolution(IndexingSolutionPtr solutionPtr, std::vector<SpotVectorPtr> existingVectors, int *failures, int added)
{
    int newFailures = 0;
    
    if (failures == NULL)
    {
        failures = &newFailures;
    }
    
    
    std::vector<SpotVectorPtr> newVectors = existingVectors;
    
    if (!solutionPtr)
    {
        logged << "Solution pointer not pointing" << std::endl;
        sendLog();
        return IndexingSolutionBranchFailure;
    }
    int newlyAdded = 1;
    int trials = 0;
    int trialLimit = FileParser::getKey("NETWORK_TRIAL_LIMIT", 3);
    
    while (newlyAdded > 0 && added < 100 && trials < trialLimit)
    {
        IndexingSolutionPtr copyPtr = solutionPtr->copy();
        
        newlyAdded = copyPtr->extendFromSpotVectors(&newVectors, 1);
        
        if (newlyAdded > 0)
        {
            trials++;
            
            if (Logger::getPriorityLevel() >= LogLevelDetailed)
            {
                logged << "Starting new branch with " << added + newlyAdded << " additions (trial " << trials << ")." << std::endl;
                sendLog(LogLevelDetailed);
            }
            
            IndexingSolutionStatus success = extendIndexingSolution(copyPtr, newVectors, failures, added + newlyAdded);
            
            if (success == IndexingSolutionBranchFailure)
            {
                if (!biggestFailedSolution || copyPtr->spotVectorCount() > biggestFailedSolution->spotVectorCount())
                {
                    biggestFailedSolution = copyPtr;
                    biggestFailedSolutionVectors = newVectors;
                }
            }

            if (success != IndexingSolutionBranchFailure)
            {
                return success;
            }
            else
            {
                if (trials >= trialLimit)
                {
                    (*failures)++;
                    if (Logger::getPriorityLevel() >= LogLevelDetailed)
                    {
                        logged << "Given up this branch, too many failures." << std::endl;
                        sendLog(LogLevelDetailed);
                    }
                    
                    return success;
                }
            }
        }
        
        if (*failures > 3)
        {
            logged << "Giving up on this thread, too many failures" << std::endl;
            sendLog(LogLevelDetailed);
            return IndexingSolutionBranchFailure;
        }
    }
    
    if (added >= minimumSolutionNetworkCount)
    {
        IndexingSolutionStatus success = tryIndexingSolution(solutionPtr);
        
        return success;
    }
    
    if (added < minimumSolutionNetworkCount)
    {
        if (Logger::getPriorityLevel() >= LogLevelDetailed)
        {
            logged << "Didn't go anywhere..." << std::endl;
            sendLog(LogLevelDetailed);
        }
    }
    
    existingVectors.clear();
    std::vector<SpotVectorPtr>().swap(existingVectors);
    
    return IndexingSolutionBranchFailure;
}

std::vector<double> Image::anglesBetweenVectorDistances(double distance1, double distance2, double tolerance)
{
    std::vector<SpotVectorPtr> firstVectors, secondVectors;
    
    for (int i = 0; i < spotVectors.size(); i++)
    {
        double vecDistance = spotVectors[i]->distance();
        double diff1 = fabs(vecDistance - distance1);
        double diff2 = fabs(vecDistance - distance2);
        
        if (diff1 < 1 / tolerance)
        {
            logged << "Image " << getFilename() << ", adding vector " << i << " " << spotVectors[i]->description() << " to group 1" << std::endl;
            sendLog();
            firstVectors.push_back(spotVectors[i]);
        }
        
        if (diff2 < 1 / tolerance)
        {
            logged << "Image " << getFilename() << ", adding vector " << i << " " << spotVectors[i]->description() << " to group 2" << std::endl;
            sendLog();
            secondVectors.push_back(spotVectors[i]);
        }
    }
    
 //   if (firstVectors.size() <= 1 && secondVectors.size() <= 1)
 //       return std::vector<double>();
    
    if (firstVectors.size() && secondVectors.size())
    {
        logged << "N: Image " << getFilename() << " has " << firstVectors.size() << " and " << secondVectors.size() << "  vector distances on same image." << std::endl;
        sendLog();
    }
    else return std::vector<double>();
    
    std::vector<double> angles;
    
    for (int j = 0; j < firstVectors.size(); j++)
    {
        for (int k = 0; k < secondVectors.size(); k++)
        {
            if (firstVectors[j] == secondVectors[k])
                continue;
            
            double angle = firstVectors[j]->angleWithVector(secondVectors[k]);
            logged << "Adding angle between " << firstVectors[j]->description() << " and " << secondVectors[k]->description() << " " << angle * 180 / M_PI << std::endl;
            sendLog();
            
            angles.push_back(angle);
            angles.push_back(M_PI - angle);
        }
    }
    
    return angles;
}

IndexingSolutionStatus Image::testSeedSolution(IndexingSolutionPtr newSolution, std::vector<SpotVectorPtr> &prunedVectors, int *successes)
{
    bool similar = checkIndexingSolutionDuplicates(newSolution->createSolution(), false);
    
    if (similar)
    {
        logged << "Solution too similar to another. Continuing..." << std::endl;
        sendLog(LogLevelDetailed);
        return IndexingSolutionTrialDuplicate;
    }
    
    logged << "Starting a new solution..." << std::endl;
    sendLog(LogLevelDetailed);
    
    IndexingSolutionStatus success = extendIndexingSolution(newSolution, prunedVectors);
    
    if (success == IndexingSolutionTrialSuccess)
    {
        logged << "(" << getFilename() << ") indexing solution trial success." << std::endl;
        (*successes)++;
    }
    else if (success == IndexingSolutionTrialFailure)
    {
        logged << "(" << getFilename() << ") indexing solution trial failure." << std::endl;
    }
    else if (success == IndexingSolutionTrialDuplicate)
    {
        logged << "(" << getFilename() << ") indexing solution trial duplicate." << std::endl;
    }
    
    return success;
}

void Image::findIndexingSolutions()
{
    if (!loadedSpots)
    {
        processSpotList();
    }

    std::vector<IndexingSolutionPtr> solutions;
    
    int maxSearch = FileParser::getKey("MAX_SEARCH_NUMBER_MATCHES", 1000);
    
    if (IOMRefinerCount() > 0)
    {
        logged << "Existing solution spot removal (image " << getFilename() << "):" << std::endl;
        
        for (int i = 0; i < IOMRefinerCount(); i++)
        {
            MtzPtr mtz = getIOMRefiner(i)->newMtz(i);
            int spotCountBefore = (int)spots.size();
            mtz->removeStrongSpots(&spots);
            int spotCountAfter = (int)spots.size();
            int spotDiff = spotCountBefore - spotCountAfter;
            
            logged << "Removed " << spotDiff << " spots from existing solution leaving " << spotCountAfter << " spots." << std::endl;
        }

        sendLog();
    }
    
    bool alwaysFilterSpots = FileParser::getKey("ALWAYS_FILTER_SPOTS", false);
    compileDistancesFromSpots(0, 0, alwaysFilterSpots);
    if (spotVectors.size() == 0)
        return;
    
    sendLog();
    
    bool continuing = true;
    bool allowBiggestSolution = !FileParser::getKey("ACCEPT_ALL_SOLUTIONS", false);
    int successes = 0;
    int maxSuccesses = FileParser::getKey("SOLUTION_ATTEMPTS", 1);
    int maxLattices = FileParser::getKey("MAX_LATTICES_PER_IMAGE", 1);
    
    if (maxLattices < maxSuccesses)
        maxLattices = maxSuccesses;
    
    if (IOMRefinerCount() >= maxLattices)
        return;
    
    int indexingTimeLimit = FileParser::getKey("INDEXING_TIME_LIMIT", 1200);
    
    IndexingSolution::calculateSimilarStandardVectorsForImageVectors(spotVectors);
    
    if (spotVectors.size() == 0)
    {
        logged << "No vectors - giving up." << std::endl;
        sendLog();
        return;
    }
    
    time_t startcputime;
    time(&startcputime);
    
    bool lastWasSuccessful = true;
    
    while (lastWasSuccessful)
    {
        for (int i = 0; i < spotVectors.size() - 1 && i < maxSearch && continuing && indexingFailureCount < 10; i++)
        {
            SpotVectorPtr spotVector1 = spotVectors[i];
            
            for (int j = i + 1; j < spotVectors.size() && continuing && indexingFailureCount < 10; j++)
            {
                SpotVectorPtr spotVector2 = spotVectors[j];
                
                std::vector<IndexingSolutionPtr> moreSolutions = IndexingSolution::startingSolutionsForVectors(spotVector1, spotVector2);
                
                for (int k = 0; k < moreSolutions.size(); k++)
                {
                    IndexingSolutionStatus status = testSeedSolution(moreSolutions[k], spotVectors, &successes);
                    
                    if (status == IndexingSolutionTrialSuccess || status == IndexingSolutionTrialDuplicate)
                    {
                        if (spotVectors.size() == 0)
                        {
                            continuing = false;
                            break;
                        }
                        
                        logged << "(" << getFilename() << ") now on " << spotVectors.size() << " vectors." << std::endl;
                    }
                    
                    if (successes >= maxSuccesses || IOMRefinerCount() >= maxLattices)
                    {
                        continuing = false;
                    }
                }
                
                moreSolutions.clear();
                std::vector<IndexingSolutionPtr>().swap(moreSolutions);
            }
            
            time_t middlecputime;
            time(&middlecputime);
            
            clock_t difference = middlecputime - startcputime;
            double seconds = difference;
            
            if (seconds > indexingTimeLimit)
            {
                logged << "N: Time limit reached on image " << filename << " on " << IOMRefinerCount() << " crystals and " << spotCount() << " remaining spots." << std::endl;
                sendLog();
                
                return;
            }
        }
        
        if (continuing && allowBiggestSolution)
        {
            if (!biggestFailedSolution)
            {
                lastWasSuccessful = false;
                continue;
            }
            
            IndexingSolutionStatus status = tryIndexingSolution(biggestFailedSolution);
            
            if (status != IndexingSolutionTrialSuccess)
            {
                lastWasSuccessful = false;
            }
            
            if (status == IndexingSolutionTrialSuccess)
            {
                spotVectors = biggestFailedSolutionVectors;
                biggestFailedSolution = IndexingSolutionPtr();
            }
        }
        else
        {
            break;
        }
    }
    
    logged << "N: Finished image " << filename << " on " << IOMRefinerCount() << " crystals and " << spotCount() << " remaining spots." << std::endl;
    sendLog();
    
    dropImage();
}

std::vector<MtzPtr> Image::getLastMtzs()
{
    std::vector<MtzPtr> mtzs;
    
    for (int i = 0; i < IOMRefinerCount(); i++)
    {
        MtzPtr lastMtz = getIOMRefiner(i)->getLastMtz();
        
        if (lastMtz)
        {
            mtzs.push_back(lastMtz);
        }
    }
    
    return mtzs;
}

void Image::writeSpotsList(std::string spotFile)
{
    std::string spotDat;
    
    if (spotFile == "")
    {
        std::string tag = "remaining";
        std::string basename = getBasename();
        
        spotFile = "_" + basename + "_" + tag + "_spots.list";
        spotDat = basename + "_spots.csv";
    }
    
    std::ofstream spotList;
    spotList.open(spotFile);
    
    for (int i = 0; i < spotCount(); i++)
    {
        spotList << spot(i)->spotLine();
    }
    
    spotList.close();
    
    Spot::writeDatFromSpots(spotDat, spots);
}

void Image::reset()
{
    clearIOMRefiners();
    failedRefiners.clear();
    minimumSolutionNetworkCount = FileParser::getKey("MINIMUM_SOLUTION_NETWORK_COUNT", 20);
    
    processSpotList();
}

double Image::resolutionAtPixel(double x, double y)
{
    vec reciprocal = pixelsToReciprocalCoordinates(x, y);
    return length_of_vector(reciprocal);
}

void Image::integrateSpots()
{
    CSV integratedCSV = CSV(3, "x", "y", "intensity");
    
    for (int i = 0; i < spotCount(); i++)
    {
        double intensity = spot(i)->integrate();
        
        integratedCSV.addEntry(0, spot(i)->getX(), spot(i)->getY(), intensity);
    }
    
    std::string csvName = "spot-int-" + getBasename() + ".csv";
    
    integratedCSV.writeToFile(csvName);
    
    logged << "Written to file: " << csvName << std::endl;
    sendLog();
}

bool compareSpotIntensities(SpotPtr her, SpotPtr him)
{
    return (her->getIntensity() < him->getIntensity());
}

void Image::excludeWeakestSpots(double fraction)
{
    std::vector<SpotPtr> sortedSpots;
    
    for (int i = 0; i < spotCount(); i++)
    {
        spot(i)->integrate();
        sortedSpots.push_back(spot(i));
    }
    
    std::sort(sortedSpots.begin(), sortedSpots.end(), compareSpotIntensities);
    
    int numExclude = fraction * sortedSpots.size();
    
    if (sortedSpots.size() == 0)
    {
        return;
    }
    double minIntensity = sortedSpots[0]->getIntensity();
    double maxIntensity = sortedSpots[numExclude]->getIntensity();
    
    /* want to preserve original order of spots I think */
    for (int i = 0; i < numExclude; i++)
    {
        for (int j = 0; j < spotCount(); j++)
        {
            if (spot(j) == sortedSpots[i])
            {
                spots.erase(spots.begin() + j);
                j--;
            }
        }
    }
    
    std::cout << "Removed " << numExclude << " spots between intensities "
    << minIntensity << " and " << maxIntensity << std::endl;
    
}

void Image::radialAverage()
{
    double maxIntegratedResolution = FileParser::getKey("MAX_INTEGRATED_RESOLUTION", 3.5);
    double minIntegratedResolution = FileParser::getKey("MIN_INTEGRATED_RESOLUTION", 100);
    double binCount = FileParser::getKey("BINS", 100);
    
    std::vector<double> bins;
    StatisticsManager::generateResolutionBins(minIntegratedResolution, maxIntegratedResolution, binCount, &bins);
    
    std::vector<double> intensities = std::vector<double>(bins.size(), 0);
    std::vector<double> counts = std::vector<double>(bins.size(), 0);
    std::vector<double> subtracted = std::vector<double>(bins.size(), 0);
    
    for (int i = 0; i < xDim; i++)
    {
        for (int j = 0; j < yDim; j++)
        {
            double resolution = 1 / resolutionAtPixel(i, j);
            double value = valueAt(i, j);
            
            for (int k = 0; k < bins.size() - 1; k++)
            {
                if ((resolution < bins[k]) && (resolution > bins[k + 1]))
                {
                    intensities[k] += value;
                    counts[k] += 1;
                }
            }
        }
    }
    
    for (int i = 0; i < bins.size(); i++)
    {
        intensities[i] /= counts[i];
    }
    
    double buffer = binCount * 0.1;
    double bufferRatio = 3;
    logged << "Buffer size set to: " << buffer << std::endl;
    sendLog();
    
    for (int i = buffer; i < intensities.size() - buffer; i++)
    {
        std::vector<double> intensityPortion, binsPortion;
        
        for (int j = i - buffer; j < i + buffer; j++)
        {
            if (j < i + (buffer / bufferRatio) && j > i - (buffer / bufferRatio))
            {
                continue;
            }
            
            binsPortion.push_back(bins[j]);
            intensityPortion.push_back(intensities[j]);
        }
        
        for (int j = 0; j < binsPortion.size(); j++)
        {
            logged << binsPortion[j] << "\t" << intensityPortion[j] << std::endl;
        }
        
        sendLog();
        
        std::vector<double> coeffs = polyfit(binsPortion, intensityPortion, 2);
        
        double backgroundSum = 0;
        
        for (int j = i - buffer / bufferRatio; j < i + buffer / bufferRatio; j++)
        {
            double x = bins[j];
            backgroundSum += coeffs[2] * pow(x, 2) + coeffs[1] * x + coeffs[0];
        }
        
        backgroundSum /= buffer * 2 / bufferRatio;
        
        subtracted[i] = intensities[i] - backgroundSum;
    }
    
    CSV csv = CSV(2, "1 over d squared", "Average intensity");
    
    for (int i = 0; i < bins.size(); i++)
    {
        double bFacRes = 1 / bins[i];
        csv.addEntry(0, bFacRes, subtracted[i]);
    }
    
    csv.writeToFile("sub-" + getFilename() + ".csv");
    
    CSV csv2 = CSV(2, "1 over d squared", "Average intensity");
    
    for (int i = 0; i < bins.size(); i++)
    {
        double bFacRes = 1 / bins[i];
        csv2.addEntry(0, bFacRes, intensities[i]);
    }
    
    csv2.writeToFile("rad-" + getFilename() + ".csv");
    
    logged << "Written csv: sub-" << getFilename() << ".csv" << std::endl;
    sendLog();
    
    dropImage();
}

double Image::standardDeviationOfPixels()
{
    std::vector<double> values;
    
    for (int i = 0; i < xDim; i++)
    {
        for (int j = 0; j < yDim; j++)
        {
            int value = valueAt(i, j);
            
            if (value < 0)
                value = 0;
            
            values.push_back(value);
        }
    }
    
    return standard_deviation(&values);
}

void Image::clusterCountWithSpotNumber(int spotNum)
{
    double minDistance = IndexingSolution::getMinDistance();
    
    if (!loadedSpots)
    {
        processSpotList();
    }
    
    if (spots.size() < 3)
    {
        return;
    }
    
    std::vector<int> results = std::vector<int>(spotNum + 1, 0);
    std::vector<SpotPtr> allSpotsIncluded;
    
    logged << getFilename() << "\t";
    
    for (int i = 0; i < spots.size() - 2; i++)
    {
        SpotPtr spot1 = spots[i];
        for (int j = i + 1; j < spots.size() - 1; j++)
        {
            SpotPtr spot2 = spots[j];
            
            SpotVectorPtr vector = SpotVectorPtr(new SpotVector(spot1, spot2));
            std::vector<SpotPtr> includedSpots;
            
            if (vector->distance() < minDistance)
            {
                includedSpots.push_back(spot1);
                includedSpots.push_back(spot2);
                allSpotsIncluded.push_back(spot1);
                allSpotsIncluded.push_back(spot2);
                
                for (int k = j + 1; k < spots.size(); k++)
                {
                    SpotPtr spot3 = spots[k];
                    bool ok = true;
                    
                    std::vector<SpotPtr>::iterator possibleSpot = std::find(allSpotsIncluded.begin(), allSpotsIncluded.end(), spot3);
                    
                    if (possibleSpot != allSpotsIncluded.end())
                    {
                        continue;
                    }
                    
                    for (int l = 0; l < includedSpots.size(); l++)
                    {
                        SpotVectorPtr checkVec = SpotVectorPtr(new SpotVector(includedSpots[l], spot3));
                        
                        if (checkVec->distance() > minDistance)
                        {
                            ok = false;
                        }
                    }
                    
                    if (ok)
                    {
                        includedSpots.push_back(spot3);
                        allSpotsIncluded.push_back(spot3);
                    }
                    
                    if (includedSpots.size() == spotNum)
                    {
                        break;
                    }
                }
                
                size_t clusterSize = includedSpots.size();
                
                results[clusterSize]++;
                includedSpots.clear();
            }
        }
    }

    logged << spots.size() << "\t";
    
    for (int i = 2; i < spotNum + 1; i++)
    {
        logged << results[i] << "\t";
    }
    
    logged << std::endl;
    
    sendLog();
}

void Image::drawSpotsOnPNG(std::string filename)
{
    if (!loadedSpots)
    {
        processSpotList();
    }
    
    double stdev = standardDeviationOfPixels();
    
    if (!filename.length())
        filename = getBasename() + ".png";
    int height = FileParser::getKey("PNG_HEIGHT", 2400);
    PNGFilePtr file = PNGFilePtr(new PNGFile(filename, height, height));
    writePNG(file);
    
    for (int i = 0; i < spotCount(); i++)
    {
        Coord coord = spot(i)->getXY();
        
        file->drawCircleAroundPixel(coord.first, coord.second, 10, 1, 0, 0, 0);
    }
    
    logged << "Written file " << filename << " (stdev of pixel intensities: " << stdev << ")" << std::endl;
    sendLog();
    
    file->writeImageOutput();
}

void Image::drawMillersOnPNG(PNGFilePtr file, MtzPtr myMtz, char red, char green, char blue)
{
    bool addShoebox = FileParser::getKey("PNG_SHOEBOX", false);
    
    for (int i = 0; i < myMtz->reflectionCount(); i++)
    {
        for (int j = 0; j < myMtz->reflection(i)->millerCount(); j++)
        {
            MillerPtr myMiller = myMtz->reflection(i)->miller(j);
            
            Coord pngCoord;
            
            if (Detector::isActive())
            {
                vec intersection = myMiller->getShiftedRay();
                
                pngCoord = std::make_pair(intersection.h, intersection.k);
            }
            else
            {
                PanelPtr panel = Panel::panelForMiller(&*myMiller);
                
                if (!panel)
                    continue;
                
                double correctedX = myMiller->getCorrectedX();
                double correctedY = myMiller->getCorrectedY();
                
                Coord corrected = std::pair<double, double>(correctedX, correctedY);
                pngCoord = panel->spotToMillerCoord(corrected);
            }
            
            bool strong = myMiller->reachesThreshold();
            file->drawCircleAroundPixel(pngCoord.first, pngCoord.second, 14, (strong ? 2 : 1), red, green, blue, (strong ? 4 : 1));
            
            if (addShoebox)
            {
                file->drawShoeboxAroundPixel(pngCoord.first, pngCoord.second, myMiller->getShoebox());
            }
        }
    }
}

void Image::drawCrystalsOnPNG(int crystalNum)
{
    int count = mtzCount();
    std::string crystNumString = (crystalNum > 0 ? i_to_str(crystalNum) : "all" + i_to_str(count));
    
    std::string filename = getBasename() + "_" + crystNumString + ".png";
    int height = FileParser::getKey("PNG_HEIGHT", 2400);
    PNGFilePtr file = PNGFilePtr(new PNGFile(filename, height, height));
    writePNG(file);
    
    if (crystalNum >= 0)
    {
        MtzPtr myMtz = mtz(crystalNum);
        drawMillersOnPNG(file, myMtz, 0, 0, 255);
    }
    else
    {
        double hueJump = 360 / (double)count;
        double hue = 0;
        
        for (int i = 0; i < count; i++)
        {
            double brightness = 0.8;
            double saturation = 1.0;
            png_byte red, green, blue;
            PNGFile::HSB_to_RGB(hue, saturation, brightness, &red, &blue, &green);
            
            MtzPtr myMtz = mtz(i);
            drawMillersOnPNG(file, myMtz, red, green, blue);
            
            hue += hueJump;
        }
    }
    
    logged << "Written file " << filename << std::endl;
    sendLog();
    
    file->writeImageOutput();
}

void Image::writePNG(PNGFilePtr file)
{
    file->setCentre(getBeamX(), getBeamY());
    
    double defaultThreshold = standardDeviationOfPixels() * 5; // calculate
    
    double threshold = FileParser::getKey("PNG_THRESHOLD", defaultThreshold);
    
    PanelPtr lastPanel = PanelPtr();
    DetectorPtr lastDetector = DetectorPtr();
    double minZ = 0;
    double maxZ = 0;
    
    if (Detector::isActive())
    {
        Detector::getMaster()->zLimits(&minZ, &maxZ);
        double nudge = std::max(0.1 * (maxZ - minZ), 0.1);
        
        minZ -= nudge;
        maxZ += nudge;
    }
    
    for (int i = 0; i < xDim; i++)
    {
        for (int j = 0; j < yDim; j++)
        {
            double value = valueAt(i, j);
            if (value < 0)
                value = 0;
            
            unsigned char pixelValue = std::min(value, threshold) * 255 / threshold;
            pixelValue = 255 - pixelValue;
            float brightness = 1 - std::min(value, threshold) / threshold;
            
            Coord coord = std::make_pair(i, j);
            
            PanelPtr panel = lastPanel;
        //    PanelPtr mask = lastMask;
            
            DetectorPtr detector = lastDetector;
            
            if (Detector::isActive())
            {
                if (!detector || !detector->spotCoordWithinBounds(i, j))
                {
                    detector = Detector::getMaster()->findDetectorPanelForSpotCoord(i, j);
                    
                    if (detector)
                    {
                        lastDetector = detector;
                    }
                }
                
                if (detector)
                {
                    vec arranged;
                    detector->spotCoordToAbsoluteVec(i, j, &arranged);
                    double proportion_distance = (arranged.l - minZ) / (maxZ - minZ) * 120 + 120;
                    png_byte red, green, blue;
                    
                    PNGFile::HSB_to_RGB(proportion_distance, 1.0, brightness, &red, &green, &blue);
                    
                    file->setPixelColourRelative(arranged.h, arranged.k, red, green, blue);
                }
                
                continue;
            }
            
            if (panel)
            {
                Coord shifted = lastPanel->spotToMillerCoord(coord);
                
                if (!panel->isMillerCoordInPanel(shifted) || Panel::spotCoordFallsInMask(shifted))
                {
                    panel = Panel::panelForSpotCoord(coord);
                }
            }
            else
            {
                panel = Panel::panelForSpotCoord(coord);
            }
            
            if (panel)
            {
                lastPanel = panel;
          //      lastMask = mask;
                
                Coord realPosition = panel->realCoordsForImageCoord(coord);
                file->setPixelColourRelative(realPosition.first, realPosition.second, pixelValue, pixelValue, pixelValue);
            }
        }
        
    //    std::cout << std::endl;
    }
}
