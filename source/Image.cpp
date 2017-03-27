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


std::vector<DetectorPtr> Image::perPixelDetectors;
vector<signed char> Image::generalMask;
ImagePtr Image::_imageMask;
std::mutex Image::setupMutex;
bool Image::interpolate = false;

Image::Image(std::string filename, double wavelength,
             double distance)
{
    vector<double> dims = FileParser::getKey("DETECTOR_SIZE", vector<double>());
    
    xDim = 1765;
    yDim = 1765;
    
    beamX = 0;
    beamY = 0;
    
    int mss = FileParser::getKey("METROLOGY_SEARCH_SIZE", -1);
    
    if (mss == 0)
    {
        interpolate = true;
    }
    
    
    spotsFile = "";
    highScore = 0;
    fake = false;
    _isMask = false;
    averageZ = 0;
    
    if (dims.size())
    {
        xDim = dims[0];
        yDim = dims[1];
        
        logged << "Setting xDim/yDim to " << xDim << ", " << yDim << std::endl;
        sendLog(LogLevelDebug);
    }
    
    minimumSolutionNetworkCount = FileParser::getKey("MINIMUM_SOLUTION_NETWORK_COUNT", 20);
    indexingFailureCount = 0;
    data = vector<int>();
    shortData = vector<short>();
    mmPerPixel = FileParser::getKey("MM_PER_PIXEL", MM_PER_PIXEL);
    
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
    
	vector<double> beam = FileParser::getKey("BEAM_CENTRE", vector<double>());

	if (beam.size() == 2)
	{
		beamX = beam[0];
		beamY = beam[1];
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
    
    this->setFilename(filename);
    this->wavelength = wavelength;
    detectorDistance = distance;
    this->fitBackgroundAsPlane = FileParser::getKey("FIT_BACKGROUND_AS_PLANE", false);
    
    loadedSpots = false;
    
    pixelCountCutoff = FileParser::getKey("PIXEL_COUNT_CUTOFF", 0);
}

void Image::setUpCrystal(MatrixPtr matrix)
{
	MtzPtr mtz = MtzPtr(new MtzManager());
	std::string numStr = i_to_str(mtzCount());
	mtz->setImage(shared_from_this());

	mtz->setFilename("img-" + getBasename() + "_" + numStr + ".mtz");
	mtz->setMatrix(matrix);

	mtzs.push_back(mtz);
}

Image::~Image()
{
    data.clear();
    vector<int>().swap(data);
    
    overlapMask.clear();
    vector<signed char>().swap(overlapMask);

    spots.clear();
    vector<SpotPtr>().swap(spots);
    
    spotVectors.clear();
    vector<SpotVectorPtr>().swap(spotVectors);
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
}

void Image::newImage()
{
    int totalPixels = xDim * yDim;
    fake = true;
    
    data = std::vector<int>(totalPixels, 0);
    overlapMask = vector<signed char>(totalPixels, 0);
    
    checkAndSetupLookupTable();

}

void Image::loadBadPixels()
{
    std::string pixelList = FileParser::getKey("BAD_PIXELS", std::string(""));
    
    if (!pixelList.length())
    {
        return;
    }
    else if (!FileReader::exists(pixelList))
    {
        logged << "Pixel list file " << pixelList << " appears to be missing." << std::endl;
        sendLogAndExit();
    }
    
    std::string pixelContents = FileReader::get_file_contents(pixelList.c_str());
    std::vector<std::string> pixels = FileReader::split(pixelContents, '\n');
    int count = 0;
    
    for (int i = 0; i < pixels.size(); i++)
    {
        std::vector<std::string> coords = FileReader::split(pixels[i], ' ');

        if (coords.size() < 2)
        {
            coords = FileReader::split(pixels[i], '\t');

            if (coords.size() < 2)
            {
                continue;
            }
        }
        
        int xCoord = atoi(coords[0].c_str());
        int yCoord = atoi(coords[1].c_str());
        
        int pos = yCoord * xDim + xCoord;
        generalMask[pos] = 0;
        
        count++;
    }
    
    logged << "Loaded " << count << " bad pixels from " << pixelList << "." << std::endl;
    sendLog();
}

void Image::checkAndSetupLookupTable()
{
    if (!perPixelDetectors.size() && (!_isMask))
    {
        setupMutex.lock();
        
        if (!perPixelDetectors.size())
        {
            int totalSize = xDim * yDim;

            logged << "Setting up detector lookup table for size " << xDim << " " << yDim << std::endl;
            sendLog();

			generalMask = vector<signed char>(totalSize, -1);
            perPixelDetectors = vector<DetectorPtr>(totalSize, DetectorPtr());
            
            for (int y = 0; y < yDim; y++)
            {
                for (int x = 0; x < xDim; x++)
                {
                    int pos = y * xDim + x;
                    
                    DetectorPtr det = Detector::getMaster()->findDetectorPanelForSpotCoord(x, y);
                    
                    bool isDet = (det != DetectorPtr());
                    
                    perPixelDetectors[pos] = det;
                    generalMask[pos] = isDet;
                    
                }
            }
        }
        
        loadBadPixels();
        
        setupMutex.unlock();
    }
}

void Image::loadImage()
{
    if (isLoaded())
    {
        return;
    }
    
    bool asFloat = FileParser::getKey("HDF5_AS_FLOAT", false);
    
    std::streampos size;
    vector<char> memblock;
    
    std::ifstream file(getFilename().c_str(), std::ios::in | std::ios::binary);
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
		overlapMask = vector<signed char>(memblock.size(), 0);

        if (!useShortData)
            data.resize(memblock.size() / sizeof(int));
        else
            shortData.resize(memblock.size() / sizeof(short));

        logged << "Image size: " << memblock.size() << " for image: "
        << getFilename() << std::endl;
        sendLog();
        
        if (asFloat)
        {
            data.reserve(size / sizeof(int));
            
            for (int i = 0; i < memblock.size(); i += sizeof(float))
            {
                float value = (float)memblock[i];
                data.push_back(value);
            }
        }
        else if (!useShortData)
        {
            memcpy(&data[0], &memblock[0], memblock.size());
        }
        else if (useShortData)
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
    
    checkAndSetupLookupTable();
}

void Image::dropImage()
{
    data.clear();
    vector<int>().swap(data);
    
    shortData.clear();
    vector<short>().swap(shortData);
    
    overlapMask.clear();
    vector<signed char>().swap(overlapMask);
    
    for (int i = 0; i < mtzCount(); i++)
        mtz(i)->dropMillers();
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
    loadImage();

	if (!isLoaded())
	{
		return 0;
	}

    if (x < 0 || y < 0)
        return 0;
    
    if (x > xDim || y > yDim)
        return 0;
    
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

double Image::interpolateAt(double x, double y, double *total)
{
    int leftX = x;
    int leftY = y;
    int rightX = leftX + 1;
    int rightY = leftY + 1;

    double leftPropX = (double)rightX - x;
    double rightPropX = x - (double)leftX;
    
    double leftPropY = (double)rightY - y;
    double rightPropY = y - (double)leftY;
    
    double topLeftSq = leftPropX * leftPropY;
    double topRightSq = rightPropX * leftPropY;
    double bottomLeftSq = leftPropX * rightPropY;
    double bottomRightSq = rightPropX * rightPropY;
    
    *total = 1;
    
    if (!accepted(leftX, leftY))
    {
        *total -= topLeftSq;
        topLeftSq = 0;
    }
    
    if (!accepted(leftX, rightY))
    {
        *total -= bottomLeftSq;
        bottomLeftSq = 0;
    }

    if (!accepted(rightX, leftY))
    {
        *total -= topRightSq;
        topRightSq = 0;
    }
    
    if (!accepted(rightX, rightY))
    {
        *total -= bottomRightSq;
        bottomRightSq = 0;
    }
    
    double topLeftValue = valueAt(leftX, leftY);
    double topRightValue = valueAt(rightX, leftY);
    double bottomLeftValue = valueAt(leftX, rightY);
    double bottomRightValue = valueAt(rightX, rightY);
    
    double weighted = topLeftSq * topLeftValue;
    weighted += topRightSq * topRightValue;
    weighted += bottomRightSq * bottomRightValue;
    weighted += bottomLeftSq * bottomLeftValue;
    
    return weighted / *total;
}

int Image::valueAt(int x, int y)
{
    loadImage();
    ImagePtr mask = getImageMask();
    
    int pos = y * xDim + x;
    
    if (pos < 0 || pos > xDim * yDim)
    {
        return false;
    }
    
    signed char maskValue = generalMask[pos];
    
    if (maskValue == 0)
    {
        return 0;
    }
    
    if (mask && (&*mask != this))
    {
        if (mask->valueAt(x, y) > 0)
        {
            generalMask[pos] = 0;
            return 0;
        }
    }
    
    double rawValue = rawValueAt(x, y);
    
    DetectorPtr det = perPixelDetectors[pos];
    
    if (!det)
    {
        return 0;
    }
    
    double panelGain = det->getGain();
    
    return rawValue * panelGain;
}

void Image::focusOnAverageMax(double *x, double *y, int tolerance1, int tolerance2, bool even)
{
    int maxValue = 0;
    int newX = *x;
    int newY = *y;
    std::string bestPixels;
    std::vector<double> values;
    
    for (int j = *y - tolerance1; j <= *y + tolerance1; j++)
    {
        for (int i = *x - tolerance1; i <= *x + tolerance1; i++)
        {
            if (!accepted(i, j))
                continue;
            
            double newValue = valueAt(i, j);
            
            if (newValue > maxValue)
            {
                newX = i;
                newY = j;
                maxValue = newValue;
            }
        }
    }
    
    *x = (double)newX + 0.5;
    *y = (double)newY + 0.5;
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

void Image::makeMaximumFromImages(std::vector<ImagePtr> images, bool listResults)
{
    if (images.size() == 0)
    {
        return;
    }
    
    images[0]->loadImage();
    
    xDim = images[0]->getXDim();
    yDim = images[0]->getYDim();
    newImage();
    maxes = std::vector<ImagePtr>(xDim * yDim, ImagePtr());
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < yDim; j++)
        {
            for (int k = 0; k < xDim; k++)
            {
                int pos = j * xDim + k;
                double rawValue = images[i]->rawValueAt(k, j);
                double myValue = rawValueAt(k, j);
                
                if (rawValue > myValue)
                {
                    if (!images[i]->fake)
                    {
                        maxes[pos] = images[i];
                    }
                    else
                    {
                        maxes[pos] = images[i]->maxes[pos];
                    }
                    
                    addValueAt(k, j, images[i]->rawValueAt(k, j));
                }
            }
        }
        
        images[i]->dropImage();
    }
    
    findSpots();
    std::map<ImagePtr, int> imageSpotMap;
    
    for (int i = 0; i < spotCount(); i++)
    {
        int x = spot(i)->getRawXY().first;
        int y = spot(i)->getRawXY().second;
        int pos = y * xDim + x;
        
        if (maxes[pos])
        {
            if (listResults)
            {
                logged << "Spot for image " << maxes[pos]->getFilename() << " (" << x << ", " << y << ")" << std::endl;
            }
            
            if (!imageSpotMap.count(maxes[pos]))
            {
                imageSpotMap[maxes[pos]] = 0;
            }
            imageSpotMap[maxes[pos]]++;
        }
    }
    
    std::string datname = FileReader::addOutputDirectory("all-strong.dat");
    int count = 0;
    
    if (listResults)
    {
        std::ofstream strongest;
        strongest.open(datname);
        strongest << "version 3.0" << std::endl;
        
        for (std::map<ImagePtr, int>::iterator it = imageSpotMap.begin(); it != imageSpotMap.end(); it++)
        {
            logged << it->first->getBasename() << "\t" << it->second << std::endl;
        
            if (it->second > 10)
            {
                strongest << "image " << it->first->getBasename() << std::endl;
                count++;
            }
        }
        
        strongest.close();
        
        logged << "Written strongest " << count << " files (10+ peaks) to " << datname << std::endl;
        logged << "You can limit yourself to just these images by editing your input file:" << std::endl;
        logged << std::endl << "ORIENTATION_MATRIX_LIST " << datname << std::endl;
        sendLog();
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
    
    std::vector<double> corners;
    
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

double Image::integrateSimpleSummation(double x, double y, ShoeboxPtr shoebox, float *error)
{
    int centreX = 0;
    int centreY = 0;
    
    shoebox->centre(&centreX, &centreY);
    
    int slowSide = 0;
    int fastSide = 0;
    
    shoebox->sideLengths(&slowSide, &fastSide);
    
    double foreground = 0;
    int foreNum = 0;
    
    double background = 0;
    int backNum = 0;
    
    int rejects = 0;
    
    for (int i = 0; i < slowSide; i++)
    {
        double shoeX = i;
        double panelPixelX = (shoeX - centreX) + x;
        
        if (shoebox->isEven())
        {
            panelPixelX += 0.5;
        }
        
        for (int j = 0; j < fastSide; j++)
        {
            double shoeY = j;
            double panelPixelY = (shoeY - centreY) + y;
            
            if (shoebox->isEven())
            {
                panelPixelY += 0.5;
            }
            
            Mask flag = flagAtShoeboxIndex(shoebox, i, j);
            
            if (!accepted(panelPixelX, panelPixelY))
            {
                rejects++;
                
                if (flag == MaskForeground || rejects > 4)
                {
                    return std::nan(" ");
                }
            }
            
            double value = 0;
            double pixelSize = 1;
            
            if (!interpolate)
            {
                value = valueAt(panelPixelX, panelPixelY);
            }
            else
            {
                value = interpolateAt(panelPixelX, panelPixelY, &pixelSize);
                
                if (pixelSize <= 0.1)
                {
                    rejects++;
                    continue;
                }
            }
            
            if (flag == MaskForeground)
            {
                double weight = weightAtShoeboxIndex(shoebox, i, j);
                weight *= pixelSize;
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
    
	double aveBackground = 0;

	if (backNum > 0)
	{
		aveBackground = (double) background / (double) backNum;
	}
    double backgroundInForeground = aveBackground * (double) foreNum;
    
    double totalSigmaForForeground = sqrt(foreground);
    double backgroundSigmaInForeground = averageSigmaForBackground * (double)foreNum;
    *error = backgroundSigmaInForeground + totalSigmaForForeground;
    
    double intensity = (foreground - backgroundInForeground);
    
    return intensity;
}

double Image::integrateWithShoebox(double x, double y, ShoeboxPtr shoebox, float *error)
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

double Image::intensityAt(double x, double y, ShoeboxPtr shoebox, float *error, int tolerance)
{
    double x1 = x;
    double y1 = y;
    
    if (error )
    
    if (tolerance > 0)
    {
        focusOnAverageMax(&x1, &y1, tolerance, 1, shoebox->isEven());
    }
    
    double integral = integrateWithShoebox(x1, y1, shoebox, error);
    
    return integral;
}

bool Image::accepted(int x, int y)
{
    int pos = xDim * y + x;
    
    if (pos < 0 || pos > xDim * yDim)
    {
        return false;
    }
    
    if (pos < generalMask.size() && generalMask[pos] == 0)
    {
        return false;
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
    
    return true;
}

void Image::refineOrientations()
{
    if (mtzCount() == 0)
    {
        logged << "No crystals, moving on..." << std::endl;
        sendLog();
        return;
    }
    
    for (int i = 0; i < mtzCount(); i++)
    {
        mtzs[i]->refineOrientationMatrix();
    }
}

vector<MtzPtr> Image::currentMtzs()
{
    bool rejecting = FileParser::getKey("REJECT_OVERLAPPING_REFLECTIONS", true);
	int count = 0;

	if (!overlapMask.size())
	{
		overlapMask = vector<signed char>(xDim * yDim, 0);
	}

    for (int i = 0; i < mtzCount(); i++)
    {
        mtz(i)->removeStrongSpots(&spots, true);
        mtz(i)->writeToFile("", true);

        if (rejecting)
            count += mtz(i)->rejectOverlaps();
    }

	logged << "After integrating for " << getBasename() << ", spot count now " << spots.size() << std::endl;
	sendLog();

    logged << "Generated " << mtzs.size() << " mtzs with " << count << " rejected reflections." << std::endl;
    sendLog();

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

void Image::incrementOverlapMask(int x, int y)
{
    int position = y * xDim + x;
    
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
    for (int i = 0; i < slowSide; i++)
    {
        for (int j = 0; j < fastSide; j++)
        {
            double x = j + startX;
            double y = i + startY;
            
            if (flagAtShoeboxIndex(shoebox, i, j) != MaskForeground)
            {
                continue;
            }

            incrementOverlapMask(x, y);
        }
    }
}

signed char Image::maskValueAt(int x, int y)
{
    int position = y * xDim + x;
    int maxSize = xDim * yDim;
    
    if (position < 0 || position >= maxSize)
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
    for (int i = 0; i < slowSide; i++)
    {
        for (int j = 0; j < fastSide; j++)
        {
            double x = j + startX;
            double y = i + startY;
            
            if (flagAtShoeboxIndex(shoebox, i, j))
            
            if (maskValueAt(x, y) > max)
            {
                max = maskValueAt(x, y);
            }
        }
    }
    
    return max;
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

void Image::fakeSpots()
{
    spots.clear();
    spotVectors.clear();
    
    if (mtzCount() == 0)
    {
        int num = rand() % 2 + 1;
        
        for (int i = 0; i < num; i++)
        {
            MatrixPtr mat = Matrix::randomOrientationMatrix();
            MtzPtr anMtz = MtzPtr(new MtzManager());
			anMtz->setImage(shared_from_this());
			anMtz->setMatrix(mat);
            
			addMtz(anMtz);
        }
    }
    
    for (int i = 0; i < mtzCount(); i++)
    {
        mtz(i)->fakeSpots();
    }
    
    logged << "Generated " << spotCount() << " fake spots." << std::endl;
    sendLog();
    
    setSpotsFile("_" + getBasename() + "_fake.list");
    writeSpotsList();
    
    dropImage();
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

void Image::weedOutCloseSpots()
{
    bool rejectCloseSpots = FileParser::getKey("REJECT_CLOSE_SPOTS", false);
    
    if (!rejectCloseSpots)
    {
        return;
    }
    
    std::vector<size_t> indicesToDelete;
    double tooCloseDistance = IndexingSolution::getMinDistance() * 0.8;

    for (int i = 0; i < spotCount() - 1; i++)
    {
        bool keepGoingI = true;
        
        SpotPtr firstSpot = spot(i);

        for (int j = i + 1; j < spotCount(); j++)
        {
            bool keepGoingJ = true;
            
            for (int check = 0; check < indicesToDelete.size(); check++)
            {
                if (indicesToDelete[check] == i)
                {
                    keepGoingI = false;
                }
                
                if (indicesToDelete[check] == j)
                {
                    keepGoingJ = false;
                }
            }
            
            if (!keepGoingI)
            {
                break;
            }
            
            if (!keepGoingJ)
            {
                continue;
            }
            
            SpotPtr secondSpot = spot(j);

            vec firstVec = firstSpot->estimatedVector();
            vec secondVec = secondSpot->estimatedVector();
            take_vector_away_from_vector(firstVec, &secondVec);
            
            double distance = length_of_vector(secondVec);
            if (distance < tooCloseDistance)
            {
                indicesToDelete.push_back(i);
                indicesToDelete.push_back(j);
            }
        }
    }
    
    logged << "Rejected " << indicesToDelete.size() << " spots for being too close." << std::endl;
    sendLog();
    
    for (int i = (int)indicesToDelete.size() - 1; i >= 0; i--)
    {
        spots.erase(spots.begin() + i);
    }
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
        try
        {
            spotContents = FileReader::get_file_contents(spotsFile.c_str());
        }
        catch(int e)
        {
            logged << "Error reading spot file for " << getFilename() << std::endl;
            sendLog();
            
            return;
        }
        
        spots.clear();
        
        vector<std::string> spotLines = FileReader::split(spotContents, '\n');
        
        double minIndexing = FileParser::getKey("MIN_INTEGRATED_RESOLUTION", 0.0);
        
        if (minIndexing > 0)
        {
            minIndexing = 1 / minIndexing;
        }
        
        for (int i = 0; i < spotLines.size(); i++)
        {
            std::string line = spotLines[i];
            vector<std::string> components = FileReader::split(line, '\t');
            bool fake = false;
            
            if (components.size() < 2)
                continue;

            if (components.size() >= 3)
            {
                std::string flags = components[2];
                
                vector<std::string> flagComponents = FileReader::split(flags, ',');
                
                for (int j = 0; j < flagComponents.size(); j++)
                {
                    if (flagComponents[j] == "fake")
                    {
                        fake = true;
                    }
                }
            }
            
            double x = 0; double y = 0;
            x = atof(components[0].c_str());
            y = atof(components[1].c_str());
            
            vec beamXY = new_vector(getBeamX(), getBeamY(), 1);
            vec newXY = new_vector(x, y, 1);
            vec xyVec = vector_between_vectors(beamXY, newXY);
            MatrixPtr rotateMat = MatrixPtr(new Matrix());
            rotateMat->rotate(0, 0, M_PI);
            rotateMat->multiplyVector(&xyVec);
            
            SpotPtr newSpot = SpotPtr(new Spot(shared_from_this()));
            newSpot->setXY(getBeamX() - xyVec.h, getBeamY() - xyVec.k);
            newSpot->setFake(fake);
            bool add = true;
            
            vec myVec = newSpot->estimatedVector();
            
            vec aCopy = myVec;
            double dStar = length_of_vector(aCopy);
            
            if (dStar < minIndexing)
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
    
    SpotPtr newSpot = SpotPtr(new Spot(shared_from_this()));
    newSpot->setToBeamCentre();
    
    addSpotIfNotMasked(newSpot);
    
    double fractionToExclude = FileParser::getKey("EXCLUDE_WEAKEST_SPOT_FRACTION", 0.0);
    
    weedOutCloseSpots();
    
    if (fractionToExclude > 0)
    {
        excludeWeakestSpots(fractionToExclude);
    }
    
    loadedSpots = true;
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
        
        if (mtzCount() == 0)
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
            vec spotPos2 = spots[j]->estimatedVector();
            
            bool close = within_vicinity(spotPos1, spotPos2, maxReciprocalDistance);
            
            if (close)
            {
                SpotVectorPtr newVec = SpotVectorPtr(new SpotVector(spots[i], spots[j]));
                
                double distance = newVec->distance();
                
                if (distance == 0)
                    continue;
                
                if (distance > maxReciprocalDistance)
                    continue;
                
                spotVectors.push_back(newVec);
            }
        }
        
        sendLog(LogLevelDetailed);
    }
    
    
    if (filter)
    {
        filterSpotVectors();
    }
    
    std::sort(spotVectors.begin(), spotVectors.end(), SpotVector::isGreaterThan);
    
    
    logged << "(" << getFilename() << ") " << spotCount() << " spots produced " << spotVectorCount() << " spot vectors." << std::endl;
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
    for (int i = 0; i < mtzCount() - excludeLast; i++)
    {
        MatrixPtr oldSolution = mtz(i)->getMatrix();
        
        bool similar = IndexingSolution::matrixSimilarToMatrix(newSolution, oldSolution, true);
        
        if (similar)
            return true;
    }
    
    for (int i = 0; i < badSolutions.size(); i++)
    {
        MatrixPtr badSol = badSolutions[i]->createSolution();
        
        bool similar = IndexingSolution::matrixSimilarToMatrix(newSolution, badSol, true);
        
        if (similar)
            return true;
    }
    
    return false;
}

IndexingSolutionStatus Image::tryIndexingSolution(IndexingSolutionPtr solutionPtr, bool modify, int *spotsRemoved, MtzPtr *anMtz)
{
    logged << "(" << getFilename() << ") Trying solution from " << solutionPtr->spotVectorCount() << " vectors." << std::endl;
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
        badSolutions.push_back(solutionPtr);
        solutionPtr->removeSpotVectors(&spotVectors);
        
        return IndexingSolutionTrialDuplicate;
    }
    
    int minimumSpotsExplained = FileParser::getKey("MINIMUM_SPOTS_EXPLAINED", 20);
    
    logged << solutionPtr->printNetwork();
    logged << solutionPtr->getNetworkPDB();
    
    sendLog(LogLevelDetailed);
    
    setUpCrystal(solutionMatrix);
    int lastMtz = mtzCount() - 1;
    MtzPtr myMtz = mtz(lastMtz);
    *anMtz = myMtz;

	myMtz->refineOrientationMatrix(true);
	similar = checkIndexingSolutionDuplicates(myMtz->getMatrix(), true);

	if (similar)
	{
		removeCrystal(lastMtz);

		badSolutions.push_back(solutionPtr);
		logged << "Indexing solution too similar to previous solution after refinement. Continuing..." << std::endl;
		sendLog(LogLevelNormal);

		return IndexingSolutionTrialDuplicate;
	}

    bool successfulImage = myMtz->isGoodSolution();

    int minSpotsForImage = std::min(minimumSpotsExplained, (int)(spotCount() * 0.5));

	std::vector<SpotPtr> someSpots = spots;
	*spotsRemoved = myMtz->removeStrongSpots(&someSpots);

    if (*spotsRemoved < minSpotsForImage)
    {
        logged << "(" << getFilename() << ") However, does not explain enough spots (" << *spotsRemoved << " vs  " << minSpotsForImage << ")" << std::endl;
        sendLog();
        successfulImage = false;
    }
    else
    {
        logged << "(" << getFilename() << ") Enough spots are explained (" << *spotsRemoved << " vs  " << minSpotsForImage << ")" << std::endl;
        sendLog();
		spots = someSpots;
		myMtz->getWavelengthHistogram(NULL, NULL, false);
        
    }
    
    removeCrystal(lastMtz);
    
    return successfulImage ? IndexingSolutionTrialSuccess : IndexingSolutionTrialFailure;
}

IndexingSolutionStatus Image::tryIndexingSolution(IndexingSolutionPtr solutionPtr)
{
    bool acceptAllSolutions = FileParser::getKey("ACCEPT_ALL_SOLUTIONS", false);
    MtzPtr normRefiner, modRefiner;
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
    
    if (status == IndexingSolutionTrialDuplicate || modStatus == IndexingSolutionTrialDuplicate)
    {
        return IndexingSolutionTrialDuplicate;
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

    MtzPtr myMtz = modBetter ? modRefiner : normRefiner;
    myMtz->setIndexingSolution(solutionPtr);
    
    addMtz(myMtz);

    logged << "Successful crystal for " << getFilename() << std::endl;
    Logger::log(logged); logged.str("");
    goodSolutions.push_back(solutionPtr);

    myMtz->removeStrongSpots(&spots);
    compileDistancesFromSpots();
    IndexingSolution::calculateSimilarStandardVectorsForImageVectors(spotVectors);

	myMtz->writeToFile("", true);
    
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
    
    int maxSearch = FileParser::getKey("MAX_SEARCH_NUMBER_MATCHES", 2000);
    
    if (mtzCount() > 0)
    {
        logged << "Existing solution spot removal (image " << getFilename() << "):" << std::endl;
        
        for (int i = 0; i < mtzCount(); i++)
        {
			MtzPtr anMtz = mtz(i);
			anMtz->integrateChosenMillers();
			int spotCountBefore = (int)spots.size();
            anMtz->removeStrongSpots(&spots);
            int spotCountAfter = (int)spots.size();
            int spotDiff = spotCountBefore - spotCountAfter;
			anMtz->getWavelengthHistogram(NULL, NULL, false);

            logged << "Removed " << spotDiff << " spots from existing solution leaving " << spotCountAfter << " spots." << std::endl;
        }

        sendLog();
    }
    
    bool reachedTime = false;
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
    
    if (mtzCount() >= maxLattices)
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
        for (int i = 1; i < spotVectors.size() && i < maxSearch && continuing && indexingFailureCount < 10; i++)
        {
            SpotVectorPtr spotVector1 = spotVectors[i];
            
            for (int j = 0; j < i && j < spotVectors.size() && continuing && indexingFailureCount < 10; j++)
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
                        
                        if (status == IndexingSolutionTrialSuccess)
                        {
                            i = 1;
                            j = 0;
                        }
                    }
                    
                    if (successes >= maxSuccesses || mtzCount() >= maxLattices)
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
                reachedTime = true;
                break;
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
            
            if (reachedTime)
            {
                logged << "N: Time limit reached on image " << getFilename() << " on " << mtzCount() << " crystals and " << spotCount() << " remaining spots." << std::endl;
                sendLog();
                return;
            }
        }
        else
        {
            break;
        }
    }
    
    logged << "N: Finished image " << getFilename() << " on " << mtzCount() << " crystals and " << spotCount() << " remaining spots." << std::endl;
    sendLog();
}


void Image::writeSpotsList(std::string spotFile)
{
    std::string spotDat;
    
    if (spotFile == "")
    {
        if (spotsFile.length())
        {
            spotFile = spotsFile;
        }
        else
        {
            std::string tag = "remaining";
            std::string basename = getBasename();
            
            spotFile = "_" + basename + "_" + tag + "_spots.list";
            spotDat = basename + "_spots.csv";
        }
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


double Image::resolutionAtPixel(double x, double y)
{
    vec reciprocal;
    Detector::getMaster()->findDetectorAndSpotCoordToAbsoluteVec(x, y, &reciprocal);
    scale_vector_to_distance(&reciprocal, 1 / getWavelength());
    reciprocal.l -= 1 / getWavelength();
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

void Image::drawSpotsOnPNG(std::string filename, bool drawPanels)
{
    if (!loadedSpots)
    {
        processSpotList();
    }
    
    if (!filename.length())
        filename = getBasename() + ".png";
    
    int height = FileParser::getKey("PNG_HEIGHT", 2400);
    PNGFilePtr file = PNGFilePtr(new PNGFile(filename, height, height));
    writePNG(file);
    
    
    for (int i = 0; i < spotCount(); i++)
    {
        Coord coord = spot(i)->getXY();
        
        if (!spot(i)->isFake())
        {
            spot(i)->integrate();
            
            if (!spot(i)->isBeamCentre())
            {
           //     file->drawText(i_to_str(intensity), coord.first, coord.second - 20);
            }
            else
            {
                file->drawText("Beam centre", coord.first, coord.second - 40);

            }
            file->drawCircleAroundPixel(coord.first, coord.second, 10, 1, 0, 0, 0);
        }
        else
        {
            int xLeft = coord.first - 6;
            int xRight = coord.first + 6;
            int yTop = coord.second - 6;
            int yBottom = coord.second + 6;
            
            std::string text = spot(i)->getText();
            
            file->drawText(text, coord.first, coord.second - 20);
            
            file->drawLine(xLeft, yTop, xRight, yBottom, 0.0, 0, 0, 0);
            file->drawLine(xLeft, yBottom, xRight, yTop, 0.0, 0, 0, 0);
        }
    }
    
    if (drawPanels)
    {
        std::vector<DetectorPtr> detectors;
        Detector::getMaster()->getAllSubDetectors(detectors, true);
        
        for (int i = 0; i < detectors.size(); i++)
        {
            vec xyz = detectors[i]->midPointOffsetFromParent();
            double x = xyz.h;
            double y = xyz.k;
            
            file->drawText(detectors[i]->getTag(), x, y);
        }
    }
    
    file->writeImageOutput();
    
    logged << "Written file " << filename << std::endl;
    sendLog();
}

void Image::drawMillersOnPNG(PNGFilePtr file, MtzPtr myMtz, char red, char green, char blue)
{
    int count = 0;
    bool addShoebox = FileParser::getKey("PNG_SHOEBOX", false);
    bool strongOnly = FileParser::getKey("PNG_STRONG_ONLY", false);
    
    for (int i = 0; i < myMtz->reflectionCount(); i++)
    {
        for (int j = 0; j < myMtz->reflection(i)->millerCount(); j++)
        {
            MillerPtr myMiller = myMtz->reflection(i)->miller(j);
            
            Coord pngCoord;
            
			vec intersection = myMiller->getShiftedRay();

			pngCoord = std::make_pair(intersection.h, intersection.k);

            bool strong = myMiller->reachesThreshold();
            
            if (myMiller->isSpecial())
            {
                std::string text = prettyDesc(myMiller->getHKL());
                file->drawText(text, pngCoord.first, pngCoord.second - 20);
            }
            
            if (strongOnly && !strong)
            {
                continue;
            }
            
            file->drawCircleAroundPixel(pngCoord.first, pngCoord.second, 14, (strong ? 2 : 1), red, green, blue, (strong ? 4 : 1));
            count++;
            
            if (addShoebox)
            {
                file->drawShoeboxAroundPixel(pngCoord.first, pngCoord.second, myMiller->getShoebox());
            }
        }
    }
    
    logged << "Written " << count << " reflections." << std::endl;
    sendLog();
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
    
    file->writeImageOutput();
    
    logged << "Written file " << filename << std::endl;
    sendLog();
}

void Image::writePNG(PNGFilePtr file)
{
    file->setCentre(0, 0);
    
    double threshold = FileParser::getKey("PNG_THRESHOLD", 0.);
    
    if (!FileParser::hasKey("PNG_THRESHOLD"))
    {
        threshold = standardDeviationOfPixels() * 5;
    }
    
    std::cout << "Image " << getBasename() << "- (threshold for pixel intensities: " << threshold << ")" << std::endl;
    
    DetectorPtr lastDetector = DetectorPtr();
    double minZ = 0;
    double maxZ = 0;
    double halfTol = 0.5 / mmPerPixel;
    
    Detector::getMaster()->zLimits(&minZ, &maxZ);
    
    if (averageZ <= 0)
    {
        averageZ = (maxZ + minZ) / 2;
    }
    
    logged << "Image " << getFilename() << " (min/max Z: " << minZ * mmPerPixel << ", " << maxZ * mmPerPixel << ")" << std::endl;
    sendLog();
    
    minZ = averageZ - halfTol;
    maxZ = averageZ + halfTol;
    
    for (int j = 0; j < yDim; j++)
    {
        for (int i = 0; i < xDim; i++)
        {
            bool good = accepted(i, j);
            
            if (!good)
            {
                continue;
            }
            
            double value = valueAt(i, j);
            
            if (value < 0)
                value = 0;
            
            unsigned char pixelValue = std::min(value, threshold) * 255 / threshold;
            pixelValue = 255 - pixelValue;
            float brightness = 1 - std::min(value, threshold) / threshold;
            
            int pos = xDim * j + i;
            
            DetectorPtr det = perPixelDetectors[pos];
            
            vec arranged;
            det->spotCoordToAbsoluteVec(i, j, &arranged);
            double proportion_distance = (arranged.l - minZ) / (maxZ - minZ);
            
            
            proportion_distance *= 360;
            
            png_byte red, green, blue;
            
            PNGFile::HSB_to_RGB(proportion_distance, 1.0, brightness, &red, &green, &blue);
            
            file->setPixelColourRelative(arranged.h, arranged.k, red, green, blue);
        }
    }
}

void Image::plotVectorsOnPNG(std::vector<SpotVectorPtr> vectors, std::string aFilename)
{
    int height = FileParser::getKey("PNG_HEIGHT", 2400);
    
    if (aFilename == "")
    {
        aFilename = getBasename() + "_vec.png";
    }

    PNGFilePtr file = PNGFilePtr(new PNGFile(aFilename, height, height));
    writePNG(file);
    int cx, cy;
    file->getCentre(&cx, &cy);
    
    for (int k = 0; k < vectors.size(); k++)
    {
        SpotVectorPtr spotVec = vectors[k];
        
        spotVec->drawOnImage(file);
    }
    
    file->writeImageOutput();
    
    logged << "Written file " << getFilename() << std::endl;
    sendLog();
}

void Image::plotTakeTwoVectors(std::vector<ImagePtr> images)
{
    std::vector<SpotVectorPtr> vecs;
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->mtzCount(); j++)
        {
            MtzPtr anMtz = images[i]->mtz(j);
            IndexingSolutionPtr solution = anMtz->getIndexingSolution();
            
            if (!solution)
                continue;
            
            SpotVectorMap::iterator it = solution->spotVectors.begin();
            
            for (int k = 0; k < solution->spotVectorCount(); k++)
            {
                SpotVectorPtr spotVec = it->first;
                vecs.push_back(spotVec);
                it++;
            }
        }
    }
    
    plotVectorsOnPNG(vecs);
}

