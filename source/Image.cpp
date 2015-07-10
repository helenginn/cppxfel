/*
 * Image.cpp
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#include "Image.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "FileReader.h"
#include "Logger.h"
#include "FileParser.h"
#include "Panel.h"
#include "Shoebox.h"

Image::Image(std::string filename, double wavelength,
		double distance)
{
    vector<double> dims = FileParser::getKey("DETECTOR_SIZE", vector<double>());
    
	xDim = 1765;
	yDim = 1765;

    if (dims.size())
    {
        xDim = dims[0];
        yDim = dims[1];
    }
        
    data = vector<int>();
    mmPerPixel = FileParser::getKey("MM_PER_PIXEL", MM_PER_PIXEL);
    vector<double> beam = FileParser::getKey("BEAM_CENTRE", vector<double>());

    shouldMaskValue = FileParser::hasKey("IMAGE_MASKED_VALUE");

    if (shouldMaskValue)
        maskedValue = FileParser::getKey("IMAGE_MASKED_VALUE", 0);

    if (beam.size() == 0)
    {
        beam.push_back(BEAM_CENTRE_X);
        beam.push_back(BEAM_CENTRE_Y);
    }
    
    beamX = beam[0];
    beamY = beam[1];
    
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

void Image::setUpIndexer(MatrixPtr matrix)
{
	IndexerPtr indexer = IndexerPtr(new Indexer(this, matrix));
    
    if (matrix->isComplex())
        indexer->setComplexMatrix();

    indexers.push_back(indexer);
}

void Image::setUpIndexer(MatrixPtr unitcell, MatrixPtr rotation)
{
    IndexerPtr indexer = IndexerPtr(new Indexer(this, MatrixPtr(new Matrix())));
    
    indexers.push_back(indexer);
}

Image::~Image()
{
	data.clear();
	vector<int>().swap(data);
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

void Image::applyMaskToImages(vector<Image *> images, int startX,
		int startY, int endX, int endY)
{
	for (int i = 0; i < images.size(); i++)
	{
		images[i]->addMask(startX, startY, endX, endY);
	}
}

bool Image::isLoaded()
{
	return (data.size() > 0);
}

void Image::loadImage()
{
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
        bool shortInt = (bitsPerPixel == 16);
        
        if (!shortInt)
            data.resize(memblock.size() / sizeof(int));
        else
            data.resize(memblock.size() / sizeof(short));
        
        overlapMask = vector<unsigned char>(memblock.size(), 0);

        std::ostringstream logged;
		logged << "Image size: " << memblock.size() << " for image: "
				<< filename << std::endl;
        Logger::mainLogger->addStream(&logged);
        
        if (!shortInt)
            memcpy(&data[0], &memblock[0], memblock.size());
        
        if (shortInt)
        {
            for (int i = 0; i < data.size(); i++)
            {
                unsigned short point = *(&memblock[i * sizeof(short)]);
                int convertedPoint = point;
                data[i] = convertedPoint;
            }
            
            addMask(968, 950, 1441, 1079);
        }
	}
	else
        Logger::mainLogger->addString("Unable to open file");
}

void Image::dropImage()
{
	data.clear();
	vector<int>().swap(data);
    
    overlapMask.clear();
    vector<unsigned char>().swap(overlapMask);

    for (int i = 0; i < indexerCount(); i++)
		getIndexer(i)->dropMillers();
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

int Image::valueAt(int x, int y)
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

    int position = y * yDim + x;
    
	if (position < 0 || position >= data.size())
		return 0;
    
	return data[position];
}

void Image::focusOnSpot(int *x, int *y, int tolerance1, int tolerance2)
{
	Spot *spot = new Spot(this);
	spot->makeProbe(150, 5);

	double maxLift = 0;
	double maxX = *x;
	double maxY = *y;

	for (int i = *x - tolerance1; i <= *x + tolerance1; i++)
	{
		for (int j = *y - tolerance1; j <= *y + tolerance1; j++)
		{
			double lift = spot->maximumLift(this, i, j, true);

			if (lift > maxLift)
			{
				maxLift = lift;
				maxX = i;
				maxY = j;
			}
		}
	}

	*x = maxX;
	*y = maxY;
}

void Image::focusOnAverageMax(int *x, int *y, int tolerance1, int tolerance2)
{
//	focusOnSpot(x, y, tolerance1, tolerance2);
//	return;

	int maxValue = 0;
	int newX = *x;
	int newY = *y;

	for (int i = *x - tolerance1; i <= *x + tolerance1; i++)
	{
		for (int j = *y - tolerance1; j <= *y + tolerance1; j++)
		{
			double newValue = 0;

			for (int h = i - tolerance2; h <= i + tolerance2; h++)
			{
				for (int k = j - tolerance2; k <= j + tolerance2; k++)
				{
					int addition = valueAt(h, k);
					newValue += addition;
				}
			}

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

    if (tolerance > 0 && accepted(newX, newY))
        printBox(newX, newY, tolerance);
    
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
    
    if (value < 0 && value > 1)
        return 0;
    
    return value;
}


// @TODO
bool Image::checkShoebox(ShoeboxPtr shoebox, int x, int y)
{
    int slowSide = 0;
    int fastSide = 0;
    
	shoebox->sideLengths(&slowSide, &fastSide);
	int slowTol = (slowSide - 1) / 2;
    int fastTol = (fastSide - 1) / 2;
    
    int zeroCount = 0;
    int count = 0;

	for (int i = x - slowTol; i < x + slowTol; i++)
	{
		for (int j = y - fastTol; j <= y + fastTol; j++)
		{
			if (!accepted(i, j))
			{
                std::ostringstream logged;
                logged << "Rejecting miller - too many zeros" << std::endl;
                Logger::mainLogger->addStream(&logged, LogLevelDebug);
                return false;
			}
            
            count++;
            if (valueAt(i, j) == 0)
                zeroCount++;
		}
	}
    
    if (double(zeroCount) / double(count) > 0.4)
    {
        std::ostringstream logged;
        logged << "Rejecting miller - too many zeros" << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelDebug);
        return false;
    }

	return true;
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
    
   // Logger::mainLogger->addStream(&logged, LogLevelDebug);
}

double Image::integrateWithShoebox(int x, int y, ShoeboxPtr shoebox, double *error)
{
    int slowSide = 0;
    int fastSide = 0;
    
    shoebox->sideLengths(&slowSide, &fastSide);
    int slowTol = (slowSide - 1) / 2;
    int fastTol = (fastSide - 1) / 2;
    
	int startX = x - slowTol;
	int startY = y - fastTol;

	int foreground = 0;
	int foreNum = 0;

	int background = 0;
	int backNum = 0;

//	print = true;

	for (int i = startX; i < startX + slowSide; i++)
	{
		for (int j = startY; j < startY + fastSide; j++)
		{
            double sX = i - startX;
            double sY = j - startY;
            
			Mask flag = flagAtShoeboxIndex(shoebox, sX, sY);

			if (flag == MaskForeground)
			{
                double weight = weightAtShoeboxIndex(shoebox, sX, sY);
                double value = valueAt(i, j);
				foreNum += weight;
				foreground += value * weight;
                
                if (value > pixelCountCutoff && pixelCountCutoff > 0)
                {
              //      return isnan(' ');
                }
			}
			else if (flag == MaskBackground)
			{
				backNum++;
				background += valueAt(i, j);
			}
		}
	}

    double aveBackground = (double) background / (double) backNum;
	double backgroundInForeground = aveBackground * (double) foreNum;

	double totalPhotons = foreground + background;
	*error = sqrt(totalPhotons);

    
    std::ostringstream logged;
 //   logged << "Photons (back/fore/back-in-fore): " << background << "\t" << foreground << "\t" << "; weights: " << backNum << "\t" << foreNum << std::endl;
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
    
	double intensity = (foreground - backgroundInForeground);
	
    if (!std::isfinite(intensity))
	{
        Logger::mainLogger->addString("Infinite intensity", LogLevelDebug);
	}

	return intensity;
}

double Image::intensityAt(int x, int y, ShoeboxPtr shoebox, double *error, int tolerance)
{
	int x1 = x;
	int y1 = y;

	if (tolerance > 0)
	{
		if (pinPoint)
			focusOnMaximum(&x1, &y1, tolerance);
		else
			focusOnAverageMax(&x1, &y1, tolerance, 2);
	}

	if (checkShoebox(shoebox, x1, y1) == 0)
	{
		return nan("");
	}

	double integral = integrateWithShoebox(x1, y1, shoebox, error);

	return integral;
}

bool Image::accepted(int x, int y)
{
    if (shouldMaskValue)
    {
        double value = valueAt(x, y);
        
        if (value == maskedValue)
            return false;
    }
    
    Coord coord = std::make_pair(x, y);
    
    PanelPtr panel = Panel::panelForCoord(coord);
    
    return !panel ? false : true;
}

void Image::index()
{
	if (indexers.size() == 0)
    {
        std::cout << "No orientation matrices, cannot index/integrate." << std::endl;
        return;
    }

    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->checkAllMillers(indexers[i]->getMaxResolution(),
			indexers[i]->getTestBandwidth());
    }
}

void Image::refineIndexing(MtzManager *reference)
{
    if (indexers.size() == 0)
    {
        std::cout << "No orientation matrices, cannot index/integrate." << std::endl;
        return;
    }
    
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->refineDetectorAndWavelength(reference);
    }
}

void Image::refineOrientations()
{
    if (indexers.size() == 0)
    {
        std::cout << "No orientation matrices, cannot index/integrate." << std::endl;
        return;
    }
    
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->refineOrientationMatrix();
    }
}

void Image::refineDistances()
{
    if (indexers.size() == 0)
    {
        std::cout << "No orientation matrices, cannot refine distance." << std::endl;
        return;
    }
    
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->refineDetectorAndWavelength();
    }
}

vector<MtzPtr> Image::currentMtzs()
{
    vector<MtzPtr> mtzs;
    int count = 0;
    
    for (int i = 0; i < indexerCount(); i++)
    {
        MtzPtr newMtz = indexers[i]->newMtz(i);
        count += newMtz->rejectOverlaps();
        mtzs.push_back(newMtz);
    }
    
    std::cout << "Generated " << mtzs.size() << " mtzs with " << count << " rejected reflections." << std::endl;

    // fix me
    
	return mtzs;
}

void Image::setSpaceGroup(CCP4SPG *spg)
{
    for (int i = 0; i < indexerCount(); i++)
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

void Image::setIntensityThreshold(double threshold)
{
    for (int i = 0; i < indexers.size(); i++)
    {
        indexers[i]->setIntensityThreshold(threshold);
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
    for (int i = 0; i < indexerCount(); i++)
    {
        bool isOk = true;
        double *cellDims = new double[3];
        getIndexer(i)->getMatrix()->unitCellLengths(&cellDims);
        
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
    
    return indexerCount() > 0;
}

