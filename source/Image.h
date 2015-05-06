/*
 * Image.h
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#ifndef IMAGE_H_
#define IMAGE_H_

#include "Matrix.h"
#include "Indexer.h"
#include "parameters.h"
#include "Logger.h"
#include "headers/csymlib.h"

class Indexer;
class MtzManager;

class Image
{
private:
    int pixelCountCutoff;
	std::string filename;
	std::vector<int> data;
    std::vector<unsigned char> overlapMask;
	void loadImage();
    std::vector<IndexerPtr> indexers;
    bool shouldMaskValue;
    bool maskedValue;
    
	/* Shoebox must be n by n where n is an odd number */
	int shoebox[7][7];

//	int customShoeboxLength;
//	ShoeboxPtr customShoebox;

	int xDim;
	int yDim;

	int beamX;
	int beamY;
	double mmPerPixel;

	double detectorDistance; // mm
	double wavelength;
	bool pinPoint;

	std::vector<std::vector<int> > masks;
	std::vector<std::vector<int> > spotCovers;

	int shoeboxLength();
	Mask flagAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
	double integrateWithShoebox(int x, int y, ShoeboxPtr shoebox, double *error);
	bool checkShoebox(ShoeboxPtr shoebox, int x, int y);
    double weightAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
public:
    void incrementOverlapMask(int x, int y, ShoeboxPtr shoebox);
    void incrementOverlapMask(int x, int y);
    unsigned char overlapAt(int x, int y);
    unsigned char maximumOverlapMask(int x, int y, ShoeboxPtr shoebox);
	Image(std::string filename = "", double wavelength = 0,
			double distance = 0);
	void focusOnSpot(int *x, int *y, int tolerance1, int tolerance2);
	void focusOnAverageMax(int *x, int *y, int tolerance1, int tolerance2);
	void focusOnMaximum(int *x, int *y, int tolerance = 0, double shiftX = 0, double shiftY = 0);
	void dropImage();
	virtual ~Image();
	void setUpIndexer(MatrixPtr matrix);
	std::string filenameRoot();
	void printBox(int x, int y, int tolerance);
	void addMask(int startX, int startY, int endX, int endY);
	void addSpotCover(int startX, int startY, int endX, int endY);
	bool coveredBySpot(int x, int y);
	static void applyMaskToImages(std::vector<Image *> images, int startX,
			int startY, int endX, int endY);
    void refineDistances();
    
	const std::string& getFilename() const
	{
		return filename;
	}

	void setFilename(const std::string& filename)
	{
		this->filename = filename;
	}

	int valueAt(int x, int y);
	bool accepted(int x, int y);
	double intensityAt(int x, int y, ShoeboxPtr shoebox, double *error, int tolerance = 0);

	void index();
	void refineIndexing(MtzManager *reference);
	void refineOrientations();
	std::vector<MtzPtr> currentMtzs();
	bool isLoaded();
    
    void setSpaceGroup(CSym::CCP4SPG *spg);
    void setMaxResolution(double res);
    void setSearchSize(int searchSize);
    void setIntensityThreshold(double threshold);
    void setUnitCell(std::vector<double> dims);
    void setInitialStep(double step);
    void setTestSpotSize(double spotSize);
    void setTestBandwidth(double bandwidth);
    
    int indexerCount()
    {
        return (int)indexers.size();
    }
    
    IndexerPtr getIndexer(int i)
    {
        return indexers[i];
    }
    
	int getXDim() const
	{
		return xDim;
	}

	void setXDim(int dim)
	{
		xDim = dim;
	}

	int getYDim() const
	{
		return yDim;
	}

	void setYDim(int dim)
	{
		yDim = dim;
	}

	double getDetectorDistance() const
	{
		return detectorDistance;
	}

	void setDetectorDistance(double detectorDistance)
	{
		this->detectorDistance = detectorDistance;
	}

	double getWavelength() const
	{
		return wavelength;
	}

	void setWavelength(double wavelength)
	{
		this->wavelength = wavelength;
	}

	int getBeamX() const
	{
		return beamX;
	}

	void setBeamX(int beamX)
	{
		this->beamX = beamX;
	}

	int getBeamY() const
	{
		return beamY;
	}

	void setBeamY(int beamY)
	{
		this->beamY = beamY;
	}

	double getMmPerPixel() const
	{
		return mmPerPixel;
	}

	void setMmPerPixel(double mmPerPixel)
	{
		this->mmPerPixel = mmPerPixel;
	}

	bool isPinPoint() const
	{
		return pinPoint;
	}

	void setPinPoint(bool pinPoint)
	{
		this->pinPoint = pinPoint;
	}
};

#endif /* IMAGE_H_ */
