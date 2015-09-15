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
#include "CommonCircle.h"
#include "csymlib.h"

class CommonLine;
class Indexer;
class ImageCluster;

class Image
{
private:
    int pixelCountCutoff;
	std::string filename;
	vector<int> data;
    vector<unsigned char> overlapMask;
	void loadImage();
    vector<IndexerPtr> indexers;
    bool shouldMaskValue;
    bool maskedValue;
    bool fitBackgroundAsPlane;
    std::string spotsFile;
    std::ostringstream logged;
    
	/* Shoebox must be n by n where n is an odd number */
	int shoebox[7][7];

	int xDim;
	int yDim;

	int beamX;
	int beamY;
	double mmPerPixel;
    bool noCircles;

	double detectorDistance; // mm
	double wavelength;
	bool pinPoint;

    vector<SpotPtr> spots;
    vector<CommonLinePtr> commonLines;
    vector<CommonLinePair> commonLinePairs;
    vector<CommonCirclePtr> commonCircles;
    vector<CommonCirclePair> commonCirclePairs;
    double commonCircleThreshold;
    bool _hasSeeded;
    std::map<ImageCluster *, bool> unexpectedMatches;
    
	vector<vector<int> > masks;
	vector<vector<int> > spotCovers;

	int shoeboxLength();
	Mask flagAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
    double integrateFitBackgroundPlane(int x, int y, ShoeboxPtr shoebox, double *error);
    double integrateSimpleSummation(int x, int y, ShoeboxPtr shoebox, double *error);
	double integrateWithShoebox(int x, int y, ShoeboxPtr shoebox, double *error);
	bool checkShoebox(ShoeboxPtr shoebox, int x, int y);
    double weightAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
public:
    void incrementOverlapMask(int x, int y, ShoeboxPtr shoebox);
    void incrementOverlapMask(int x, int y);
    void processSpotList();
    
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
    void setUpIndexer(MatrixPtr unitcell, MatrixPtr rotation);
	std::string filenameRoot();
	void printBox(int x, int y, int tolerance);
	void addMask(int startX, int startY, int endX, int endY);
	void addSpotCover(int startX, int startY, int endX, int endY);
	bool coveredBySpot(int x, int y);
	static void applyMaskToImages(vector<Image *> images, int startX,
			int startY, int endX, int endY);
    void refineDistances();
    
    double commonCircleWithImage(Image *otherImage, CommonCirclePair *pair);
    void commonCirclesWithImages(std::vector<Image *>images);
    bool hasCommonCircleWithImage(Image *otherImage);
    void addCommonCirclePair(CommonCirclePtr thisLine, CommonCirclePtr otherLine);
    static void commonCirclesWithImagesThreaded(Image *currentImage, std::vector<Image *> images, std::vector<CommonCirclePair> *pairs, int offset);
    void commonCirclesWithImagesOneThread(std::vector<Image *> images, std::vector<CommonCirclePair> *pairs, int offset);
    CommonCirclePair getCommonCircleWithImage(Image *otherImage);
    void rotatedSpotPositions(MatrixPtr rotationMatrix, std::vector<vec> *spotPositions, std::vector<std::string> *spotElements);

    void generateCommonCircles();
    
    void addUnexpectedMatch(ImageCluster *theCluster)
    {
        unexpectedMatches[theCluster] = true;
    }
    
    void removeUnexpectedMatch(ImageCluster *theCluster)
    {
        if (unexpectedMatches.count(theCluster) > 0)
        {
            unexpectedMatches.erase(theCluster);
        }
    }
    
    bool isUnexpectedMatch(ImageCluster *theCluster)
    {
        return (unexpectedMatches.count(theCluster) > 0) && (unexpectedMatches[theCluster] == true);
    }
    
    bool hasAnyUnexpectedMatch()
    {
        return unexpectedMatches.size() > 0;
    }
    
    unsigned long circlePairCount()
    {
        return commonCirclePairs.size();
    }
    
    CommonCirclePair circlePair(int i)
    {
        return commonCirclePairs[i];
    }
    
	const std::string& getFilename() const
	{
		return filename;
	}

	void setFilename(const std::string& filename)
	{
		this->filename = filename;
	}
    
    std::string getSpotsFile()
    {
        return spotsFile;
    }
    
    void setSpotsFile(std::string newFile)
    {
        spotsFile = newFile;
    }
    
    int commonCircleCount()
    {
        return (int)commonCircles.size();
    }
    
    CommonCirclePtr getCommonCircle(int i)
    {
        return commonCircles[i];
    }
    
    void addProtoCircle(double x, double y)
    {
        CommonCirclePtr newCircle = CommonCirclePtr(new CommonCircle(this, x, y));
        commonCircles.push_back(newCircle);
    }

	int valueAt(int x, int y);
	bool accepted(int x, int y);
	double intensityAt(int x, int y, ShoeboxPtr shoebox, double *error, int tolerance = 0);

	void index();
	void refineIndexing(MtzManager *reference);
	void refineOrientations();
	vector<MtzPtr> currentMtzs();
	bool isLoaded();
    
    void setSpaceGroup(CSym::CCP4SPG *spg);
    void setMaxResolution(double res);
    void setSearchSize(int searchSize);
    void setIntensityThreshold(double threshold);
    void setUnitCell(vector<double> dims);
    void setInitialStep(double step);
    void setTestSpotSize(double spotSize);
    void setTestBandwidth(double bandwidth);
    void setOrientationTolerance(double newTolerance);
    
    bool checkUnitCell(double trueA, double trueB, double trueC, double tolerance);
    
    void setNoCircles(bool no = true)
    {
        noCircles = no;
    }
    
    bool hasSeeded()
    {
        return _hasSeeded;
    }
    
    void setSeeded(bool seed = true)
    {
        _hasSeeded = seed;
    }
    
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
   
    void setImageData(vector<int> newData);
    
    void sendLog(LogLevel priority = LogLevelNormal);

};

#endif /* IMAGE_H_ */
