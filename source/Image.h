/*
 * Image.h
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#ifndef IMAGE_H_
#define IMAGE_H_

#include "Matrix.h"
#include "IOMRefiner.h"
#include "parameters.h"
#include "Logger.h"
#include "csymlib.h"
#include "SpotVector.h"
#include "LoggableObject.h"

typedef enum
{
    ImageClassCppxfel,
    ImageClassDIALS,
    ImageClassHdf5,
} ImageClass;

class IOMRefiner;
class ImageCluster;

typedef enum
{
    IndexingSolutionTrialSuccess,
    IndexingSolutionTrialFailure,
    IndexingSolutionTrialDuplicate,
    IndexingSolutionBranchFailure,
} IndexingSolutionStatus;

class Image : protected LoggableObject, public boost::enable_shared_from_this<Image>
{
private:
    int pixelCountCutoff;
	std::string filename;
    std::vector<MtzPtr> mtzs;
    
    vector<unsigned char> overlapMask;
	virtual void loadImage();
    vector<IOMRefinerPtr> indexers;
    vector<IOMRefinerPtr> failedRefiners;
    double metrologyMoveThreshold;
    bool shouldMaskValue;
    bool shouldMaskUnderValue;
    int maskedValue;
    int maskedUnderValue;
    bool fitBackgroundAsPlane;
    std::string spotsFile;
    IndexingSolutionStatus extendIndexingSolution(IndexingSolutionPtr solutionPtr, std::vector<SpotVectorPtr> existingVectors, int *failures = NULL, int added = 0);
    
	/* Shoebox must be n by n where n is an odd number */
	int shoebox[7][7];

	double beamX;
	double beamY;
	double mmPerPixel;
    double detectorGain;

    static double globalDetectorDistance;
    static double globalBeamX;
    static double globalBeamY;
    
	double detectorDistance; // mm
	double wavelength;
	bool pinPoint;

    int indexingFailureCount;
    
    std::vector<IndexingSolutionPtr> goodSolutions;
    std::vector<IndexingSolutionPtr> badSolutions;
    std::vector<SpotVectorPtr> spotVectors;
    bool _hasSeeded;
    
	vector<vector<int> > masks;
	vector<vector<int> > spotCovers;

	int shoeboxLength();
	Mask flagAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
    double integrateFitBackgroundPlane(int x, int y, ShoeboxPtr shoebox, float *error);
    double integrateSimpleSummation(int x, int y, ShoeboxPtr shoebox, float *error);
	double integrateWithShoebox(int x, int y, ShoeboxPtr shoebox, float *error);
	bool checkShoebox(ShoeboxPtr shoebox, int x, int y);
    double weightAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
    IndexingSolutionStatus testSeedSolution(IndexingSolutionPtr newSolution, std::vector<SpotVectorPtr> &prunedVectors, int *successes);
    IndexingSolutionPtr biggestFailedSolution;
    std::vector<SpotVectorPtr> biggestFailedSolutionVectors;
    double resolutionAtPixel(double x, double y);
protected:
    int xDim;
    int yDim;
    
    double standardDeviationOfPixels();
    std::vector<SpotPtr> spots;
    virtual IndexingSolutionStatus tryIndexingSolution(IndexingSolutionPtr solutionPtr, bool modify, int *spotsRemoved, IOMRefinerPtr *refiner);
    virtual IndexingSolutionStatus tryIndexingSolution(IndexingSolutionPtr solutionPtr);
    virtual bool checkIndexingSolutionDuplicates(MatrixPtr newSolution, bool excludeLast = false);
    int minimumSolutionNetworkCount;
    bool loadedSpots;

    // this really ought to be a template
    vector<short> shortData;
    vector<int> data;
    bool useShortData;
    // end of should be a template
    void writePNG(PNGFilePtr file);

public:
    void incrementOverlapMask(int x, int y, ShoeboxPtr shoebox);
    void incrementOverlapMask(int x, int y);
    virtual void processSpotList();
    virtual void writeSpotsList(std::string spotFile = "");
    virtual std::pair<double, double> reciprocalCoordinatesToPixels(vec hkl, double myWavelength = 0);
    virtual vec pixelsToReciprocalCoordinates(double xPix, double yPix);
    virtual vec millimetresToReciprocalCoordinates(double xmm, double ymm);
    
    unsigned char overlapAt(int x, int y);
    unsigned char maximumOverlapMask(int x, int y, ShoeboxPtr shoebox);
	Image(std::string filename = "", double wavelength = 0,
			double distance = 0);
	void focusOnSpot(int *x, int *y, int tolerance1, int tolerance2);
	void focusOnAverageMax(int *x, int *y, int tolerance1, int tolerance2 = 1, bool even = false);
    void focusOnMaximum(int *x, int *y, int tolerance = 0, double shiftX = 0, double shiftY = 0);
	void dropImage();
    void newImage();
	virtual ~Image();
	void setUpIOMRefiner(MatrixPtr matrix);
    void setUpIOMRefiner(MatrixPtr unitcell, MatrixPtr rotation);
	std::string filenameRoot();
	void printBox(int x, int y, int tolerance);
	void addMask(int startX, int startY, int endX, int endY);
	void addSpotCover(int startX, int startY, int endX, int endY);
	bool coveredBySpot(int x, int y);
	static void applyMaskToImages(vector<ImagePtr> images, int startX,
			int startY, int endX, int endY);
    std::vector<double> anglesBetweenVectorDistances(double distance1, double distance2, double tolerance);
    void reset();
    void findSpots();

    void rotatedSpotPositions(MatrixPtr rotationMatrix, std::vector<vec> *spotPositions, std::vector<std::string> *spotElements);
    void radialAverage();
    void integrateSpots();
    void drawMillersOnPNG(PNGFilePtr file, MtzPtr myMtz, char red = 0, char green = 0, char blue = 0);
    void drawCrystalsOnPNG(int crystalNum);
    void drawSpotsOnPNG(std::string filename = "");
    void dumpImage();
    void makeMaximumFromImages(std::vector<ImagePtr> images);
    void excludeWeakestSpots(double fraction);
    
	const std::string& getFilename() const
	{
		return filename;
	}

	void setFilename(const std::string& filename)
	{
		this->filename = filename;
	}
    
    std::string getBasename()
    {
        int fullStopIndex = (int)filename.rfind(".");
        if (fullStopIndex == std::string::npos)
            return filename;
        
        std::string basename = filename.substr(0, fullStopIndex);
        
        return basename;
    }
    
    std::string getSpotsFile()
    {
        return spotsFile;
    }
    
    void setSpotsFile(std::string newFile)
    {
        spotsFile = newFile;
    }

	int valueAt(int x, int y);
    int rawValueAt(int x, int y);
    void addValueAt(int x, int y, int addedValue);
	bool accepted(int x, int y);
	double intensityAt(int x, int y, ShoeboxPtr shoebox, float *error, int tolerance = 0);

	void index();
	void refineIndexing(MtzManager *reference);
	void refineOrientations();
	vector<MtzPtr> currentMtzs();
    std::vector<MtzPtr> getLastMtzs();
	bool isLoaded();
    
    void setSpaceGroup(CSym::CCP4SPG *spg);
    void setMaxResolution(double res);
    void setSearchSize(int searchSize);
    void setUnitCell(vector<double> dims);
    void setInitialStep(double step);
    void setTestSpotSize(double spotSize);
    void setTestBandwidth(double bandwidth);
    void setOrientationTolerance(double newTolerance);
    
    bool checkUnitCell(double trueA, double trueB, double trueC, double tolerance);
    
    virtual void findIndexingSolutions();
    void compileDistancesFromSpots(double maxReciprocalDistance = 0, double tooCloseDistance = 0, bool filter = false);
    void filterSpotVectors();
    void plotTakeTwoVectors(std::vector<ImagePtr> images);
    int throwAwayIntegratedSpots(std::vector<MtzPtr> mtzs);
    void updateAllSpots();
    void clusterCountWithSpotNumber(int spotNum);
    bool acceptableSpotCount();

    void addSpotIfNotMasked(SpotPtr newSpot);
    
    static void setGlobalDetectorDistance(void *object, double distance)
    {
        globalDetectorDistance = distance;
    }
    
    static double getGlobalDetectorDistance(void *object)
    {
        return globalDetectorDistance;
    }
    
    static double getGlobalBeamX(void *object)
    {
        return globalBeamX;
    }
    
    static void setGlobalBeamX(void *object, double newX)
    {
        globalBeamX = newX;
    }
    
    static double getGlobalBeamY(void *object)
    {
        return globalBeamY;
    }
    
    static void setGlobalBeamY(void *object, double newY)
    {
        globalBeamY = newY;
    }
    
    void removeRefiner(int j)
    {
        indexers.erase(indexers.begin() + j);
    }
    
    void setIOMRefiners(std::vector<IOMRefinerPtr> refiners)
    {
        indexers = refiners;
    }
    
    int spotVectorCount()
    {
        return (int)spotVectors.size();
    }
    
    SpotPtr spot(int i)
    {
        return spots[i];
    }
    
    int spotCount()
    {
        return (int)spots.size();
    }
    SpotVectorPtr spotVector(int i)
    {
        return spotVectors[i];
    }
    
    bool hasSeeded()
    {
        return _hasSeeded;
    }
    
    void setSeeded(bool seed = true)
    {
        _hasSeeded = seed;
    }
    
    int IOMRefinerCount()
    {
        return (int)indexers.size();
    }
    
    IOMRefinerPtr getIOMRefiner(int i)
    {
        return indexers[i];
    }
    
    void addIOMRefiner(IOMRefinerPtr newRefiner)
    {
        indexers.push_back(newRefiner);
    }
    
    void clearIOMRefiners()
    {
        indexers.clear();
        std::vector<IOMRefinerPtr>().swap(indexers);
    }
    
    void clearMtzs()
    {
        mtzs.clear();
        std::vector<MtzPtr>().swap(mtzs);
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

	double getDetectorDistance()
	{
        if (detectorDistance == 0)
        {
            return globalDetectorDistance;
        }
        
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

	double getBeamX() const
	{
        if (beamX == INT_MAX)
        {
            return globalBeamX;
        }
        
		return beamX;
	}

	void setBeamX(double beamX)
	{
		this->beamX = beamX;
	}

	double getBeamY() const
	{
        if (beamY == INT_MAX)
        {
            return globalBeamY;
        }
        
		return beamY;
	}

	void setBeamY(double beamY)
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
    
    double getDetectorGain()
    {
        return detectorGain;
    }
    
    void setDetectorGain(double newGain)
    {
        detectorGain = newGain;
    }
    
    IndexingSolutionPtr getGoodOrBadSolution(int i, bool good = true)
    {
        if (good) return goodSolutions[i];
        else return badSolutions[i];
    }
    
    int goodOrBadSolutionCount(bool good = true)
    {
        if (good) return (int)goodSolutions.size();
        else return (int)badSolutions.size();
    }
    
    int failedRefinerCount()
    {
        return (int)failedRefiners.size();
    }
    
    IOMRefinerPtr failedRefiner(int i)
    {
        return failedRefiners[i];
    }
    
    virtual ImageClass getClass()
    {
        return ImageClassCppxfel;
    }
    
    void addMtz(MtzPtr newMtz)
    {
        mtzs.push_back(newMtz);
    }
    
    
    int mtzCount()
    {
        return (int)mtzs.size();
    }
    
    MtzPtr mtz(int i)
    {
        return mtzs[i];
    }
    
    short int *getShortDataPtr()
    {
        if (!shortData.size())
        {
            return NULL;
        }
        
        return &shortData[0];
    }
    
    int *getDataPtr()
    {
        if (!data.size())
        {
            return NULL;
        }
        
        return &data[0];
    }
};

#endif /* IMAGE_H_ */
