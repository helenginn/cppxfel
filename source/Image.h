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
#include "hasFilename.h"

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

class Image : protected LoggableObject, public hasFilename, public boost::enable_shared_from_this<Image>
{
private:
    int pixelCountCutoff;
    std::vector<MtzPtr> mtzs;
    int highScore;
    bool fake;
    static ImagePtr _imageMask;
    static bool interpolate;
    
    virtual void loadImage();
    vector<IOMRefinerPtr> indexers;
    vector<IOMRefinerPtr> failedRefiners;
    vector<ImagePtr> maxes;

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
    double averageZ;

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

	Mask flagAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
    double integrateFitBackgroundPlane(int x, int y, ShoeboxPtr shoebox, float *error);
    double integrateSimpleSummation(double x, double y, ShoeboxPtr shoebox, float *error);
	double integrateWithShoebox(double x, double y, ShoeboxPtr shoebox, float *error);
	double weightAtShoeboxIndex(ShoeboxPtr shoebox, int x, int y);
    IndexingSolutionStatus testSeedSolution(IndexingSolutionPtr newSolution, std::vector<SpotVectorPtr> &prunedVectors, int *successes);
    IndexingSolutionPtr biggestFailedSolution;
    std::vector<SpotVectorPtr> biggestFailedSolutionVectors;
protected:
    int xDim;
    int yDim;
    bool _isMask;

    void addHighScore(int score)
    {
        highScore += score;
    }
    
    
    void checkAndSetupLookupTable();
    static std::mutex setupMutex;
    double standardDeviationOfPixels();
    std::vector<SpotPtr> spots;
    virtual IndexingSolutionStatus tryIndexingSolution(IndexingSolutionPtr solutionPtr, bool modify, int *spotsRemoved, IOMRefinerPtr *refiner);
    virtual IndexingSolutionStatus tryIndexingSolution(IndexingSolutionPtr solutionPtr);
    virtual bool checkIndexingSolutionDuplicates(MatrixPtr newSolution, bool excludeLast = false);
    int minimumSolutionNetworkCount;
    bool loadedSpots;
    vector<signed char> overlapMask;
    static vector<signed char> generalMask;
    static vector<DetectorPtr> perPixelDetectors;
    
    // this really ought to be a template
    vector<short> shortData;
    vector<int> data;
    bool useShortData;
    // end of should be a template
    void writePNG(PNGFilePtr file);
    double spotVectorWeight;

public:
    void incrementOverlapMask(int x, int y, ShoeboxPtr shoebox);
    void incrementOverlapMask(int x, int y);
    virtual void processSpotList();
    virtual void writeSpotsList(std::string spotFile = "");
    
    signed char maskValueAt(signed char *firstByte, int x, int y);
    unsigned char maximumOverlapMask(int x, int y, ShoeboxPtr shoebox);
	Image(std::string filename = "", double wavelength = 0,
			double distance = 0);
	void focusOnSpot(int *x, int *y, int tolerance1, int tolerance2);
	void focusOnAverageMax(double *x, double *y, int tolerance1, int tolerance2 = 1, bool even = false);
    void dropImage();
    void newImage();
	virtual ~Image();
	void setUpIOMRefiner(MatrixPtr matrix);
	void addMask(int startX, int startY, int endX, int endY);
	static void applyMaskToImages(vector<ImagePtr> images, int startX,
			int startY, int endX, int endY);
    std::vector<double> anglesBetweenVectorDistances(double distance1, double distance2, double tolerance);
    void findSpots();
    double resolutionAtPixel(double x, double y);

    void loadBadPixels();
    void fakeSpots();
    void integrateSpots();
    void drawMillersOnPNG(PNGFilePtr file, MtzPtr myMtz, char red = 0, char green = 0, char blue = 0);
    void drawCrystalsOnPNG(int crystalNum);
    void drawSpotsOnPNG(std::string filename = "", bool drawPanels = false);
    void dumpImage();
    void makeMaximumFromImages(std::vector<ImagePtr> images, bool listResults = false);
    void excludeWeakestSpots(double fraction);
    void plotVectorsOnPNG(std::vector<SpotVectorPtr> vectors, std::string filename = "");
    
    static void setImageMask(ImagePtr mask)
    {
        _imageMask = mask;
        std::ostringstream logged;
        logged << "Loaded mask from HDF5" << std::endl;
        Logger::log(logged);
    }
    
    static ImagePtr getImageMask()
    {
        return _imageMask;
    }
    
    int getHighScore()
    {
        return highScore;
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
    double interpolateAt(double x, double y, double *total);
    int rawValueAt(int x, int y);
    void addValueAt(int x, int y, int addedValue);
	bool accepted(int x, int y);
	double intensityAt(double x, double y, ShoeboxPtr shoebox, float *error, int tolerance = 0);

	void refineIndexing(MtzManager *reference);
	void refineOrientations();
	vector<MtzPtr> currentMtzs();
	bool isLoaded();
    
    void setSpaceGroup(CSym::CCP4SPG *spg);
    void setUnitCell(vector<double> dims);
    void weedOutCloseSpots();
    
    virtual void findIndexingSolutions();
    void compileDistancesFromSpots(double maxReciprocalDistance = 0, double tooCloseDistance = 0, bool filter = false);
    void filterSpotVectors();
    void plotTakeTwoVectors(std::vector<ImagePtr> images);
    int throwAwayIntegratedSpots(std::vector<MtzPtr> mtzs);
    void updateAllSpots();
    bool acceptableSpotCount();

    void addSpotIfNotMasked(SpotPtr newSpot);
    
    double getSpotVectorWeight()
    {
        return spotVectorWeight;
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
		return beamX;
	}

	void setBeamX(double beamX)
	{
		this->beamX = beamX;
	}

	double getBeamY() const
	{
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
