#ifndef mtz_manager
#define mtz_manager

#include <string>
#include <iostream>
#include "cmtzlib.h"
#include "csymlib.h"
#include <vector>

#include "definitions.h"
#include "parameters.h"
#include <tuple>
#include "Logger.h"
#include "Matrix.h"
#include "hasFilename.h"
#include "hasSymmetry.h"

using namespace CMtz;
using namespace CSym;

typedef enum
{
	ScoreTypeMinimizeRSplit = 0,
	ScoreTypeCorrelation = 1,
	ScoreTypeReward = 2,
    ScoreTypeSymmetry = 8,
    ScoreTypeRSplitIntensity = 11,
} ScoreType;

typedef enum
{
	TrustLevelGood, TrustLevelAverage, TrustLevelBad
} TrustLevel;

class Miller;

class MtzManager : public LoggableObject, public hasFilename, public hasSymmetry, public boost::enable_shared_from_this<MtzManager>
{

protected:
/* From IOMRefiner */
	IndexingSolutionPtr indexingSolution;
	double testBandwidth;
	double testSpotSize;
	double lastTotal;
	double lastStdev;
	double initialStep;
	double searchSize;
	double orientationTolerance;
	bool refineOrientations;

	vector<MillerPtr> nearbyMillers;

	void calculateNearbyMillers();
	void checkNearbyMillers();

	double getReflectionWavelengthStdev();
	bool millerWithinBandwidth(MillerPtr miller);
	int getTotalReflections();

	double lScore(bool silent);
	double hkScore(bool finalise);
	static double hkScoreFinalise(void *object);
	static double hkScoreStdevWrapper(void *object);
	static double lScoreWrapper(void *object);

/* End */

    std::string parentImage;
	static bool reflection_comparison(ReflectionPtr i, ReflectionPtr j);

	ImageWeakPtr image;
    bool dropped;
    void getWavelengthFromHDF5();
    
    BeamPtr beam;
    
	// minimisation stuff

    int millerCount();
    
	vector<ReflectionPtr> reflections;
	vector<ReflectionPtr> refReflections;
	vector<ReflectionPtr> matchReflections;
    MtzManager *previousReference;
    int previousAmbiguity;
    bool allowTrust;
    bool setInitialValues;
    
	double hRot;
	double kRot;
	double lRot;
	double mosaicity;
	double spotSize;
	double wavelength;
	double bandwidth;
	double exponent;
	double refCorrelation;
    double refPartCorrel;
	double scale;
	
    double timeDelay;
    
    uint32 activeAmbiguity;

	bool optimisingWavelength;
	bool optimisingBandwidth;
	bool optimisingMosaicity;
	bool optimisingRlpSize;
	bool optimisingOrientation;
	bool optimisingExponent;

	double stepSizeWavelength;
	double stepSizeBandwidth;
	double stepSizeMosaicity;
	double stepSizeRlpSize;
	double stepSizeOrientation;
    double stepSizeOrientABC;
	double stepSizeExponent;

	double toleranceWavelength;
	double toleranceBandwidth;
	double toleranceMosaicity;
	double toleranceRlpSize;
	double toleranceOrientation;
	double toleranceExponent;

	ScoreType defaultScoreType;
	double minResolutionAll;
	double maxResolutionAll;
    float lastRSplit;

    static double superGaussianScale;
    double *params;
    double detectorDistance;
    
	bool rejected;

	TrustLevel trust;
	ScoreType scoreType;

    MatrixPtr baseMatrix;
    MatrixPtr rotatedMatrix;
	static MtzManager *referenceManager;
	MtzManager *lastReference;

public:
/* From IOMRefiner */
	void dropMillers();
	bool isGoodSolution();
	void refineOrientationMatrix(bool force = false);
	void fakeSpots();
	void getWavelengthHistogram(vector<double> *wavelengths = NULL,
								vector<int> *frequencies = NULL, bool silent = false);
	void integrateChosenMillers(bool quick = false);

/* END */

    static vector<double> superGaussianTable;
    static bool setupSuperGaussian;
    static std::mutex tableMutex;
    double bFactor;
    double externalScale;
    int removeStrongSpots(std::vector<SpotPtr> *spots, bool actuallyDelete = true);
    int rejectOverlaps();
    
	MtzManager(void);
	virtual ~MtzManager(void);
	void loadParametersMap();

    void addMiller(MillerPtr miller);
    void clearReflections();
	void addReflection(ReflectionPtr reflection);
	void removeReflection(int i);
	void excludeFromLogCorrelation();

    MatrixPtr matrix;
    
    void description(void);
    void resetDefaultParameters();
    
    void setDefaultMatrix();
	void setMatrix(MatrixPtr newMat);
	void sortReflections();
	void applyUnrefinedPartiality();
    void incrementActiveAmbiguity();
    double maxResolution();

	virtual void loadReflections();
    void dropReflections();
	static void setReference(MtzManager *reference);
    ReflectionPtr findReflectionWithId(ReflectionPtr exampleRefl, size_t *lowestId = NULL);
	int findReflectionWithId(long unsigned int refl_id, ReflectionPtr *reflection, bool insertionPoint = false);
	void findCommonReflections(MtzManager *other,
			vector<ReflectionPtr> &reflectionVector1, vector<ReflectionPtr> &reflectionVector2,
			int *num = NULL, bool acceptableOnly = false, bool preserve = false);
	void scaleToMtz(MtzManager *otherManager, bool withCutoff = false, double lowRes = 0, double highRes = 0);
	void bFactorAndScale(double *scale, double *bFactor, double exponent = 1);
	void applyBFactor(double bFactor);
	void applyScaleFactor(double scaleFactor, double lowRes = 0, double highRes = 0, bool absolute = false);
	double averageIntensity();
	void setSigmaToUnity();
	void replaceBeamWithSpectrum();

	std::vector<double> getDifferencesWith(MtzPtr other, std::vector<double> &refls);
	
	void copySymmetryInformationFromManager(MtzPtr toCopy);
	void applyPolarisation(void);

	virtual void writeToFile(std::string newFilename, bool announce = false, bool plusAmbiguity = false);
    void writeToHdf5();
    
	double correlationWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false, bool freeOnly = false);

	double rSplitWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false, bool freeOnly = false);

	double statisticsWithManager(MtzManager *otherManager, StatisticsFunction *function,
			bool printHits, bool silent, double lowRes, double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog, bool freeOnly = false);

	double statisticsWithManager(MtzManager *otherManager,
			StatisticsFunction *function, RFactorFunction *rFactorFunction,
			RFactorType rFactor, bool printHits, bool silent, double lowRes,
			double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog, bool freeOnly = false);

	double rFactorWithManager(RFactorType rFactor, bool printHits = false, bool silent = true, double lowRes = 0,
			double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL,
			bool shouldLog = false);

// minimisation stuff

    void recalculateWavelengths();

    double medianWavelength(double lowRes, double highRes);
    double bestWavelength(double lowRes = 0.0, double highRes = 0, bool usingReference = false);
	int accepted(void);
	static double exclusionScoreWrapper(void *object, double lowRes = 0,
			double highRes = 0);
    double leastSquaresPartiality(double low = 0, double high = 0);
    double correlation(bool silent = true, double lowResolution = 0, double highResolution = -1);
	double rSplit(double low, double high);
	double rewardAgreement(double low, double high);
	std::string describeScoreType();
    double refinePartialitiesOrientation(int ambiguity, bool reset = true);
    
    void refinePartialities();

    void refreshCurrentPartialities();
    // delete
	void refreshPartialities(double hRot, double kRot, double mosaicity,
                             double spotSize, double wavelength,
							 double bandwidth, double exponent);

	std::vector<MillerPtr> strongMillers();
// more grid search

	std::string getOrientationMatrixListScript();
    void setParamLine(std::string line);
    std::string getParamLine();
    
    static void makeSuperGaussianLookupTable(double exponent);
    
    static std::string parameterHeaders();
    std::string writeParameterSummary();
    
    void millersToDetector();
    
    int ambiguityCount();
    void flipToActiveAmbiguity();
    void resetFlip();
    
    float getLastRSplit()
    {
        return lastRSplit;
    }
    
    static bool setupGaussianTable()
    {
        return setupSuperGaussian;
    }
    
    int getActiveAmbiguity()
    {
        return activeAmbiguity;
    }
    
    void setActiveAmbiguity(int newAmbiguity);

    virtual int reflectionCount()
    {
        return (int)reflections.size();
    }
    
    virtual ReflectionPtr reflection(int i)
    {
        return reflections[i];
    }
    
    void setDefaultScoreType(ScoreType scoreType)
    {
        defaultScoreType = scoreType;
    }
    
	double getBandwidth() const
	{
		return bandwidth;
	}

	void setBandwidth(double bandwidth)
	{
		this->bandwidth = bandwidth;
	}

	double getMosaicity() const
	{
		return mosaicity;
	}

	void setMosaicity(double mosaicity)
	{
		this->mosaicity = mosaicity;
	}

	double getSpotSize() const
	{
		return spotSize;
	}

	void setSpotSize(double spotSize)
	{
		this->spotSize = spotSize;
	}

	double getWavelength()
	{
		if (wavelength == 0)
		{
			getWavelengthFromHDF5();
		}
		return wavelength;
	}

	void setWavelength(double wavelength)
	{
		this->wavelength = wavelength;
	}

	double getRefCorrelation() const
	{
		return refCorrelation;
	}

	void setRefCorrelation(double refCorrelation)
	{
		this->refCorrelation = refCorrelation;
	}
    
    double getRefPartCorrel() const
    {
        return refPartCorrel;
    }
    
    void setRefPartCorrel(double refPartCorrel)
    {
        this->refPartCorrel = refPartCorrel;
    }

	double getExponent() const
	{
		return exponent;
	}

	void setExponent(double exponent)
	{
		this->exponent = exponent;
	}

	ScoreType getScoreType()
	{
		return scoreType;
	}

	void setScoreType(ScoreType newScoreType)
	{
		scoreType = newScoreType;
	}

	TrustLevel getTrust() const
	{
		return trust;
	}

	void setTrust(TrustLevel trust)
	{
		this->trust = trust;
	}

	MatrixPtr getMatrix()
	{
		return matrix;
	}

	static MtzManager*& getReferenceManager()
	{
        return referenceManager;
	}

	bool isRejected()
	{
		if (!reflectionCount())
		{
			return true;
		}

		return rejected;
	}

	void setRejected(bool rejected)
	{
		this->rejected = rejected;
	}

	double getScale() const
	{
		return scale;
	}

	void setScale(double scale)
	{
		this->scale = scale;
	}
    
    void addParameters(RefinementStrategyPtr map);
    static double refineParameterScore(void *object);
    
    static double getSpotSizeStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getSpotSize();
    }
    
    static void setSpotSizeStatic(void *object, double newSize)
    {
        static_cast<MtzManager *>(object)->setSpotSize(newSize);
    }

    void refineParametersStepSearch(RefinementStepSearch &map);
    
    static double getMosaicityStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getMosaicity();
    }
    
    static void setMosaicityStatic(void *object, double newMosaicity)
    {
        static_cast<MtzManager *>(object)->setMosaicity(newMosaicity);
    }
    
    void updateLatestMatrix()
    {
        if (matrix->isComplex())
            this->matrix->changeOrientationMatrixDimensions(getUnitCell());
    
        rotatedMatrix = matrix->copy();
        rotatedMatrix->rotate(hRot * M_PI / 180, kRot * M_PI / 180, 0);
    }

	static double getScaleStatic(void *object)
	{
		return static_cast<MtzManager *>(object)->scale;
	}

	static void setScaleStatic(void *object, double newScale)
	{
		static_cast<MtzManager *>(object)->setScale(newScale);
	}

    static double getHRot(void *object)
    {
        return static_cast<MtzManager *>(object)->hRot;
    }
    
    static void setHRot(void *object, double newHRot)
    {
        static_cast<MtzManager *>(object)->hRot = newHRot;
    }
    
    static double getKRot(void *object)
    {
        return static_cast<MtzManager *>(object)->kRot;
    }
    
    static void setKRot(void *object, double newKRot)
    {
        static_cast<MtzManager *>(object)->kRot = newKRot;
    }

	static double getLRot(void *object)
	{
		return static_cast<MtzManager *>(object)->lRot;
	}

	static void setLRot(void *object, double newLRot)
	{
		static_cast<MtzManager *>(object)->lRot = newLRot;
	}

    static double getWavelengthStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getWavelength();
    }

    static void setWavelengthStatic(void *object, double newWavelength)
    {
        static_cast<MtzManager *>(object)->setWavelength(newWavelength);
    }
    
    static double getBandwidthStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getBandwidth();
    }
    
    static void setBandwidthStatic(void *object, double newBandwidth)
    {
        static_cast<MtzManager *>(object)->setBandwidth(newBandwidth);
    }
    
    static double getExponentStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getExponent();
    }
    
    static void setExponentStatic(void *object, double newExponent)
    {
        static_cast<MtzManager *>(object)->setExponent(newExponent);
    }


	static double getUnitCellAStatic(void *object)
	{
		MtzManager *mtz = static_cast<MtzManager *>(object);
		double unitCellVal = mtz->getUnitCell()[0];
		return unitCellVal;
	}

	static void setUnitCellAStatic(void *object, double newDim)
	{
		static_cast<MtzManager *>(object)->_unitCell[0] = newDim;
	}

	static double getUnitCellBStatic(void *object)
	{
		MtzManager *mtz = static_cast<MtzManager *>(object);
		double unitCellVal = mtz->getUnitCell()[1];
		return unitCellVal;
	}

	static void setUnitCellBStatic(void *object, double newDim)
	{
		static_cast<MtzManager *>(object)->_unitCell[1] = newDim;
	}

	static double getUnitCellCStatic(void *object)
	{
		MtzManager *mtz = static_cast<MtzManager *>(object);
		double unitCellVal = mtz->getUnitCell()[2];
		return unitCellVal;
	}

	static void setUnitCellCStatic(void *object, double newDim)
	{
		static_cast<MtzManager *>(object)->_unitCell[2] = newDim;
	}

    ImagePtr getImagePtr()
    {
        return image.lock();
    }
    
    void setImage(ImageWeakPtr imagePtr)
    {
        image = imagePtr;
    }
    
    void setTimeDelay(double delay)
    {
        timeDelay = delay;
    }
    
    double getTimeDelay()
    {
        return timeDelay;
    }
    
    void setParentImageName(std::string pImage)
    {
        parentImage = pImage;
	}

	IndexingSolutionPtr getIndexingSolution()
	{
		return indexingSolution;
	}

	void setIndexingSolution(IndexingSolutionPtr solution)
	{
		indexingSolution = solution;
	}

	double getSearchSize()
	{
		return searchSize;
	}

	void setRefineOrientations(bool newValue)
	{
		refineOrientations = newValue;
	}
};

#endif

