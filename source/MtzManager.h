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

using namespace CMtz;
using namespace CSym;

typedef enum
{
	ScoreTypeMinimizeRSplit = 0,
	ScoreTypeCorrelation = 1,
    ScoreTypeSymmetry = 8,
    ScoreTypeMinimizeRMeas = 10,
    ScoreTypeRSplitIntensity = 11,
} ScoreType;

typedef enum
{
	TrustLevelGood, TrustLevelAverage, TrustLevelBad
} TrustLevel;

class Miller;

class MtzManager : public LoggableObject, public hasFilename, public boost::enable_shared_from_this<MtzManager>
{

protected:
    std::string parentImage;
	static CCP4SPG *low_group;
    static std::mutex spaceGroupMutex;
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
	double mosaicity;
	double spotSize;
	double wavelength;
	double bandwidth;
	double exponent;
	double refCorrelation;
    double refPartCorrel;
	double scale;
	double cellDim[3];
	double cellAngles[3];
    
    double timeDelay;
    
    uint32 activeAmbiguity;
    
	bool usingFixedWavelength;

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
	double maxResolutionAll;
    float lastRSplit;

    static double superGaussianScale;
    double *params;
    double detectorDistance;
    
	bool finalised;
	bool rejected;

	TrustLevel trust;
	ScoreType scoreType;

    MatrixPtr baseMatrix;
    MatrixPtr rotatedMatrix;
	static MtzManager *referenceManager;
	MtzManager *lastReference;

    std::ostringstream logged;
public:
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
    void addReflections(vector<ReflectionPtr>reflections, bool assumeSorted = false);
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

	void setSpaceGroup(int spgnum);
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
	void applyScaleFactorsForBins(int binCount = 20);
	void clearScaleFactor();
	double averageIntensity();
	void setSigmaToUnity();
	void replaceBeamWithSpectrum();

	void setUnitCell(double a, double b, double c, double alpha, double beta,
			double gamma);
	void setUnitCell(vector<double> unitCell);
	void getUnitCell(double *a, double *b, double *c, double *alpha, double *beta,
			double *gamma);
    std::vector<double> getUnitCell()
    {
        std::vector<double> cell;
        cell.push_back(cellDim[0]);
        cell.push_back(cellDim[1]);
        cell.push_back(cellDim[2]);
        cell.push_back(cellAngles[0]);
        cell.push_back(cellAngles[1]);
        cell.push_back(cellAngles[2]);

        return cell;
    }
    
	void copySymmetryInformationFromManager(MtzPtr toCopy);
	void applyPolarisation(void);

	virtual void writeToFile(std::string newFilename, bool announce = false, bool plusAmbiguity = false);
    void writeToHdf5();
    void writeToDat(std::string prefix = "");
    
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

	int refinedParameterCount();
    
    void recalculateWavelengths();

    double medianWavelength(double lowRes, double highRes);
    double bestWavelength(double lowRes = 0.0, double highRes = 0, bool usingReference = false);
	int accepted(void);
	static double exclusionScoreWrapper(void *object, double lowRes = 0,
			double highRes = 0);
    double leastSquaresPartiality(double low = 0, double high = 0);
    double correlation(bool silent = true, double lowResolution = 0, double highResolution = -1);
	double rSplit(double low, double high);
	std::string describeScoreType();
    double refinePartialitiesOrientation(int ambiguity, bool reset = true);
    
    void refinePartialities();

    void refreshCurrentPartialities();
    // delete
	void refreshPartialities(double hRot, double kRot, double mosaicity,
                             double spotSize, double wavelength, double bandwidth, double exponent,
                             double a, double b, double c);
	void refreshPartialities(double parameters[]);
	
// more grid search

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

    double getDetectorDistance()
    {
        return detectorDistance;
    }
    
    void setDetectorDistance(double distance)
    {
        detectorDistance = distance;
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
    
    double getSuperGaussianScale() const
    {
        return superGaussianScale;
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

	double getWavelength() const
	{
		return wavelength;
	}

	void setWavelength(double wavelength)
	{
		this->wavelength = wavelength;
	}

	double getHRot() const
	{
		return hRot;
	}

	void setHRot(double hRot)
	{
		this->hRot = hRot;
	}

	double getKRot() const
	{
		return kRot;
	}

	void setKRot(double kRot)
	{
		this->kRot = kRot;
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

	CCP4SPG*& getLowGroup()
	{
		return low_group;
	}

	void setLowGroup(CCP4SPG*& lowGroup)
	{
		low_group = lowGroup;
	}

	TrustLevel getTrust() const
	{
		return trust;
	}

	void setTrust(TrustLevel trust)
	{
		this->trust = trust;
	}

	bool isFinalised() const
	{
		return finalised;
	}

	void setFinalised(bool finalised)
	{
		this->finalised = finalised;
	}

	MatrixPtr getMatrix()
	{
		return matrix;
	}

	static MtzManager*& getReferenceManager()
	{
        return referenceManager;
	}

	bool isUsingFixedWavelength() const
	{
		return usingFixedWavelength;
	}

	void setUsingFixedWavelength(bool usingFixedWavelength)
	{
		this->usingFixedWavelength = usingFixedWavelength;
	}

	bool isOptimisingBandwidth() const
	{
		return optimisingBandwidth;
	}

	void setOptimisingBandwidth(bool optimisingBandwidth)
	{
		this->optimisingBandwidth = optimisingBandwidth;
	}

	bool isOptimisingExponent() const
	{
		return optimisingExponent;
	}

	void setOptimisingExponent(bool optimisingExponent)
	{
		this->optimisingExponent = optimisingExponent;
	}

	bool isOptimisingMosaicity() const
	{
		return optimisingMosaicity;
	}

	void setOptimisingMosaicity(bool optimisingMosaicity)
	{
		this->optimisingMosaicity = optimisingMosaicity;
	}

	bool isOptimisingOrientation() const
	{
		return optimisingOrientation;
	}

	void setOptimisingOrientation(bool optimisingOrientation)
	{
		this->optimisingOrientation = optimisingOrientation;
	}

	bool isOptimisingRlpSize() const
	{
		return optimisingRlpSize;
	}

	void setOptimisingRlpSize(bool optimisingRlpSize)
	{
		this->optimisingRlpSize = optimisingRlpSize;
	}

	bool isOptimisingWavelength() const
	{
		return optimisingWavelength;
	}

	void setOptimisingWavelength(bool optimisingWavelength)
	{
		this->optimisingWavelength = optimisingWavelength;
	}

	double getStepSizeBandwidth() const
	{
		return stepSizeBandwidth;
	}

	void setStepSizeBandwidth(double stepSizeBandwidth)
	{
		this->stepSizeBandwidth = stepSizeBandwidth;
	}

	double getStepSizeExponent() const
	{
		return stepSizeExponent;
	}

	void setStepSizeExponent(double stepSizeExponent)
	{
		this->stepSizeExponent = stepSizeExponent;
	}

	double getStepSizeMosaicity() const
	{
		return stepSizeMosaicity;
	}

	void setStepSizeMosaicity(double stepSizeMosaicity)
	{
		this->stepSizeMosaicity = stepSizeMosaicity;
	}

	double getStepSizeOrientation() const
	{
		return stepSizeOrientation;
	}

	void setStepSizeOrientation(double stepSizeOrientation)
	{
		this->stepSizeOrientation = stepSizeOrientation;
	}

	double getStepSizeRlpSize() const
	{
		return stepSizeRlpSize;
	}

	void setStepSizeRlpSize(double stepSizeRlpSize)
	{
		this->stepSizeRlpSize = stepSizeRlpSize;
	}

	double getStepSizeWavelength() const
	{
		return stepSizeWavelength;
	}

	void setStepSizeWavelength(double stepSizeWavelength)
	{
		this->stepSizeWavelength = stepSizeWavelength;
	}

	bool isRejected() const
	{
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
            this->matrix->changeOrientationMatrixDimensions(cellDim[0], cellDim[1], cellDim[2], cellAngles[0], cellAngles[1], cellAngles[2]);
    
        rotatedMatrix = matrix->copy();
        rotatedMatrix->rotate(hRot * M_PI / 180, kRot * M_PI / 180, 0);
    }
    
    static double getHRotStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getHRot();
    }
    
    static void setHRotStatic(void *object, double newHRot)
    {
        static_cast<MtzManager *>(object)->setHRot(newHRot);
    }
    
    static double getKRotStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getKRot();
    }
    
    static void setKRotStatic(void *object, double newKRot)
    {
        static_cast<MtzManager *>(object)->setKRot(newKRot);
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
        return static_cast<MtzManager *>(object)->getUnitCell()[0];
    }
    
    static void setUnitCellAStatic(void *object, double newDim)
    {
        static_cast<MtzManager *>(object)->getUnitCell()[0] = newDim;
    }
    
    static double getUnitCellBStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getUnitCell()[1];
    }
    
    static void setUnitCellBStatic(void *object, double newDim)
    {
        static_cast<MtzManager *>(object)->getUnitCell()[1] = newDim;
    }

    static double getUnitCellCStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getUnitCell()[2];
    }
    
    static void setUnitCellCStatic(void *object, double newDim)
    {
        static_cast<MtzManager *>(object)->getUnitCell()[2] = newDim;
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
    
    std::string getParentImageName()
    {
        return parentImage;
    }
    
    void setParentImageName(std::string pImage)
    {
        parentImage = pImage;
    }
};

#endif

