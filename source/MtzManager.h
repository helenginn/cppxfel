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

using namespace CMtz;
using namespace CSym;

class Miller;
class Reflection;
class Matrix;

typedef enum
{
	ScoreTypeMinimizeRSplit = 0,
	ScoreTypeCorrelation = 1,
	ScoreTypePartialityCorrelation = 2,
	ScoreTypeMinimizeRSplitLog = 3,
	ScoreTypeCorrelationLog = 4,
	ScoreTypeAngle = 5,
	ScoreTypePartialityLeastSquares = 6,
	ScoreTypePartialityGradient = 7,
    ScoreTypeSymmetry = 8,
    ScoreTypeStandardDeviation = 9,
    ScoreTypeMinimizeRMeas = 10,
    ScoreTypeRSplitIntensity = 11,
} ScoreType;

typedef enum
{
	TrustLevelGood, TrustLevelAverage, TrustLevelBad
} TrustLevel;

class Miller;

class MtzManager : public boost::enable_shared_from_this<MtzManager>
{

protected:
    std::string filename;
	static CCP4SPG *low_group;
    static std::mutex spaceGroupMutex;
	static bool reflection_comparison(ReflectionPtr i, ReflectionPtr j);

	double extreme_index(MTZ *mtz, int max);
	void hkls_for_reflection(MTZ *mtz, float *adata, int *h, int *k, int *l,
			int *multiplier, int *offset);
	int index_for_reflection(int h, int k, int l, bool inverted);
	void findMultiplier(MTZ *mtz, int *multiplier, int *offset);
    ImageWeakPtr image;
    bool dropped;
    void getWavelengthFromHDF5();
    
    BeamPtr beam;
    
	// minimisation stuff

    int millerCount();
    
	vector<ReflectionPtr> reflections;
	vector<Reflection *> refReflections;
	vector<Reflection *> matchReflections;
    MtzManager *previousReference;
    int previousAmbiguity;
    bool allowTrust;
    bool setInitialValues;
    double penaltyWeight;
    double penaltyResolution;
    
    RotationMode rotationMode;
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
	MtzPtr copy();
	void loadParametersMap();

    void addMiller(MillerPtr miller);
    void addReflections(vector<ReflectionPtr>reflections);
	void clearReflections();
	void addReflection(ReflectionPtr reflection);
	void removeReflection(int i);
	void excludeFromLogCorrelation();
	void excludePartialityOutliers();  // delete

    MatrixPtr matrix;
    bool checkUnitCell(double trueA, double trueB, double trueC, double tolerance);
    
	void setFilename(std::string name);  // classify
    std::string getFilename(void);  // classify
	void description(void);
    void resetDefaultParameters();
    void chooseAppropriateTarget();
    
    void setDefaultMatrix();
	void setMatrix(double *components);
	void setMatrix(MatrixPtr newMat);
	void insertionSortReflections(void);
	void sortLastReflection(void);
    void applyUnrefinedPartiality();
    void incrementActiveAmbiguity();
    void cutToResolutionWithSigma(double acceptableSigma);
    double maxResolution();

    std::string filenameRoot();
	void setSpaceGroup(int spgnum);
	virtual void loadReflections(PartialityModel model, bool special = false); // FIXME: remove unnecessary parameters
    void dropReflections();
	void loadReflections(int partiality);
	static void setReference(MtzManager *reference);
	int findReflectionWithId(long unsigned int refl_id, ReflectionPtr *reflection, bool insertionPoint = false);
	void findCommonReflections(MtzManager *other,
			vector<ReflectionPtr> &reflectionVector1, vector<ReflectionPtr> &reflectionVector2,
			int *num = NULL, bool acceptableOnly = false);
	double gradientAgainstManager(MtzManager *otherManager, bool withCutoff = true, double lowRes = 0, double highRes = 0);
	void bFactorAndScale(double *scale, double *bFactor, double exponent = 1, vector<std::pair<double, double> > *dataPoints = NULL);
	void applyBFactor(double bFactor);
	void applyScaleFactor(double scaleFactor, double lowRes = 0, double highRes = 0, bool absolute = false);
	void applyScaleFactorsForBins(int binCount = 20);
	void clearScaleFactor();
	void makeScalesPermanent();
	double averageIntensity(void);
	void setSigmaToUnity();
	void setPartialityToUnity();
	double partialityRatio(ReflectionPtr imgReflection, ReflectionPtr refReflection);
	void reallowPartialityOutliers();
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

	virtual void writeToFile(std::string newFilename, bool announce = false);
    void writeToHdf5();
    void writeToDat(std::string prefix = "");
    void sendLog(LogLevel priority = LogLevelNormal);

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

    
    int symmetryRelatedReflectionCount();
    
// minimisation stuff

	int refinedParameterCount();
    
    void recalculateWavelengths();
	void getParams(double *parameters[], int paramCount = PARAM_NUM);
    void getParamPointers(double ***parameters, int paramCount = PARAM_NUM);
    void getSteps(double *parameters[], int paramCount = PARAM_NUM);
    void setParams(double parameters[], int paramCount = PARAM_NUM);

    double medianWavelength(double lowRes, double highRes);
    double wavelengthFromTopReflections(double lowRes, double highRes, int refNum);
	double bestWavelength(double lowRes = 0.0, double highRes = 0, bool usingReference = false);
	int accepted(void);
	static double exclusionScoreWrapper(void *object, double lowRes = 0,
			double highRes = 0);
    static double bFactorScoreWrapper(void *object);
    static double scoreNelderMead(void *object);
	double exclusionScore(double lowRes, double highRes, ScoreType scoreType);
    double leastSquaresPartiality(ScoreType typeOfScore);
    double leastSquaresPartiality(double low, double high, ScoreType typeOfScore = ScoreTypePartialityCorrelation);
	double correlation(bool silent = true, double lowResolution = 0, double highResolution = -1);
	double rSplit(double low, double high, bool withCutoff = false, bool set = false);
	std::string describeScoreType();
    double refinePartialitiesOrientation(int ambiguity, int cycles = 30);
    
    void refinePartialities();

    void refreshCurrentPartialities();
	void refreshPartialities(double hRot, double kRot, double mosaicity,
                             double spotSize, double wavelength, double bandwidth, double exponent,
                             double a, double b, double c);
	void refreshPartialities(double parameters[]);
	double minimizeParameter(double *meanStep, double **params, int paramNum,
			double (*score)(void *object, double lowRes, double highRes),
			void *object, double lowRes, double highRes);
	
// more grid search

    void setParamLine(std::string line);
    std::string getParamLine();
    void findSteps(int param1, int param2, std::string csvName);
	void gridSearch(bool silent = false);
	double minimize();
	double minimize(double (*score)(void *object, double lowRes, double highRes),
			void *object);
	double minimizeTwoParameters(double *meanStep1, double *meanStep2,
			double **params, int paramNum1, int paramNum2,
			double (*score)(void *object, double lowRes, double highRes),
			void *object, double lowRes, double highRes, double low);

    static void makeSuperGaussianLookupTable(double exponent);
    
    static std::string parameterHeaders();
    std::string writeParameterSummary();
    static MtzPtr trajectoryMtz(MtzPtr one, MtzPtr two, double fractionExcited);
    
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
    
    MatrixPtr getLatestMatrix();
    void addParameters(GetterSetterMapPtr map);
    static double refineParameterScore(void *object);
    
    static double getSpotSizeStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getSpotSize();
    }
    
    static void setSpotSizeStatic(void *object, double newSize)
    {
        static_cast<MtzManager *>(object)->setSpotSize(newSize);
    }

    void refineParametersStepSearch(GetterSetterMap &map);
    
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
 //       static_cast<MtzManager *>(object)->updateLatestMatrix();
    }
    
    static double getKRotStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getKRot();
    }
    
    static void setKRotStatic(void *object, double newKRot)
    {
        static_cast<MtzManager *>(object)->setKRot(newKRot);
  //      static_cast<MtzManager *>(object)->updateLatestMatrix();
    }

    static double getWavelengthStatic(void *object)
    {
        return static_cast<MtzManager *>(object)->getWavelength();
    }

    static void setWavelengthStatic(void *object, double newWavelength)
    {
        static_cast<MtzManager *>(object)->setWavelength(newWavelength);
    }

    ImagePtr getImagePtr()
    {
        return image.lock();
    }
    
    void setImage(ImageWeakPtr imagePtr)
    {
        image = imagePtr;
    }
};

#endif

