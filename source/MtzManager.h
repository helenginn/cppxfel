#ifndef mtz_manager
#define mtz_manager

#include <string>
#include <iostream>
#include "headers/cmtzlib.h"
#include "headers/csymlib.h"
#include <vector>
#include "Matrix.h"
#include "liblbfgs.h"

#include "Holder.h"
#include "Miller.h"
//#include <tr1/memory>
#include "definitions.h"
#include "parameters.h"
#include <tuple>
#include "Logger.h"

using namespace std;
using namespace CMtz;
using namespace CSym;

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
} ScoreType;

typedef enum
{
	TrustLevelGood, TrustLevelAverage, TrustLevelBad
} TrustLevel;

class Miller;
class Holder;

class MtzManager
{

protected:
	string filename;
	CCP4SPG *low_group;
	static bool holder_comparison(Holder *i, Holder *j);

	double extreme_index(MTZ *mtz, int max);
	void hkls_for_reflection(MTZ *mtz, float *adata, int *h, int *k, int *l,
			int *multiplier, int *offset);
	int index_for_reflection(int h, int k, int l, bool inverted);
	void findMultiplier(MTZ *mtz, int *multiplier, int *offset);

	// minimisation stuff

	vector<Holder *> holders;
	vector<Holder *> refHolders;
	vector<Holder *> matchHolders;
    MtzManager *previousReference;
    int previousAmbiguity;
    bool allowTrust;
    
	double hRot;
	double kRot;
	double mosaicity;
	double spotSize;
	double wavelength;
	double bandwidth;
	double exponent;
	double refCorrelation;
	double scale;
	double cellDim[3];
	double cellAngles[3];
    int activeAmbiguity;
    
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
	double stepSizeExponent;

	double toleranceWavelength;
	double toleranceBandwidth;
	double toleranceMosaicity;
	double toleranceRlpSize;
	double toleranceOrientation;
	double toleranceExponent;

	ScoreType defaultScoreType;
	bool usePartialityFunction;
	double maxResolutionAll;
	double maxResolutionRlpSize;

    double trialUnitCell;
    double superGaussianScale;
    double lastExponent;
    double *params;
    double detectorDistance;
    
	bool finalised;
	bool inverse;
	bool flipped;
	static double lowRes;
	static double highRes;
	bool freePass;
	bool rejected;

	TrustLevel trust;
	ScoreType scoreType;

    MatrixPtr baseMatrix;
	static MtzManager *referenceManager;
	MtzManager *lastReference;

    static double unitCellScore(void *object);
    ostringstream logged;
public:
    vector<double> superGaussianTable;
    double bFactor;
    int rejectOverlaps();

	MtzManager(void);
	virtual ~MtzManager(void);
	MtzManager *copy();
	void loadParametersMap();

    void addHolders(vector<Holder *>holders);
	void clearHolders();
	void addHolder(Holder *holder);
	void removeHolder(int i);
	void excludeFromLogCorrelation();
	void excludePartialityOutliers();

    MatrixPtr matrix;

	void setFilename(string name);
	string getFilename(void);
	void description(void);

    void setDefaultMatrix();
	void setMatrix(double *components);
	void setMatrix(MatrixPtr newMat);
	void insertionSortHolders(void);
	void sortLastHolder(void);
    void applyUnrefinedPartiality();
    void incrementActiveAmbiguity();
    void cutToResolutionWithSigma(double acceptableSigma);
    double maxResolution();

    std::string filenameRoot();
	void setSpaceGroup(int spgnum);
	virtual void loadReflections(PartialityModel model);
	void loadReflections(int partiality);
	static void setReference(MtzManager *reference);
	int findHolderWithId(int refl_id, Holder **holder, bool insertionPoint = false);
	void findCommonReflections(MtzManager *other,
			vector<Holder *> &holderVector1, vector<Holder *> &holderVector2,
			int *num = NULL, bool force = false);
	double meanCorrectedIntensity(Holder *holder);
	double gradientAgainstManager(MtzManager &otherManager, bool leastSquares = false, double lowRes = 0, double highRes = 0);
	double minimizeGradient(MtzManager *otherManager, bool leastSquares);
	void bFactorAndScale(double *scale, double *bFactor, double exponent = 1, vector<pair<double, double> > *dataPoints = NULL);
	double minimizeRFactor(MtzManager *otherManager);
    void applyBFactor(double bFactor);
	void applyScaleFactor(double scaleFactor, double lowRes = 0, double highRes = 0, bool absolute = false);
	void applyScaleFactorsForBins();
	void clearScaleFactor();
	void makeScalesPermanent();
	double averageIntensity(void);
	void setSigmaToUnity();
	void setPartialityToUnity();
	double partialityRatio(Holder *imgHolder, Holder *refHolder);
	void getRefHolders(vector<Holder *> *refPointer,
			vector<Holder *> *matchPointer);
	void reallowPartialityOutliers();

	void setUnitCell(double a, double b, double c, double alpha, double beta,
			double gamma);
	void setUnitCell(vector<double> unitCell);
	void getUnitCell(double *a, double *b, double *c, double *alpha, double *beta,
			double *gamma);
	void copySymmetryInformationFromManager(MtzPtr toCopy);
	void applyPolarisation(void);

	void writeToFile(string newFilename, bool announce = false, bool shifts = false, bool includeAmbiguity = false);
	void writeToDat();
    void sendLog(LogLevel priority = LogLevelNormal);

	double correlationWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false);

	double rSplitWithManager(MtzManager *otherManager, bool printHits = false,
			bool silent = true, double lowRes = 0, double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL, bool shouldLog = false);

	double statisticsWithManager(MtzManager *otherManager, StatisticsFunction *function,
			bool printHits, bool silent, double lowRes, double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog);

	double statisticsWithManager(MtzManager *otherManager,
			StatisticsFunction *function, RFactorFunction *rFactorFunction,
			RFactorType rFactor, bool printHits, bool silent, double lowRes,
			double highRes, int bins,
			vector<boost::tuple<double, double, double, int> > *correlations,
			bool shouldLog);

	double rFactorWithManager(RFactorType rFactor, bool printHits = false, bool silent = true, double lowRes = 0,
			double highRes = 0, int bins = 20,
			vector<boost::tuple<double, double, double, int> > *correlations = NULL,
			bool shouldLog = false);

    int symmetryRelatedReflectionCount();
    
// minimisation stuff

	static MtzManager *currentManager;

    int refinedParameterCount();
    
	static int paramMult, spotMult;
	void getParams(double *parameters[], int paramCount = PARAM_NUM);
	void setParams(double parameters[], int paramCount = PARAM_NUM);

	double bestWavelength(double lowRes = 0.0, double highRes = 0);
	double weightedBestWavelength(double lowRes, double highRes);
	int accepted(void);
	static double exclusionScoreWrapper(void *object, double lowRes,
			double highRes);
    static double bFactorScoreWrapper(void *object);
	double exclusionScore(double lowRes, double highRes, ScoreType scoreType);
	double leastSquaresPartiality(double low, double high, ScoreType typeOfScore = ScoreTypePartialityCorrelation);
	double correlation(bool silent = true, double lowResolution = 0, double highResolution = -1);
	double rSplit(double low, double high, bool square = false);
	string describeScoreType();
    
	void refreshPartialities(double hRot, double kRot, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent);
	void refreshPartialities(double parameters[]);
	double minimizeParameter(double *meanStep, double **params, int paramNum,
			double (*score)(void *object, double lowRes, double highRes),
			void *object, double lowRes, double highRes);
	double partialHits(vector<double> *partials, vector<double> *percentages);
	double partialityGradient();

// more grid search

    void findSteps();
	void gridSearch(bool spotSize);
	void gridSearch(double (*score)(void *object, double lowRes, double highRes),
			void *object);
	static void gridSearchWrapper(MtzManager *image, bool spotSize);
	double leastSquaresPartiality(ScoreType typeOfScore);
	double minimize(bool minimizeSpotSize = true, bool suppressOutput = true);
	double minimize(bool minimizeSpotSize, bool suppress,
			double (*score)(void *object, double lowRes, double highRes),
			void *object);
	double minimizeTwoParameters(double *meanStep1, double *meanStep2,
			double **params, int paramNum1, int paramNum2,
			double (*score)(void *object, double lowRes, double highRes),
			void *object, double lowRes, double highRes, double low);

    void makeSuperGaussianLookupTable(double exponent);
    
    int ambiguityCount();
    
    void flipToActiveAmbiguity();
    void resetFlip();
    
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

    virtual int holderCount()
    {
        return (int)holders.size();
    }
    
    virtual Holder *holder(int i)
    {
        return holders[i];
    }
    
    void setDefaultScoreType(ScoreType scoreType)
    {
        defaultScoreType = scoreType;
    }
    
    double getTrialUnitCell()
    {
        return trialUnitCell;
    }
    
    double getSuperGaussianScale() const
    {
        return superGaussianScale;
    }
    
    double getLastExponent() const
    {
        return lastExponent;
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

	static double getHighRes()
	{
		return highRes;
	}

	static void setHighRes(double _highRes)
	{
		highRes = _highRes;
	}

	static double getLowRes()
	{
		return lowRes;
	}

	static void setLowRes(double _lowRes)
	{
		lowRes = _lowRes;
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

	bool isFreePass() const
	{
		return freePass;
	}

	void setFreePass(bool freePass)
	{
		this->freePass = freePass;
	}

	bool isFlipped() const
	{
		return flipped;
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
};

#endif

