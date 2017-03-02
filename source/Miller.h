/*
 * Miller.h
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#ifndef MILLER_H_
#define MILLER_H_

#include "Matrix.h"
#include <string>

#include <map>
#include "definitions.h"
#include "parameters.h"

class MtzManager;
class Reflection;
class Image;

typedef enum
{
	CalculationTypeOriginal, CalculationTypeIntegrate
} CalculationType;

typedef enum
{
    RlpModelUniform, RlpModelGaussian,
} RlpModel;

class Miller : public LoggableObject, public boost::enable_shared_from_this<Miller>
{
private:
    static bool normalised;
    static bool correctingPolarisation;
    static bool correctedPartiality;
    static double polarisationFactor;
    static int maxSlices;
    static short int slices;
    static float trickyRes;
    static bool setupStatic;
    static int peakSize;
    static bool individualWavelength;
    
    short int h;
    short int k;
    short int l;
    bool free;
    float phase;
	RlpModel rlpModel;
    double polarisationCorrection;
    double getPolarisationCorrection();
    unsigned int rejectedReasons;
	float partialCutoff; // could/should be a float
	float bFactor;
	float scale; // should be extracted from Mtz. Maybe?
	float lastX;
	float lastY;
    float lastPeakX;
    float lastPeakY;
    float correctedX;
    float correctedY;
    double predictedWavelength;
    vec shiftedRay;
    
    float bFactorScale; // should be extracted from Mtz.
    
	std::pair<float, float> shift;
    float resol;
 
    double superGaussian(double bandwidth, double mean,
                        double sigma, double exponent);
    double integrate_beam_slice(double pBandwidth, double qBandwidth, double mean,
                               double sigma, double exponent);
    double integrate_special_beam_slice(double pBandwidth, double qBandwidth);
    double sliced_integral(double low_wavelength, double high_wavelength,
                          double spot_size_radius, double maxP, double maxQ, double mean, double sigma,
                          double exponent, bool binary = false, bool withBeamObject = false, bool fixPredicted = false);
    
    double integrate_sphere_uniform(double p, double q);
    double integrate_sphere_gaussian(double p, double q);
    double integrate_sphere(double p, double q);
    void recalculatePredictedWavelength();

    double expectedRadius(double spotSize, double mosaicity, vec *hkl);
    
    BeamPtr beam;
    ImageWeakPtr image;
    DetectorWeakPtr lastDetector;
    IOMRefinerWeakPtr refiner;
    ShoeboxPtr shoebox;
    unsigned char flipMatrix;
    static bool absoluteIntensity;
    static double intensityThreshold;
public:
    int getH();
    int getK();
    int getL();
    
    bool is(int _h, int _k, int _l)
    {
        return (h == _h && k == _k && l == _l);
    }
    
    static void setupStaticVariables();
    vec hklVector(bool shouldFlip = true);
    void setFlipMatrix(int i);
    
    MatrixPtr getFlipMatrix();
    
	MatrixPtr matrix;
	ReflectionWeakPtr parentReflection;
    MtzManager *mtzParent;
    bool crossesBeamRoughly(MatrixPtr rotatedMatrix, double mosaicity,
                            double spotSize, double wavelength, double bandwidth);

	Miller(MtzManager *parent, int _h = 0, int _k = 0, int _l = 0, bool calcFree = true);
	MillerPtr copy(void);
	static double scaleForScaleAndBFactor(double scaleFactor, double bFactor, double resol, double exponent_exponent = 1);
    void limitingEwaldWavelengths(vec hkl, double mosaicity, double spotSize, double wavelength, double *limitLow, double *limitHigh, vec *inwards = NULL, vec *outwards = NULL);
    double slicedIntegralWithVectors(vec low_wl_pos, vec high_wl_pos, double rlpSize, double mean, double sigma, double exponent);
    
    bool isOverlappedWithSpots(std::vector<SpotPtr> *spots, bool actuallyDelete = true);
    void setPartialityModel(PartialityModel model);
	void setData(double _intensity, double _sigma, double _partiality,
			double _wavelength);
	void setParent(ReflectionPtr reflection);
	bool positiveFriedel(bool *positive, int *isym = NULL);
	void setRejected(RejectReason reason, bool rejection);
	bool isRejected(RejectReason reason);
	void makeScalesPermanent();
    void integrateIntensity(MatrixPtr transformedMatrix);
    vec getTransformedHKL(double hRot, double kRot);
    
	bool accepted(void);
	bool isFree()
    {
        return free;
    }
    
    void setFree(bool newFree)
    {
        free = newFree;
    }

    bool isRejected();
    double getBFactorScale();
	double intensity(bool withCutoff = true);
	double getSigma();
	double getPartiality();
	double getWavelength();
	double getWavelength(double hRot, double kRot);
    double getWavelength(MatrixPtr transformedMatrix);
    double recalculateWavelength();
	double getWeight(bool cutoff = true, WeightType weighting = WeightTypePartialitySigma);
	double resolution();
    double twoTheta(bool horizontal);

    void incrementOverlapMask(double hRot = 0, double kRot = 0);
    bool isOverlapped();
	void positionOnDetector(MatrixPtr transformedMatrix, int *x,
			int *y, bool search = true);
    void recalculateBetterPartiality();
    
    void setHorizontalPolarisationFactor(double newFactor);
	void recalculatePartiality(MatrixPtr rotatedMatrix, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent, bool binary = false, bool no_norm = false);
	double partialityForHKL(vec hkl, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent, bool binary = false, bool fixPredicted = false);
	void applyScaleFactor(double scaleFactor);
	double calculateNormPartiality(MatrixPtr rotatedMatrix, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent);
	double calculateNormFromResolution(MatrixPtr rotatedMatrix, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent,
			double d);
    double observedPartiality(double reference);
    double observedPartiality(MtzManager *reference);
    
    static void refreshMillerPositions(std::vector<MillerWeakPtr> millers);
    static void refreshMillerPositions(std::vector<MillerPtr> millers);
    vec getTransformedHKL(MatrixPtr matrix = MatrixPtr());
    void makeComplexShoebox(double wavelength, double bandwidth, double mosaicity, double rlpSize);
    
    void setHKL(int _h, int _k, int _l)
    {
        h = _h;
        k = _k;
        l = _l;
    }
    
    DetectorPtr getDetector()
    {
        return lastDetector.lock();
    }
    
    void setDetector(DetectorPtr newD);
    
    ShoeboxPtr getShoebox()
    {
        return shoebox;
    }
    
    static double averageRawIntensity(vector<MillerPtr> millers);
    RejectReason getRejectedReason();

    virtual ~Miller();
    
    bool reachesThreshold();
    
    void setBeam(BeamPtr newBeam)
    {
        beam = newBeam;
    }
    
    MtzManager *&getMtzParent()
    {
        return mtzParent;
    }

    void setMtzParent(MtzManager *mtz)
    {
        mtzParent = mtz;
    }
    
	void setPartiality(double partiality)
	{
		this->partiality = partiality;
	}
    
    double getRawestIntensity();
    

	double getRawIntensity() const
	{
		return rawIntensity * scale;
	}

	void setRawIntensity(double rawIntensity)
	{
		this->rawIntensity = rawIntensity;
	}

	void setSigma(double sigma)
	{
		this->sigma = sigma;
	}

	void applyPolarisation(double wavelength);

	double getCountingSigma() const
	{
		return countingSigma * scale;
	}
    
    double getRawCountingSigma() const
    {
        return countingSigma;
    }

	void setCountingSigma(double countingSigma)
	{
		this->countingSigma = countingSigma;
	}

	bool isNormalised() const
	{
		return normalised;
	}

	void setNormalised(bool normalised)
	{
		this->normalised = normalised;
	}

	MatrixPtr getMatrix()
	{
		return matrix;
	}

	void setMatrix(MatrixPtr matrix)
	{
		this->matrix = matrix;
	}

	void setPolarisationCorrection(double polarisationCorrection)
	{
		this->polarisationCorrection = polarisationCorrection;
	}

	void setRejected(int rejected)
	{
        rejectedReasons = rejected;
	}
    
    Coord getLastXY()
    {
        Coord lastXY = std::make_pair(lastX, lastY);
        return lastXY;
    }

	double getLastX() const
	{
		return lastX;
	}

	double getLastY() const
	{
		return lastY;
	}
    
    double getCorrectedX() const
    {
        return correctedX;
    }

    double getCorrectedY() const
    {
        return correctedY;
    }

    std::pair<double, double> position()
    {
        return std::make_pair(lastX, lastY);
    }

    void setCorrectedX(double lastX)
    {
        this->correctedX = lastX;
    }
    
    void setCorrectedY(double lastY)
    {
        this->correctedY = lastY;
    }
    
	void setLastX(double lastX)
	{
		this->lastX = lastX;
	}

	void setLastY(double lastY)
	{
		this->lastY = lastY;
	}

	double getPartialCutoff() const
	{
		return partialCutoff;
	}

	void setPartialCutoff(double partialCutoff)
	{
		this->partialCutoff = partialCutoff;
	}

	double getBFactor() const
	{
		return bFactor;
	}

	void setBFactor(double factor)
	{
		if (factor == factor)
			bFactor = factor;
        
        bFactorScale = 0;
	}

	double getScale() const
	{
		return scale;
	}

	void setScale(double scale)
	{
		if (scale == scale)
			this->scale = scale;
	}

	double getResolution()
	{
        return resolution();
	}

	void setResolution(double resol)
	{
		this->resol = resol;
	}

	std::pair<float, float>& getShift()
	{
		return shift;
	}
    
    float *getXShiftPointer()
    {
        return &(shift.first);
    }
    
    float *getYShiftPointer()
    {
        return &(shift.second);
    }

	void setShift(const std::pair<int, int>& shift)
	{
		this->shift = shift;
	}
    
    void setImageAndIOMRefiner(ImagePtr newImage, IOMRefinerPtr indexer)
    {
        if (newImage) this->image = newImage;
        if (indexer) this->refiner = indexer;
    }
    
    IOMRefinerPtr getIOMRefiner()
    {
        return refiner.lock();
    }
    
    ImagePtr getImage()
    {
        return image.lock();
    }
    
    void setCorrectingPolarisation(bool on)
    {
        correctingPolarisation = on;
    }
    
    void setPhase(double newPhase)
    {
        phase = newPhase;
    }
    
    double getPhase()
    {
        return phase;
    }
    
    void setWavelength(double wave)
    {
        wavelength = wave;
    }
    
    int getRejectionFlags()
    {
        return rejectedReasons;
    }
    
    ReflectionPtr getParentReflection()
    {
        return parentReflection.lock();
    }
    
    double getRawSigma()
    {
        return sigma;
    }
    
    double getPredictedWavelength();
    
    vec getShiftedRay();
    
    static void rotateMatrixHKL(double hRot, double kRot, double lRot, MatrixPtr oldMatrix, MatrixPtr *newMatrix);

protected:
	static PartialityModel model;
	float rawIntensity;
	float sigma;
	float countingSigma;
	double partiality;
	float wavelength;

};

#endif /* MILLER_H_ */
