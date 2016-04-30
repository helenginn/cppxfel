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
#include <cctbx/miller/asu.h>

using cctbx::sgtbx::reciprocal_space::asu;

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
    static asu p1_asu;
    static cctbx::sgtbx::space_group p1_spg;
    static bool initialised_p1;
    static bool normalised;
    static bool correctingPolarisation;
    static bool correctedPartiality;
    static double polarisationFactor;
    static int maxSlices;
    static short int slices;
    static float trickyRes;
    static bool setupStatic;
    static int peakSize;
    
    short int h;
    short int k;
    short int l;
    bool free;
    float phase;
    char fakeFriedel; // delete?
	RlpModel rlpModel;
    double polarisationCorrection;
    double getPolarisationCorrection();
	std::map<RejectReason, bool> rejectedReasons;
	double partialCutoff; // could/should be a float
	float bFactor;
	float scale; // should be extracted from Mtz.
	float lastX;
	float lastY;
    
    double latestHRot;
    double latestKRot;
    float bFactorScale; // should be extracted from Mtz.
    bool excluded; // don't think this is used anymore
    bool rejected;
    bool calculatedRejected;
    
	std::pair<double, double> shift;
    float resol;
 
    double superGaussian(double bandwidth, double mean,
                        double sigma, double exponent);
    double integrate_beam_slice(double pBandwidth, double qBandwidth, double mean,
                               double sigma, double exponent);
    double integrate_special_beam_slice(double pBandwidth, double qBandwidth);
    double sliced_integral(double low_wavelength, double high_wavelength,
                          double spot_size_radius, double maxP, double maxQ, double mean, double sigma,
                          double exponent, bool binary = false, bool withBeamObject = false);
    
    double integrate_sphere_uniform(double p, double q);
    double integrate_sphere_gaussian(double p, double q);
    double integrate_sphere(double p, double q);
    
    double expectedRadius(double spotSize, double mosaicity, vec *hkl);
    
    BeamPtr beam;
    ImageWeakPtr image;
    IOMRefiner *indexer;
    ShoeboxPtr shoebox;
    unsigned char flipMatrix;
public:
    int getH();
    int getK();
    int getL();
    
    static void setupStaticVariables();
    vec hklVector(bool shouldFlip = true);
    void setFlipMatrix(int i);
    
    MatrixPtr getFlipMatrix();
    
	MatrixPtr matrix;
	Reflection *parentReflection;
    MtzManager *mtzParent;
    bool crossesBeamRoughly(MatrixPtr rotatedMatrix, double mosaicity,
                            double spotSize, double wavelength, double bandwidth);

	Miller(MtzManager *parent, int _h = 0, int _k = 0, int _l = 0);
	MillerPtr copy(void);
	void printHkl(void);
	static double scaleForScaleAndBFactor(double scaleFactor, double bFactor, double resol, double exponent_exponent = 1);
    void limitingEwaldWavelengths(vec hkl, double mosaicity, double spotSize, double wavelength, double *limitLow, double *limitHigh, vec *inwards = NULL, vec *outwards = NULL);
    double slicedIntegralWithVectors(vec low_wl_pos, vec high_wl_pos, double rlpSize, double mean, double sigma, double exponent);
    
    bool isOverlappedWithSpots(std::vector<SpotPtr> *spots, bool actuallyDelete = true);
    double calculateDefaultNorm();
    void setPartialityModel(PartialityModel model);
	void setData(double _intensity, double _sigma, double _partiality,
			double _wavelength);
	void setParent(Reflection *reflection);
	bool positiveFriedel(bool *positive, int *isym = NULL);
	void setRejected(RejectReason reason, bool rejection);
	bool isRejected(RejectReason reason);
	void makeScalesPermanent();
    void integrateIntensity(MatrixPtr transformedMatrix);
    vec getTransformedHKL(double hRot, double kRot);
    double getEwaldWeight(double hRot, double kRot, bool isH);
    
	bool accepted(void);
	bool isFree()
    {
        return free;
    }
    
    void setFree(bool newFree)
    {
        free = newFree;
    }
    
	void flip(void);

    bool isRejected();
    double getBFactorScale();
	double intensity(bool withCutoff = true);
	double getSigma(void);
	double getPartiality();
	double getWavelength(void);
	double getWavelength(double hRot, double kRot);
    double getWavelength(MatrixPtr transformedMatrix);
	double getWeight(bool cutoff = true, WeightType weighting = WeightTypePartialitySigma);
	double resolution();
    double twoTheta(bool horizontal);
	double scatteringAngle(ImagePtr image);

    void incrementOverlapMask(double hRot = 0, double kRot = 0);
    bool isOverlapped();
	void positionOnDetector(MatrixPtr transformedMatrix, int *x,
			int *y);
    void recalculateBetterPartiality();
    
    void setHorizontalPolarisationFactor(double newFactor);
	void recalculatePartiality(MatrixPtr rotatedMatrix, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent, bool binary = false);
	double partialityForHKL(vec hkl, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent, bool binary = false);
	void applyScaleFactor(double scaleFactor);
	double calculateNormPartiality(MatrixPtr rotatedMatrix, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent);
	double calculateNormFromResolution(MatrixPtr rotatedMatrix, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent,
			double d);
    double observedPartiality(double reference);
    double observedPartiality(MtzManager *reference);
    
    vec getTransformedHKL(MatrixPtr matrix);
    void makeComplexShoebox(double wavelength, double bandwidth, double mosaicity, double rlpSize);
    
    static double averageRawIntensity(vector<MillerPtr> millers);
    RejectReason getRejectedReason();

    virtual ~Miller();
    
    void setBeam(BeamPtr newBeam)
    {
        beam = newBeam;
    }
    
    void setExcluded(bool exc = true)
    {
        excluded = exc;
    }
    
    bool isExcluded()
    {
        return excluded;
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
    
	const std::string& getFilename() const
	{
		return filename;
	}

	void setFilename(const std::string& filename)
	{
		this->filename = filename;
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

	void setRejected(bool rejected)
	{
        setRejected(RejectReasonMerge, rejected);
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
    
    std::pair<double, double> position()
    {
        return std::make_pair(lastX, lastY);
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

	std::pair<double, double>& getShift()
	{
		return shift;
	}

	void setShift(const std::pair<int, int>& shift)
	{
		this->shift = shift;
	}
    
    void setImageAndIOMRefiner(ImagePtr newImage, IOMRefiner *indexer)
    {
        this->image = newImage;
        this->indexer = indexer;
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

    static void rotateMatrixABC(double aRot, double bRot, double cRot, MatrixPtr oldMatrix, MatrixPtr *newMatrix);
    static void rotateMatrixHKL(double hRot, double kRot, double lRot, MatrixPtr oldMatrix, MatrixPtr *newMatrix);

protected:
	static PartialityModel model;
	double rawIntensity;
	double sigma;
	double countingSigma;
	double partiality;
	double wavelength;
	std::string filename;

};

#endif /* MILLER_H_ */
