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
#include "MtzManager.h"
#include "Image.h"
#include "definitions.h"
#include "parameters.h"
#include <deque>

class Holder;
class Image;

#include "Holder.h"

typedef enum
{
	CalculationTypeOriginal, CalculationTypeIntegrate
} CalculationType;

class Miller
{
private:
    double h;
    double k;
    double l;
    short int fakeFriedel;
	bool normalised;
    bool correctingPolarisation;
	double polarisationCorrection;
    double polarisationFactor;
	double getPolarisationCorrection();
	std::map<std::string, bool> rejectedReasons;
	double partialCutoff;
	double bFactor;
	double scale;
    double gainScale;
	double lastX;
	double lastY;
    double lastBandwidth;
    double lastWavelength;
    double lastRlpSize;
    double lastMosaicity;
    double latestHRot;
    double latestKRot;
    double bFactorScale;
    bool rejected;
    bool calculatedRejected;
    
	std::pair<double, double> shift;
    double resol;
 
    double superGaussian(double bandwidth, double mean,
                        double sigma, double exponent);
    double integrate_beam_slice(double pBandwidth, double qBandwidth, double mean,
                               double sigma, double exponent);
    double sliced_integral(double low_wavelength, double high_wavelength,
                          double spot_size_radius, double maxP, double maxQ, double mean, double sigma,
                          double exponent);
    
    double expectedRadius(double spotSize, double mosaicity, vec *hkl);
    
    Image *image;
    Indexer *indexer;
    ShoeboxPtr shoebox;
    std::weak_ptr<Miller> selfPtr;
    MatrixPtr flipMatrix;
public:
    int getH();
    int getK();
    int getL();
    
    vec hklVector(bool shouldFlip = true);
    void setFlipMatrix(MatrixPtr flipMat);
	MatrixPtr matrix;
	Holder *parentHolder;
    MtzManager *mtzParent;

	Miller(MtzManager *parent, int _h = 0, int _k = 0, int _l = 0);
	MillerPtr copy(void);
	void printHkl(void);
	static double scaleForScaleAndBFactor(double scaleFactor, double bFactor, double resol, double exponent_exponent = 1);

	void setPartialityModel(PartialityModel model);
	void setData(double _intensity, double _sigma, double _partiality,
			double _wavelength);
	void setParent(Holder *holder);
	void setFree(bool newFree);
	bool positiveFriedel(bool *positive);
	bool activeWithAngle(double degrees);
	void setRejected(std::string reason, bool rejection);
	bool isRejected(std::string reason);
	void makeScalesPermanent();
    void integrateIntensity(double hRot, double kRot);

	bool accepted(void);
	bool free(void);
	void flip(void);

    bool isRejected();
    double getBFactorScale();
	double intensity(void);
	double getSigma(void);
	double getPartiality(void);
	double getWavelength(void);
	double getWavelength(double hRot, double kRot);
	double getWeight(WeightType weighting = WeightTypePartialitySigma);
	double resolution();
    double twoTheta(bool horizontal);
	double scatteringAngle(Image *image);

    void incrementOverlapMask(double hRot = 0, double kRot = 0);
    bool isOverlapped();
	void positionOnDetector(double hRot, double kRot, int *x,
			int *y);
    void calculatePosition(double distance, double wavelength, double beamX, double beamY, double mmPerPixel, double hRot, double kRot, double *x, double *y);

	void setHorizontalPolarisationFactor(double newFactor);
	void recalculatePartiality(double hRot, double kRot, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent,
			bool normalise = true);
	double partialityForHKL(vec hkl, double hRot, double kRot, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent);
	void applyScaleFactor(double scaleFactor);
	double calculateNormPartiality(double hRot, double kRot, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent);
	double calculateNormFromResolution(double hRot, double kRot, double mosaicity,
			double spotSize, double wavelength, double bandwidth, double exponent,
			double d);
    double observedPartiality(double reference);
    double observedPartiality(MtzManager *reference);
    
    void makeComplexShoebox(double wavelength, double bandwidth, double mosaicity, double rlpSize);
    
    static double averageRawIntensity(std::vector<MillerPtr> millers);

	virtual ~Miller();
    
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

	double getRawIntensity() const
	{
		return rawIntensity * scale * gainScale;
	}

	void setRawIntensity(double rawIntensity)
	{
		this->rawIntensity = rawIntensity;
	}

	void setSigma(double sigma)
	{
		this->sigma = sigma;
	}
    
    void setGainScale(double newScale)
    {
        if (newScale == newScale)
            gainScale = newScale;
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
		return countingSigma * gainScale * scale;
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

	double getLastWavelength() const
	{
		return lastWavelength;
	}

	double getLastBandwidth()
    {
        return lastBandwidth;
    }
    
    double getLastMosaicity()
    {
        return lastMosaicity;
    }
    
    double getLastRlpSize()
    {
        return lastRlpSize;
    }

	void setRejected(bool rejected)
	{
        setRejected("merge", rejected);        
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
        
    /*    if ((h == 3 && k == -4 && l == 12) || (h == -4 && k == 3 && l == -12))
        {
            std::cout << "Changing (" << h << ", " << k << ", " << l << ") to " << scale << std::endl;
        }*/
	}

	double getResolution() const
	{
		return resol;
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
    
    void setImageAndIndexer(Image *newImage, Indexer *indexer)
    {
        this->image = newImage;
        this->indexer = indexer;
    }
    
    Image *getImage()
    {
        return image;
    }
    
    void setSelf(MillerPtr ptr)
    {
        selfPtr = ptr;
    }
    
    void setCorrectingPolarisation(bool on)
    {
        correctingPolarisation = on;
    }

protected:
	PartialityModel model;
	double rawIntensity;
	double sigma;
	double countingSigma;
	double partiality;
	double wavelength;
	std::string filename;

	void rotateMatrix(double hRot, double kRot, MatrixPtr *newMatrix);
};

#endif /* MILLER_H_ */
