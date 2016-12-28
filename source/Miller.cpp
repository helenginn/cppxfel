/*
 * Miller.cpp
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#include "Panel.h"
#include "Miller.h"
#include "MtzManager.h"

#include "Reflection.h"
#include "definitions.h"
#include "Vector.h"
#include <cmath>
#include "parameters.h"
#include "Shoebox.h"
#include <memory>
#include <cctbx/miller/asu.h>
#include <cctbx/miller.h>
#include "FileParser.h"
#include "Beam.h"
#include "IOMRefiner.h"
#include "Image.h"
#include "Spot.h"
#include "FreeMillerLibrary.h"
#include "polyfit.hpp"
#include "Detector.h"

bool Miller::normalised = true;
bool Miller::correctingPolarisation = false;
double Miller::polarisationFactor = false;
int Miller::maxSlices = 100;
short int Miller::slices = 8;
float Miller::trickyRes = 8.0;
bool Miller::setupStatic = false;
PartialityModel Miller::model = PartialityModelScaled;
int Miller::peakSize = 0;
bool Miller::correctedPartiality = false;
bool Miller::absoluteIntensity = false;
double Miller::intensityThreshold = 0;

using cctbx::sgtbx::reciprocal_space::asu;

asu Miller::p1_asu = asu();
bool Miller::initialised_p1 = false;
space_group Miller::p1_spg = space_group();

void Miller::setupStaticVariables()
{
    if (setupStatic)
        return;
    
    model = FileParser::getKey("BINARY_PARTIALITY", false) ? PartialityModelBinary : PartialityModelScaled;
    
    normalised = FileParser::getKey("NORMALISE_PARTIALITIES", true);
    correctingPolarisation = FileParser::getKey("POLARISATION_CORRECTION", false);
    polarisationFactor = FileParser::getKey("POLARISATION_FACTOR", 0.0);
    slices = FileParser::getKey("PARTIALITY_SLICES", 8);
    trickyRes = FileParser::getKey("CAREFUL_RESOLUTION", 8.0);
    maxSlices = FileParser::getKey("MAX_SLICES", 100);
    peakSize = FileParser::getKey("FOCUS_ON_PEAK_SIZE", 0);
    correctedPartiality = FileParser::getKey("CORRECTED_PARTIALITY_MODEL", false);
    absoluteIntensity = FileParser::getKey("ABSOLUTE_INTENSITY", false);
    intensityThreshold = FileParser::getKey("INTENSITY_THRESHOLD", INTENSITY_THRESHOLD);
    
    setupStatic = true;
}

double Miller::integrate_sphere_gaussian(double p, double q)
{
    double mean = 0.5;
    double sigma = 0.25;
    double exponent = 1.5;
    
    double partiality_add = super_gaussian(q, mean, sigma, exponent);
    double partiality_remove = super_gaussian(p, mean, sigma, exponent);
    
    double proportion = (partiality_add + partiality_remove) / 2;
    proportion *= (q - p);
    
    
    return proportion;
}

double Miller::integrate_sphere_uniform(double p, double q)
{
    double partiality_add = 4 * q * (1 - q);
    double partiality_remove = 4 * p * (1 - p);
    
    double proportion = (partiality_add + partiality_remove) / 2;
    proportion *= (q - p);
    
    return proportion;
}

double Miller::integrate_sphere(double p, double q)
{
    if (rlpModel == RlpModelGaussian)
        return integrate_sphere_gaussian(p, q);
    if (rlpModel == RlpModelUniform)
        return integrate_sphere_uniform(p, q);
    
    else return integrate_sphere_uniform(p, q);
}

double Miller::superGaussian(double bandwidth, double mean,
                            double sigma, double exponent)
{
    if (mtzParent == NULL || !mtzParent->setupGaussianTable())
    {
        return super_gaussian(bandwidth, mean, sigma, exponent);
    }
    else
    {
        if (bandwidth != bandwidth || mean != mean)
            return 0;
        
        double standardisedX = fabs((bandwidth - mean) / sigma);
        
        if (standardisedX > MAX_SUPER_GAUSSIAN)
            return 0;
        
        if (!std::isfinite(standardisedX) || standardisedX != standardisedX)
            return 0;
        
        const double step = SUPER_GAUSSIAN_STEP;
        
        int lookupInt = standardisedX / step;
        return mtzParent->superGaussianTable[lookupInt];
    }
}

double Miller::integrate_special_beam_slice(double pBandwidth, double qBandwidth)
{
    return beam->integralBetweenEwaldWavelengths(pBandwidth, qBandwidth);
}

double Miller::integrate_beam_slice(double pBandwidth, double qBandwidth, double mean,
                                   double sigma, double exponent)
{
    if (mean != mean || sigma != sigma)
        return 0;
    
    double pValue = superGaussian(pBandwidth, mean, sigma, exponent);
    double qValue = superGaussian(qBandwidth, mean, sigma, exponent);
    
    double width = fabs(pBandwidth - qBandwidth) / sigma;
    
    double area = (pValue + qValue) / 2 * width;
    
    return area;
}

double integrate_beam(double pBandwidth, double qBandwidth, double mean,
                     double sigma)
{
    double normal_add = cdf(pBandwidth, mean, sigma);
    double normal_remove = cdf(qBandwidth, mean, sigma);
    
    double normal_slice = (normal_add - normal_remove);
    
    return normal_slice;
}

double Miller::sliced_integral(double low_wavelength, double high_wavelength,
                              double spot_size_radius, double maxP, double maxQ, double mean, double sigma,
                              double exponent, bool binary, bool withBeamObject)
{
    double mySlices = slices;
    
    if (resolution() < 1 / trickyRes)
    {
        mySlices = maxSlices;
    }
    
    double bandwidth_span = high_wavelength - low_wavelength;
    
    double fraction_total = 1;
    double currentP = 0;
    double currentQ = (double)1 / (double)mySlices;
    
    double total_integral = 0;
    
    double total_normal = 0;
    double total_sphere = 0;
    
    for (int i = 0; i < mySlices; i++)
    {
        double pBandwidth = high_wavelength - bandwidth_span * currentP;
        double qBandwidth = high_wavelength - bandwidth_span * currentQ;
        
        double normalSlice = 0;
        
        if (!withBeamObject)
        {
            normalSlice = integrate_beam_slice(pBandwidth, qBandwidth, mean,
                                                 sigma, exponent);
        }
        else
        {
            normalSlice = integrate_special_beam_slice(qBandwidth, pBandwidth);
        }
        
        double sphereSlice = 0;
        
        if (normalSlice > 0 && !binary)
        {
            sphereSlice = integrate_sphere(currentP, currentQ);
        }
        
        double totalSlice = normalSlice * sphereSlice;
        
        if (binary && normalSlice > 0.01)
            return 1.0;
        
        total_integral += totalSlice;
        total_normal += normalSlice;
        total_sphere += sphereSlice;
        
        currentP = currentQ;
        currentQ += fraction_total / mySlices;
    }
    
    return total_integral;
}

vec incrementedPosition(vec low_wl_pos, vec high_wl_pos, double proportion)
{
    vec diff = copy_vector(high_wl_pos);
    take_vector_away_from_vector(low_wl_pos, &diff);
    
    double currLength = length_of_vector(diff);
    currLength *= proportion;
    
    scale_vector_to_distance(&diff, currLength);
    add_vector_to_vector(&diff, low_wl_pos);
    
    return diff;
}

double Miller::slicedIntegralWithVectors(vec low_wl_pos, vec high_wl_pos, double rlpSize, double mean, double sigma, double exponent)
{
    double p = 0;
    double q = 1;
    
    double mySlices = slices;
    
    if (resolution() < 1 / trickyRes)
    {
        mySlices = maxSlices;
    }
    
    double interval = 1 / mySlices;
    
    double total_integral = 0;
    double total_normal = 0;
    double total_sphere = 0;
    
    vec oldIncrement = low_wl_pos;
    
    for (double slice = p + interval; slice < q; slice += interval)
    {
        vec newIncrement = incrementedPosition(low_wl_pos, high_wl_pos, slice);
        
        double pBandwidth = getEwaldSphereNoMatrix(oldIncrement);
        double qBandwidth = getEwaldSphereNoMatrix(newIncrement);
        
        double normalSlice = integrate_beam_slice(qBandwidth, pBandwidth, mean, sigma, exponent);
        double sphereSlice = 0;
        
        if (normalSlice > 0)
        {
            sphereSlice = integrate_sphere(slice - interval, slice);
        }
        
        double totalSlice = normalSlice * sphereSlice;
        
        total_integral += totalSlice;
        total_normal += normalSlice;
        total_sphere += sphereSlice;
    }
    
    return total_integral;
}

Miller::Miller(MtzManager *parent, int _h, int _k, int _l, bool calcFree)
{
    h = _h;
    k = _k;
    l = _l;
    lastX = 0;
    lastY = 0;
    rawIntensity = 0;
    sigma = 1;
    wavelength = 0;
    partiality = -1;
    countingSigma = 0;
    polarisationCorrection = 0;
    rejectedReasons = 0;
    scale = 1;
    bFactor = 0;
    bFactorScale = 0;
    resol = 0;
    shift = std::make_pair(0, 0);
    shoebox = ShoeboxPtr();
    flipMatrix = 0;
    lastPeakX = 0;
    lastPeakY = 0;
    
    partialCutoff = FileParser::getKey("PARTIALITY_CUTOFF",
                                       PARTIAL_CUTOFF);
    
    int rlpInt = FileParser::getKey("RLP_MODEL", 0);
    rlpModel = (RlpModel)rlpInt;
  
    if (calcFree)
    {
        free = FreeMillerLibrary::isMillerFree(this);
    }
    
    shiftedRay = new_vector(0, 0, 0);
    mtzParent = parent;
    matrix = MatrixPtr();
}

MatrixPtr Miller::getFlipMatrix()
{
    return Reflection::getFlipMatrix(flipMatrix);
}

vec Miller::hklVector(bool shouldFlip)
{
    vec newVec = new_vector(h, k, l);
    
    if (shouldFlip)
    {
        getFlipMatrix()->multiplyVector(&newVec);
    }
    
    return newVec;
}

int Miller::getH()
{
    return hklVector().h;
}

int Miller::getK()
{
    return hklVector().k;

}

int Miller::getL()
{
    return hklVector().l;
}

void Miller::setFlipMatrix(int i)
{
    flipMatrix = i;
}

void Miller::setParent(ReflectionPtr reflection)
{
    parentReflection = reflection;
}

void Miller::setPartialityModel(PartialityModel _model)
{
    model = _model;
}

void Miller::setData(double _intensity, double _sigma, double _partiality,
                     double _wavelength)
{
    rawIntensity = _intensity;
    sigma = _sigma;
    partiality = _partiality;
    wavelength = _wavelength;
}

bool Miller::accepted(void)
{
    if (this->model == PartialityModelNone)
        return true;

    if (this->partiality < partialCutoff)
    {
        return false;
    }
    
    if (this->isRejected())
        return false;
    
    if (this->rawIntensity != this->rawIntensity)
        return false;
    
    double sigma = this->getSigma();
    
    if (sigma == 0)
    {
        return false;
    }
    
    return true;
}

bool Miller::isRejected()
{
    return (rejectedReasons > 0);
}

double Miller::getBFactorScale()
{
    if (bFactor == 0)
    {
        return 1;
    }
    
    double resn = getResolution();
    
    double four_d_squared = 4 * pow(1 / resn, 2); // = 1 / (4d^2) = (sinTheta / lambda) ^ 2
    
    double factor = exp(-2 * bFactor * 1 / four_d_squared);
    
    bFactorScale = factor;
    
    return factor;
}

double Miller::scaleForScaleAndBFactor(double scaleFactor, double bFactor,
                                      double resn, double exponent_exponent)
{
    double four_d_squared = 4 * pow(resn, 2); // = 1 / (4d^2) = (sinTheta / lambda) ^ 2
    
    double right_exp = exp(-2 * bFactor * four_d_squared);
    
    double scale = scaleFactor * right_exp;
    
    return scale;
}

double Miller::intensity(bool withCutoff)
{
    if ((this->accepted() && withCutoff) || !withCutoff)
    {
        double modifier = scale;
        double bFacScale = getBFactorScale();
        
        modifier /= bFacScale;
        modifier /= correctingPolarisation ? getPolarisationCorrection() : 1;
        
        if (model == PartialityModelScaled)
        {
            modifier /= partiality;
        }
        else if (model == PartialityModelBinary)
        {
            modifier = (partiality < partialCutoff) ? 0 : modifier;
        }
        
        if (partiality == 0)
            modifier = 0;
        
        return modifier * rawIntensity;
    }
    else
    {
        std::cout << "Warning, allowing an unaccepted intensity!" << std::endl;
        return 0;
    }
}

void Miller::makeScalesPermanent()
{
    double scale = scaleForScaleAndBFactor(this->scale, this->bFactor,
                                          1 / this->resol);
    
    this->rawIntensity *= scale;
    this->sigma *= scale;
    
    this->scale = 1;
    bFactor = 0;
}


double Miller::getSigma(void)
{
    // bigger sigma - larger error
    
    double correction = scale;
    /*
    double bFactorScale = getBFactorScale();
    correction /= bFactorScale;
    
    if (bFactorScale == 0) correction = 0;*/
    
    return sigma * correction;
}

double Miller::getPartiality()
{
    return partiality;
}

double Miller::getWavelength()
{
    return wavelength;
}

vec Miller::getTransformedHKL(MatrixPtr matrix)
{
    vec hkl = new_vector(h, k, l);
    
    matrix->multiplyVector(&hkl);
    
    return hkl;
}

vec Miller::getTransformedHKL(double hRot, double kRot)
{
    MatrixPtr newMatrix = MatrixPtr();
    rotateMatrixHKL(hRot, kRot, 0, matrix, &newMatrix);
    
    vec hkl = getTransformedHKL(newMatrix);
    return hkl;
}

double Miller::recalculateWavelength()
{
    double hRot = 0;
    double kRot = 0;
    
    if (getMtzParent() != NULL)
    {
        hRot = getMtzParent()->getHRot();
        kRot = getMtzParent()->getKRot();
    }
    
    getWavelength(hRot, kRot);
    
    return getWavelength();
}

double Miller::getWavelength(double hRot, double kRot)
{
    vec hkl = getTransformedHKL(hRot, kRot);
    
    wavelength = getEwaldSphereNoMatrix(hkl);
    
    return wavelength;
}

double Miller::getWavelength(MatrixPtr transformedMatrix)
{
    vec hkl = getTransformedHKL(transformedMatrix);
    wavelength = getEwaldSphereNoMatrix(hkl);
    
    return wavelength;
}

double Miller::getEwaldWeight(double hRot, double kRot, bool isH)
{
    double weight = 0;
    vec hkl = getTransformedHKL(hRot, kRot);
    weight = getEwaldWeightForAxis(hkl, isH);
    
    return weight;
}

double Miller::getWeight(bool cutoff, WeightType weighting)
{
    if (!this->accepted() && cutoff)
    {
        return 0;
    }
    
    double weight = 1;
    double partiality = this->getPartiality();
    
    if (weighting == WeightTypePartiality
        || weighting == WeightTypePartialityCorrelation
        || weighting == WeightTypePartialitySigma)
    {
        weight = partiality;
    }
    
    if (weighting == WeightTypePartialitySigma)
    {
        double sigma = this->getSigma();
        
        if (sigma > 0)
            weight /= sigma;
    }
    
    if (weighting == WeightTypeISigI)
    {
        double intensity = this->getRawIntensity();
        double sigma = this->getSigma();
        
        weight = intensity / sigma;
    }
    
    weight *= getBFactorScale();

    return weight;
}

void Miller::rotateMatrixHKL(double hRot, double kRot, double lRot, MatrixPtr oldMatrix, MatrixPtr *newMatrix)
{
    (*newMatrix) = oldMatrix->copy();
    
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    double lRad = lRot * M_PI / 180;
    
    (*newMatrix)->rotate(hRad, kRad, lRad);
}

double Miller::expectedRadius(double spotSize, double mosaicity, vec *hkl)
{
    vec usedHKL;
    
    if (hkl == NULL)
    {
        MatrixPtr newMatrix = MatrixPtr();
        rotateMatrixHKL(0, 0, 0, matrix, &newMatrix);
        
        usedHKL = new_vector(h, k, l);
        hkl = &usedHKL;
        
        newMatrix->multiplyVector(hkl);
    }
    
    
    spotSize = fabs(spotSize);
    double radMos = fabs(mosaicity) * M_PI / 180;
    
    double distanceFromOrigin = length_of_vector(*hkl);
    
    double spotSizeIncrease = fabs(radMos * distanceFromOrigin);
    
    double radius = (spotSize + spotSizeIncrease);
    
    return radius;
}

void Miller::limitingEwaldWavelengths(vec hkl, double mosaicity, double spotSize, double wavelength, double *limitLow, double *limitHigh, vec *inwards, vec *outwards)
{
    double radius = expectedRadius(spotSize, mosaicity, &hkl);
    
    vec mean_wavelength_ewald_centre = new_vector(0, 0, -1 / wavelength);
    
    vec centred_spot_position = vector_between_vectors(
                                                       mean_wavelength_ewald_centre, hkl);
    
    double length = length_of_vector(centred_spot_position);
    
    double radiusOverLength = radius / length;
    double inwardsScalar = 1 - radiusOverLength;
    double outwardsScalar = 1 + radiusOverLength;
    double newL = hkl.l + 1 / wavelength;
    
    vec move_inwards_position = new_vector(inwardsScalar * hkl.h, inwardsScalar * hkl.k, hkl.l - radiusOverLength * newL);
    vec move_outwards_position = new_vector(outwardsScalar * hkl.h, outwardsScalar * hkl.k, hkl.l + radiusOverLength * newL);
    
    if (inwards != NULL)
    {
        *inwards = move_inwards_position;
    }
    
    if (outwards != NULL)
    {
        *outwards = move_outwards_position;
    }
    
    double inwards_bandwidth = getEwaldSphereNoMatrix(move_inwards_position);
    double outwards_bandwidth = getEwaldSphereNoMatrix(move_outwards_position);
    
    *limitHigh = inwards_bandwidth;
    *limitLow = outwards_bandwidth;
}

double Miller::partialityForHKL(vec hkl, double mosaicity,
                               double spotSize, double wavelength, double bandwidth, double exponent, bool binary)
{
    double radius = expectedRadius(spotSize, mosaicity, &hkl);
    
    bandwidth = fabs(bandwidth);
    
    double inwards_bandwidth = 0; double outwards_bandwidth = 0;
    vec inwards_vec, outwards_vec;
    
    limitingEwaldWavelengths(hkl, mosaicity, spotSize, wavelength, &outwards_bandwidth, &inwards_bandwidth, &inwards_vec, &outwards_vec);
    
    double stdev = wavelength * bandwidth / 2;
    
    double limit = 4 * stdev;
    double inwardsLimit = wavelength + limit;
    double outwardsLimit = wavelength - limit;
    
    if ((outwards_bandwidth > inwardsLimit || inwards_bandwidth < outwardsLimit) && resolution() > (2 / trickyRes))
    {
        return 0;
    }
    
    double thisPartiality = 0;
    
    if (correctedPartiality)
    {
        thisPartiality = slicedIntegralWithVectors(inwards_vec, outwards_vec, radius, wavelength, stdev, exponent);
    }
    else
    {
        thisPartiality = sliced_integral(outwards_bandwidth, inwards_bandwidth, radius, 0, 1, wavelength, stdev, exponent);
    }
    
    return thisPartiality;
    
}

double Miller::calculateNormPartiality(MatrixPtr rotatedMatrix, double mosaicity,
                                      double spotSize, double wavelength, double bandwidth, double exponent)
{
    vec hkl = new_vector(h, k, l);
    
    rotatedMatrix->multiplyVector(&hkl);
    
    double d = length_of_vector(hkl);
    
    double newH = 0;
    double newK = sqrt((4 * pow(d, 2) - pow(d, 4) * pow(wavelength, 2)) / 4);
    double newL = 0 - pow(d, 2) * wavelength / 2;
    
    vec newHKL = new_vector(newH, newK, newL);
    
    double normPartiality = partialityForHKL(newHKL, mosaicity,
                                            spotSize, wavelength, bandwidth, exponent, false);
    
    return normPartiality;
}

double Miller::resolution()
{
    if (resol == 0)
    {
        vec newVec = getTransformedHKL(0, 0);
        
        resol = length_of_vector(newVec);
    }
    
    return resol;
}

bool Miller::crossesBeamRoughly(MatrixPtr rotatedMatrix, double mosaicity,
                         double spotSize, double wavelength, double bandwidth)
{
    vec hkl = new_vector(h, k, l);
    
    rotatedMatrix->multiplyVector(&hkl);
    
    double inwards_bandwidth = 0; double outwards_bandwidth = 0;
    
    limitingEwaldWavelengths(hkl, mosaicity, spotSize, wavelength, &outwards_bandwidth, &inwards_bandwidth);
    
    double limit = 2 * bandwidth;
    double inwardsLimit = wavelength + limit;
    double outwardsLimit = wavelength - limit;
    
    if (outwards_bandwidth > inwardsLimit || inwards_bandwidth < outwardsLimit)
    {
        partiality = 0;
        return false;
    }
    
    partiality = 1;
    return true;
}


void Miller::recalculatePartiality(MatrixPtr rotatedMatrix, double mosaicity,
                                   double spotSize, double wavelength, double bandwidth, double exponent, bool binary)
{
    if (model == PartialityModelFixed)
    {
        return;
    }
    
    vec hkl = new_vector(h, k, l);
    
    rotatedMatrix->multiplyVector(&hkl);
    
    this->wavelength = getEwaldSphereNoMatrix(hkl);
    
    double tempPartiality = partialityForHKL(hkl, mosaicity,
                                            spotSize, wavelength, bandwidth, exponent, true);
    
    if (binary && tempPartiality == 1.0)
        return;
    
    double normPartiality = 1;
    
    if (normalised && tempPartiality > 0)
    {
        normPartiality = calculateNormPartiality(rotatedMatrix, mosaicity, spotSize, wavelength, bandwidth, exponent);
    }
    else
    {
        partiality = 0;
        return;
    }
    
    partiality = tempPartiality / normPartiality;
    
    if (partiality > 2.0)
        partiality = 0;
    if (partiality > 2.0)
        partiality = 2;
    
    if ((!std::isfinite(partiality)) || (partiality != partiality))
    {
        partiality = 0;
    }
}

double Miller::twoTheta(bool horizontal)
{
    double usedWavelength = wavelength;
    
    double rlpCoordinates[2];
    double beamCoordinates[2];
    
    vec hkl = new_vector(h, k, l);
    matrix->multiplyVector(&hkl);
    
    
    vec crystalCentre = new_vector(0, 0, -1 / usedWavelength);
    vec crystalToRlp = vector_between_vectors(crystalCentre, hkl);
    
    rlpCoordinates[0] = (horizontal ? crystalToRlp.h : crystalToRlp.k);
    rlpCoordinates[1] = crystalToRlp.l;
    
    beamCoordinates[0] = 0;
    beamCoordinates[1] = 1 / usedWavelength;
    
  //  double cosTwoTheta = dotProduct / (rlpMagnitude * beamMagnitude);
    
    double hOrK = (horizontal ? hkl.h : hkl.k);
    double effectiveResolution = 1 / sqrt(hOrK * hOrK + hkl.l * hkl.l);
    double sinTheta = usedWavelength / (2 * effectiveResolution);
    double theta = asin(sinTheta);
    
    return 2 * theta;
}

void Miller::setHorizontalPolarisationFactor(double newFactor)
{
    polarisationCorrection = 0;
    
    // not finished
}

double Miller::getPolarisationCorrection()
{
    if (!matrix)
    {
        return 1;
    }
    
    if (polarisationCorrection == 0)
    {
        // P = (1 + cos2 2θM cos2 2θ) / (1 + cos2 2θM)
        
        // calculate polarisation factor

        double horizontalFactor = polarisationFactor;
        double verticalFactor = 1 - polarisationFactor;
        
        double horizontalTwoTheta = twoTheta(true);
        double verticalTwoTheta = twoTheta(false);
        
        double cosSquaredHorizontal = pow(cos(horizontalTwoTheta), 2);
        double cosSquaredVertical = pow(cos(verticalTwoTheta), 2);
        
        double numerator1 = 1 * verticalFactor + cosSquaredHorizontal * horizontalFactor;
        double denominator1 = verticalFactor + horizontalFactor;
        
        double numerator2 = 1 * horizontalFactor + cosSquaredVertical * verticalFactor;
        double denominator2 = horizontalFactor + verticalFactor;
        
        double numerator = (numerator1 + numerator2);
        double denominator = (denominator1 + denominator2);
        
        polarisationCorrection = numerator / denominator;
    }
    
    return polarisationCorrection;
}

MillerPtr Miller::copy(void)
{
    MillerPtr newMiller = MillerPtr(new Miller(mtzParent));
    
    newMiller->h = h;
    newMiller->k = k;
    newMiller->l = l;
    newMiller->rawIntensity = rawIntensity;
    newMiller->model = model;
    newMiller->sigma = sigma;
    newMiller->partiality = partiality;
    newMiller->wavelength = wavelength;
    newMiller->matrix = matrix;
    newMiller->mtzParent = mtzParent;
    newMiller->countingSigma = countingSigma;
    newMiller->rejectedReasons = rejectedReasons;
    newMiller->lastX = lastX;
    newMiller->lastY = lastY;
    newMiller->shift = shift;
    
    return newMiller;
}

void Miller::applyScaleFactor(double scaleFactor)
{
    setScale(scale * scaleFactor);
}

void Miller::applyPolarisation(double wavelength)
{
    double components[] =
    { -0.0046462837103, 0.00749971430969, -0.00347981101231, 0.00531577207698,
        0.00576709258158, 0.00533160033255, 0.00633224815006,
        0.000661573996986, -0.00702906145495 };
    
    Matrix *mat = new Matrix(components);
    
    vec hkl = new_vector(h, k, l);
    mat->multiplyVector(&hkl);
    double resolution = length_of_vector(hkl);
    
    double sintheta = wavelength / (2 / resolution);
    double theta = asin(sintheta);
    
    double l = sin(0.5 * theta);
    
    if (l == l && l != 0 && std::isfinite(l))
        rawIntensity *= l + 1;
}

bool Miller::positiveFriedel(bool *positive, int *_isym)
{
    int h = getH();
    int k = getK();
    int l = getL();
    
    cctbx::miller::index<> newMiller = cctbx::miller::index<>(h, k, l);
    space_group *spaceGroup = getParentReflection()->getSpaceGroup();
    asu *asymmetricUnit = getParentReflection()->getAsymmetricUnit();
    
    cctbx::miller::asym_index asymmetricMiller = cctbx::miller::asym_index(*spaceGroup, *asymmetricUnit, newMiller);
    
    int isym = asymmetricMiller.isym();
    
    //  std::cout << isym << std::endl;
    
    *positive = ((isym > 0) == 1);
    
    return isym != 0;
}

void Miller::positionOnDetector(MatrixPtr transformedMatrix, int *x,
                                int *y, bool shouldSearch)
{
    double x_coord = 0;
    double y_coord = 0;
    
    if (!transformedMatrix)
    {
        if (!getIOMRefiner())
            throw 1;
        
        transformedMatrix = getIOMRefiner()->getMatrix();
    }
    
    vec hkl = new_vector(h, k, l);
    transformedMatrix->multiplyVector(&hkl);
    double tmp = hkl.k;
    hkl.k = hkl.h;
    hkl.h = -tmp;
    
    std::pair<double, double> coord = getImage()->reciprocalCoordinatesToPixels(hkl);
    
    lastX = int(coord.first);
    lastY = int(coord.second);
    
    *x = -INT_MAX;
    *y = -INT_MAX;
    
    PanelPtr panel = Panel::panelForMiller(this);
    
    x_coord = coord.first;
    y_coord = coord.second;
    
    
    bool even = shoebox->isEven();

    int intLastX = int(x_coord);
    int intLastY = int(y_coord);
    
    if (even)
    {
        intLastX = round(x_coord);
        intLastY = round(y_coord);
    }
    
    lastX = x_coord;
    lastY = y_coord;
    
    if (Detector::isActive())
    {
        double imageWavelength = getImage()->getWavelength();
        hkl.k *= -1;
        hkl.l += 1 / imageWavelength;
        double xSpot, ySpot;
        
        DetectorPtr detector = Detector::getMaster()->spotCoordForRayIntersection(hkl, &xSpot, &ySpot);
        
        if (!detector)
        {
            return;
        }
        
        lastX = xSpot + 0.5;
        lastY = ySpot + 0.5;
        
        int xInt = xSpot;
        int yInt = ySpot;
        
        if (shouldSearch || (lastPeakX == 0 && lastPeakY == 0))
        {
            int search = getIOMRefiner()->getSearchSize();
            getImage()->focusOnAverageMax(&xInt, &yInt, search, peakSize, even);
        }
        else
        {
            xInt = lastPeakX;
            yInt = lastPeakY;
        }
        
        shift = std::make_pair(xInt + 0.5 - lastX, yInt + 0.5 - lastY);
        
        detector->spotCoordToAbsoluteVec(xInt + 0.5, yInt + 0.5, &shiftedRay);
        
        lastPeakX = xInt;
        lastPeakY = yInt;
        
        *x = xInt;
        *y = yInt;
    }
    else if (!Panel::shouldUsePanelInfo())
    {
        if (!panel)
            return;
        
        int search = getIOMRefiner()->getSearchSize();
        getImage()->focusOnAverageMax(&intLastX, &intLastY, search, peakSize, even);
        
        shift = std::make_pair(intLastX + 0.5 - x_coord, intLastY + 0.5 - y_coord);
        
        *x = intLastX;
        *y = intLastY;
        
        lastPeakX = intLastX;
        lastPeakY = intLastY;
    }
    else
    {
        int search = getIOMRefiner()->getSearchSize();
        Coord bestShift = panel->shiftForMiller(this);
        
        if (bestShift.first != FLT_MAX)
        {
            double shiftedX = lastX + bestShift.first;
            double shiftedY = lastY + bestShift.second;
            
            int xInt, yInt;
            
            if (shouldSearch || (lastPeakX == 0 && lastPeakY == 0))
            {
                xInt = shiftedX;
                yInt = shiftedY;
            
                getImage()->focusOnAverageMax(&xInt, &yInt, search, peakSize, even);
            }
            else
            {
                xInt = lastPeakX;
                yInt = lastPeakY;
            }
            
            shift = std::make_pair(xInt + 0.5 - shiftedX, yInt + 0.5 - shiftedY);
            
            *x = xInt;
            *y = yInt;
            
            lastPeakX = xInt;
            lastPeakY = yInt;
        }
        else
        {
            int search = getIOMRefiner()->getSearchSize();
            
            getImage()->focusOnAverageMax(&intLastX, &intLastY, search, peakSize, even);
            
            lastPeakX = intLastX;
            lastPeakY = intLastY;
        }
    }
    
    correctedX = *x;
    correctedY = *y;
}

void Miller::makeComplexShoebox(double wavelength, double bandwidth, double mosaicity, double rlpSize)
{
    double radius = expectedRadius(rlpSize, mosaicity, NULL);
    
    shoebox->complexShoebox(wavelength, bandwidth, radius);
}

void Miller::integrateIntensity(MatrixPtr transformedMatrix)
{
    if (!getImage())
        throw 1;
    
    std::ostringstream logged;
    
    if (!shoebox)
    {
        shoebox = ShoeboxPtr(new Shoebox(shared_from_this()));
        
        int foregroundLength = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",
                                                  SHOEBOX_FOREGROUND_PADDING);
        int neitherLength = FileParser::getKey("SHOEBOX_NEITHER_PADDING",
                                               SHOEBOX_NEITHER_PADDING);
        int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",
                                                  SHOEBOX_BACKGROUND_PADDING);
        bool shoeboxEven = FileParser::getKey("SHOEBOX_MAKE_EVEN", false);
        
        logged << "Shoebox created from values " << foregroundLength << ", " << neitherLength << ", " << backgroundLength << std::endl;
        
        shoebox->simpleShoebox(foregroundLength, neitherLength, backgroundLength, shoeboxEven);
    }
    
    int x = 0;
    int y = 0;
    
    positionOnDetector(transformedMatrix, &x, &y);
    
    if (x == -INT_MAX || y == -INT_MAX)
    {
        rawIntensity = nan(" ");
        return;
    }
    
    rawIntensity = getImage()->intensityAt(x, y, shoebox, &countingSigma, 0);
}

void Miller::incrementOverlapMask(double hRot, double kRot)
{
    int x = lastX;
    int y = lastY;
    
    if (!shoebox)
        return;
    
    getImage()->incrementOverlapMask(x, y, shoebox);
}


bool Miller::isOverlappedWithSpots(std::vector<SpotPtr> *spots, bool actuallyDelete)
{
    Coord fullShift = Panel::shiftForMiller(this);
    
    double x = fullShift.first + lastX + shift.first;
    double y = fullShift.second + lastY + shift.second;
    int count = 0;
    double tolerance = 2.5;
    
    for (int i = 0; i < spots->size(); i++)
    {
        double x2 = (*spots)[i]->getRawXY().first;
        double y2 = (*spots)[i]->getRawXY().second;
        
        double xDiff = fabs(x2 - x);
        double yDiff = fabs(y2 - y);

        if (xDiff < tolerance && yDiff < tolerance)
        {
            if (actuallyDelete)
            {
                spots->erase(spots->begin() + i);
                i--;
            }
            count++;
        }
    }
    
    return (count > 0);
}

bool Miller::isOverlapped()
{
    int x = lastX;
    int y = lastY;
    
    unsigned char max = getImage()->maximumOverlapMask(x, y, shoebox);
    
    return (max >= 2);
}

void Miller::setRejected(RejectReason reason, bool rejection)
{
    if (rejection)
    {
        rejectedReasons = rejectedReasons | reason;
    }
    else
    {
        // removing a flag. Flip the bits so 00001000 becomes 11110111.
        // Then & the current information with the new information
        
        unsigned int flipped = ~reason;
        rejectedReasons = rejectedReasons & flipped;
    }
}

bool Miller::isRejected(RejectReason reason)
{   
    return (rejectedReasons > 0);
}

double Miller::averageRawIntensity(vector<MillerPtr> millers)
{
    double allIntensities = 0;
    int num = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        
        double anIntensity = miller->getRawIntensity();
        
        if (anIntensity != anIntensity)
            continue;
        
        allIntensities += anIntensity;
        num++;
    }
    
    allIntensities /= num;
    
    return allIntensities;
}

double Miller::observedPartiality(double reference)
{
    return rawIntensity / reference;
}

double Miller::observedPartiality(MtzManager *reference)
{
    return getParentReflection()->observedPartiality(reference, this);
}

Miller::~Miller()
{

}

void Miller::recalculateBetterPartiality()
{
    MatrixPtr bestMatrix = mtzParent->getLatestMatrix();
    double rlpSize = mtzParent->getSpotSize();
    double mosaicity = mtzParent->getMosaicity();
    double limitLow = 0;
    double limitHigh = 0;
    double wavelength = beam->getNominalWavelength();
    
    vec hkl = getTransformedHKL(bestMatrix);
    
    limitingEwaldWavelengths(hkl, mosaicity, rlpSize, wavelength, &limitLow, &limitHigh);
    
    if (beam->nonZeroPartialityExpected(limitLow, limitHigh))
    {
        partiality = 0;
        return;
    }
    
    partiality = sliced_integral(limitLow, limitHigh, rlpSize, 0, 1, 0, 0, 0, false, true);
    
    if (partiality != partiality)
        partiality = 0;
    
    if (partiality == 0)
        return;
    
    // normalise
    
    double d = length_of_vector(hkl);
    
    double newH = 0;
    double newK = sqrt((4 * pow(d, 2) - pow(d, 4) * pow(wavelength, 2)) / 4);
    double newL = 0 - pow(d, 2) * wavelength / 2;
    
    vec newHKL = new_vector(newH, newK, newL);
    
    limitingEwaldWavelengths(newHKL, mosaicity, rlpSize, wavelength, &limitLow, &limitHigh);
    double normPartiality = sliced_integral(limitLow, limitHigh, rlpSize, 0, 1, 0, 0, 0, false, true);
    
    partiality /= normPartiality;
    
    if ((!std::isfinite(partiality)) || (partiality != partiality))
    {
        partiality = 0;
    }
}

RejectReason Miller::getRejectedReason()
{
    if (RejectReasonMerge & rejectedReasons)
    {
        return RejectReasonMerge;
    }
    else if (RejectReasonCorrelation & rejectedReasons)
    {
        return RejectReasonCorrelation;
    }
    else if (RejectReasonPartiality & rejectedReasons)
    {
        return RejectReasonPartiality;
    }
    else
    {
        return RejectReasonNone;
    }
}

double Miller::getRawestIntensity()
{
    return rawIntensity;
}

bool Miller::reachesThreshold()
{
    double iSigI = getRawIntensity() / getCountingSigma();
    
    if (absoluteIntensity)
    {
        return (getRawIntensity() > intensityThreshold);
    }
    
    return (iSigI > intensityThreshold);
}
