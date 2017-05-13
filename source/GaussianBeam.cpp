//
//  GaussianBeam.cpp
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "GaussianBeam.h"
#include "Vector.h"
#include <math.h>
#include "RefinementStepSearch.h"
#include "FileParser.h"
#include <cmath>

bool GaussianBeam::setupSuperGaussian = false;
vector<double> GaussianBeam::superGaussianTable;
std::mutex GaussianBeam::tableMutex;
double GaussianBeam::superGaussianScale = 0;

double GaussianBeam::superGaussianFromTable(double x, double mean, double sigma, double exponent)
{
    if (x != x || mean != mean)
        return 0;
    
    double standardisedX = fabs((x - mean) / sigma);
    
    if (standardisedX > MAX_SUPER_GAUSSIAN)
        return 0;
    
    if (!std::isfinite(standardisedX) || standardisedX != standardisedX)
        return 0;
    
    const double step = SUPER_GAUSSIAN_STEP;
    
    int lookupInt = standardisedX / step;
    return superGaussianTable[lookupInt];
}

void GaussianBeam::makeSuperGaussianLookupTable(double exponent)
{
    if (setupSuperGaussian)
        return;
    
    tableMutex.lock();
    
    if (setupSuperGaussian)
    {
        tableMutex.unlock();
        return;
    }
    
    superGaussianScale = pow((M_PI / 2), (2 / exponent - 1));
    
    superGaussianTable.clear();
    
    const double min = 0;
    const double max = MAX_SUPER_GAUSSIAN;
    const double step = SUPER_GAUSSIAN_STEP;
    int count = 0;
    
    for (double x = min; x < max; x += step)
    {
        double value = super_gaussian(x, 0, superGaussianScale, exponent);
        
        superGaussianTable.push_back(value);
        
        count++;
    }
    
    setupSuperGaussian = true;
    tableMutex.unlock();
}

GaussianBeam::GaussianBeam(double aMeanWavelength, double aBandwidth, double anExponent)
{
    meanWavelengths.push_back(aMeanWavelength);
    bandwidths.push_back(aBandwidth);
    exponents.push_back(anExponent);
    heights.push_back(1);
    
    optimisingWavelength = FileParser::getKey("OPTIMISING_WAVELENGTH",
                                              !OPTIMISED_WAVELENGTH);
    optimisingBandwidth = FileParser::getKey("OPTIMISING_BANDWIDTH",
                                             !OPTIMISED_BANDWIDTH);
    
    optimisingExponent = FileParser::getKey("OPTIMISING_EXPONENT",
                                            !OPTIMISED_EXPONENT);
    
    stepSizeWavelength = FileParser::getKey("STEP_SIZE_WAVELENGTH",
                                            MEAN_STEP);
    stepSizeBandwidth = FileParser::getKey("STEP_SIZE_BANDWIDTH",
                                           BANDWIDTH_STEP);
    
    stepSizeExponent = FileParser::getKey("STEP_SIZE_EXPONENT",
                                          EXPONENT_STEP);
    
    toleranceWavelength = FileParser::getKey("TOLERANCE_WAVELENGTH",
                                             MEAN_TOLERANCE);
    toleranceBandwidth = FileParser::getKey("TOLERANCE_BANDWIDTH",
                                            BANDWIDTH_TOLERANCE);
    
    toleranceExponent = FileParser::getKey("TOLERANCE_EXPONENT",
                                           EXPO_TOLERANCE);
    
    if (!optimisingExponent)
    {
        makeSuperGaussianLookupTable(anExponent);
    }
}

double GaussianBeam::integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength)
{
    double pValue = valueAtWavelength(lowWavelength);
    double qValue = valueAtWavelength(highWavelength);
 
    double width = fabs(highWavelength - lowWavelength);
    
    double area = (pValue + qValue) / 2 * width;
    
    return area;
}

double GaussianBeam::valueAtWavelength(double aWavelength)
{
    double value = 0;
    
    for (int i = 0; i < meanWavelengths.size(); i++)
    {
        double addition = 0;
        
        if (setupSuperGaussian && exponents.size() == 1)
        {
            addition = superGaussianFromTable(aWavelength, meanWavelengths[i], bandwidths[i], exponents[i]);
        }
        else
        {
            addition = super_gaussian(aWavelength, meanWavelengths[i], bandwidths[i], exponents[i]);
        }
        
        value += addition;
    }
    
    return value;
}

double GaussianBeam::getNominalWavelength()
{
    double sum = 0;
    double weights = 0;
    
    for (int i = 0; i < meanWavelengths.size(); i++)
    {
        sum += meanWavelengths[i] * heights[i];
        weights += heights[i];
    }
    
    sum /= weights;
    
    return sum;
}

void GaussianBeam::addParameters(RefinementStrategyPtr map)
{
    if (optimisingWavelength)
    {
        map->addParameter(this, getWavelengthA, setWavelengthA, stepSizeWavelength, toleranceWavelength);
        
        if (meanWavelengths.size() >= 2)
        {
            map->addParameter(this, getWavelengthB, setWavelengthB, stepSizeWavelength, toleranceWavelength);
        }
    }
    
    if (optimisingBandwidth)
    {
        map->addParameter(this, getBandwidthA, setBandwidthA, stepSizeBandwidth, toleranceBandwidth);
        
        if (meanWavelengths.size() >= 2)
        {
            map->addParameter(this, getBandwidthB, setBandwidthB, stepSizeBandwidth, toleranceBandwidth);
        }
    }
    
    if (optimisingExponent)
    {
        map->addParameter(this, getExponentB, setExponentB, stepSizeExponent, toleranceExponent);
        
        if (meanWavelengths.size() >= 2)
        {
            map->addParameter(this, getExponentB, setExponentB, stepSizeExponent, toleranceExponent);
        }
    }
    
    if (meanWavelengths.size() >= 2)
    {
        map->addParameter(this, getHeightA, setHeightA, stepSizeExponent, toleranceExponent);
        map->addParameter(this, getHeightB, setHeightB, stepSizeExponent, toleranceExponent);
    }
}

bool GaussianBeam::nonZeroPartialityExpected(double lowWavelength, double highWavelength)
{
    if (highWavelength < getNominalWavelength() - bandwidths[0] * 3)
        return true;
    if (lowWavelength > getNominalWavelength() + bandwidths[0] * 3)
        return true;
    
    return false;
}