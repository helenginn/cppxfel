//
//  GaussianBeam.h
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

#ifndef __cppxfel__GaussianBeam__
#define __cppxfel__GaussianBeam__
#include "Beam.h"
#include <mutex>

#include <stdio.h>

class GaussianBeam : public Beam
{
private:
    double optimisingWavelength;
    double optimisingBandwidth;
    double optimisingExponent;
    double stepSizeWavelength;
    double stepSizeBandwidth;
    double stepSizeExponent;
    double toleranceWavelength;
    double toleranceBandwidth;
    double toleranceExponent;

    std::vector<double> meanWavelengths;
    std::vector<double> bandwidths;
    std::vector<double> exponents;
    std::vector<double> heights;

    static std::mutex tableMutex;
    static bool setupSuperGaussian;
    static double superGaussianScale;
    static vector<double> superGaussianTable;

    void makeSuperGaussianLookupTable(double exponent);
    double superGaussianFromTable(double x, double mean, double sigma, double exponent);

public:
    GaussianBeam(double meanWavelength, double bandwidth, double exponent);
    double integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength);
    void addParameters(RefinementStrategyPtr map);
    double getNominalWavelength();
    double valueAtWavelength(double aWavelength);
    bool nonZeroPartialityExpected(double lowWavelength, double highWavelength);

    static void setWavelength(void *object, double aWavelength, int which)
    {
        static_cast<GaussianBeam *>(object)->meanWavelengths[which] = aWavelength;
    }

    static double getWavelength(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->meanWavelengths[which];
    }

    static void setBandwidth(void *object, double aBandwidth, int which)
    {
        static_cast<GaussianBeam *>(object)->bandwidths[which] = aBandwidth;
    }

    static double getBandwidth(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->bandwidths[which];
    }

    static void setExponent(void *object, double anExponent, int which)
    {
        static_cast<GaussianBeam *>(object)->exponents[which] = anExponent;
    }

    static double getExponent(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->exponents[which];
    }

    static void setHeight(void *object, double aHeight, int which)
    {
        static_cast<GaussianBeam *>(object)->heights[which] = aHeight;
    }

    static double getHeight(void *object, int which)
    {
        return static_cast<GaussianBeam *>(object)->heights[which];
    }

    // individual parameter getters

    static double getWavelengthA(void *object)
    {
        return getWavelength(object, 0);
    }

    static double getWavelengthB(void *object)
    {
        return getWavelength(object, 1);
    }

    static double getBandwidthA(void *object)
    {
        return getBandwidth(object, 0);
    }

    static double getBandwidthB(void *object)
    {
        return getBandwidth(object, 1);
    }

    static double getExponentA(void *object)
    {
        return getExponent(object, 0);
    }

    static double getExponentB(void *object)
    {
        return getExponent(object, 1);
    }

    static double getHeightA(void *object)
    {
        return getHeight(object, 0);
    }

    static double getHeightB(void *object)
    {
        return getHeight(object, 1);
    }

    // individual parameter setters

    static void setWavelengthA(void *object, double aWavelength)
    {
        return setWavelength(object, aWavelength, 0);
    }

    static void setWavelengthB(void *object, double bWavelength)
    {
        return setWavelength(object, bWavelength, 1);
    }

    static void setBandwidthA(void *object, double aBandwidth)
    {
        return setBandwidth(object, aBandwidth, 0);
    }

    static void setBandwidthB(void *object, double aBandwidth)
    {
        return setBandwidth(object, aBandwidth, 1);
    }

    static void setExponentA(void *object, double aExponent)
    {
        return setExponent(object, aExponent, 0);
    }

    static void setExponentB(void *object, double aExponent)
    {
        return setExponent(object, aExponent, 1);
    }

    static void setHeightA(void *object, double aHeight)
    {
        return setHeight(object, aHeight, 0);
    }

    static void setHeightB(void *object, double aHeight)
    {
        return setHeight(object, aHeight, 1);
    }

};

#endif /* defined(__cppxfel__GaussianBeam__) */
