//
//  SuperGaussianBeam.h
//  cppxfel
//
//  Created by Helen Ginn on 18/12/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__SuperGaussianBeam__
#define __cppxfel__SuperGaussianBeam__

#include <stdio.h>
#include "Beam.h"

class SuperGaussianBeam : public Beam
{
private:
    double wavelength;
    double stdev;
    double exponent;
    
    double wavelengthStep;
    double stdevStep;
    double exponentStep;

    bool optimisingExponent;
    bool optimisingStdev;
    bool optimisingWavelength;

    double superGaussianScale;
    std::vector<double> superGaussianTable;
    void makeSuperGaussianLookupTable();
    
    bool exponentChanged;
    
public:
    SuperGaussianBeam(double wavelength, double stdev, double exponent);
    virtual double brightnessAtWavelength(double queryWavelength);
    
    virtual std::vector<SetterFunction> getParamSetters();
    virtual std::vector<GetterFunction> getParamGetters();
    virtual std::vector<double> getParamSteps();

    double getWavelength()
    {
        return wavelength;
    }
    
    double getExponent()
    {
        return exponent;
    }
    
    double getStdev()
    {
        return stdev;
    }
    
    virtual void setWavelength(double length)
    {
        wavelength = length;
    }
    
    void setStdev(double dev)
    {
        stdev = dev;
    }
    
    void setExponent(double expo)
    {
        exponent = expo;
        exponentChanged = true;
    }
    
    static void setWavelength(void *object, double length)
    {
        static_cast<SuperGaussianBeam *>(object)->setWavelength(length);
    }
    
    static void setStdev(void *object, double dev)
    {
        static_cast<SuperGaussianBeam *>(object)->setStdev(dev);
    }
    
    static void setExponent(void *object, double expo)
    {
        static_cast<SuperGaussianBeam *>(object)->setExponent(expo);
    }
    
    static double getWavelength(void *object)
    {
        return static_cast<SuperGaussianBeam *>(object)->getWavelength();
    }
    
    static double getExponent(void *object)
    {
        return static_cast<SuperGaussianBeam *>(object)->getExponent();
    }

    static double getStdev(void *object)
    {
        return static_cast<SuperGaussianBeam *>(object)->getStdev();
    }

};

#endif /* defined(__cppxfel__SuperGaussianBeam__) */
