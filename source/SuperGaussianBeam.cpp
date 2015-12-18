//
//  SuperGaussianBeam.cpp
//  cppxfel
//
//  Created by Helen Ginn on 18/12/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "SuperGaussianBeam.h"
#include <math.h>
#include "FileParser.h"

SuperGaussianBeam::SuperGaussianBeam(double wavelength, double stdev, double exponent)
{
    this->setWavelength(wavelength);
    this->setExponent(exponent);
    this->setStdev(stdev);
    
    wavelengthStep = FileParser::getKey("STEP_SIZE_WAVELENGTH", 0.01);
    exponentStep = FileParser::getKey("STEP_SIZE_EXPONENT", 0.1);
    stdevStep = FileParser::getKey("STEP_SIZE_BANDWIDTH", 0.0001);

    optimisingExponent = FileParser::getKey("OPTIMISING_EXPONENT", false);
    optimisingWavelength = FileParser::getKey("OPTIMISING_WAVELENGTH", false);
    optimisingStdev = FileParser::getKey("OPTIMISING_BANDWIDTH", false);

    if (!optimisingExponent) exponentStep = 0;
    if (!optimisingWavelength) wavelengthStep = 0;
    if (!optimisingStdev) stdevStep = 0;
    
}

void SuperGaussianBeam::makeSuperGaussianLookupTable()
{
    if (!exponentChanged || optimisingExponent)
        return;
    
    superGaussianScale = pow((M_PI / 2), (2 / exponent - 1));
    
    superGaussianTable.clear();
    
    const double min = 0;
    const double max = MAX_SUPER_GAUSSIAN;
    const double step = SUPER_GAUSSIAN_STEP;
    int count = 0;
    
    for (double x = min; x < max; x += step)
    {
        double value = super_gaussian(x, 0, superGaussianScale, exponent);
        
        //      std::cout << x << "\t" << value << std::endl;
        
        superGaussianTable.push_back(value);
        
        count++;
    }
    
    exponentChanged = false;
}

double SuperGaussianBeam::brightnessAtWavelength(double queryWavelength)
{
    if (superGaussianTable.size())
    {
        double standardisedX = fabs((queryWavelength - wavelength) / stdev);
        
        if (standardisedX > MAX_SUPER_GAUSSIAN)
            return 0;
        
        const double step = SUPER_GAUSSIAN_STEP;
        
        int lookupInt = standardisedX / step;
        return superGaussianTable[lookupInt];
    }
    
    return super_gaussian(queryWavelength, wavelength, stdev, exponent);
}

std::vector<SetterFunction> SuperGaussianBeam::getParamSetters()
{
    std::vector<SetterFunction> parameterSetters;
    
    parameterSetters.push_back(setWavelength);
    parameterSetters.push_back(setExponent);
    parameterSetters.push_back(setStdev);

    return parameterSetters;
}

std::vector<GetterFunction> SuperGaussianBeam::getParamGetters()
{
    std::vector<GetterFunction> parameterGetters;
    
    parameterGetters.push_back(getWavelength);
    parameterGetters.push_back(getExponent);
    parameterGetters.push_back(getStdev);
    
    return parameterGetters;
}

std::vector<double> SuperGaussianBeam::getParamSteps()
{
    std::vector<double> steps;
    
    steps.push_back(wavelengthStep);
    steps.push_back(exponentStep);
    steps.push_back(stdevStep);
    
    return steps;
}
