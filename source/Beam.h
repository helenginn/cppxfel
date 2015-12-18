//
//  Beam.h
//  cppxfel
//
//  Created by Helen Ginn on 18/12/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Beam__
#define __cppxfel__Beam__

#include <stdio.h>
#include "parameters.h"
#include "Vector.h"
#include "RefineableObject.h"

class Beam : public RefineableObject
{
private:
    
public:
    virtual void setWavelength(double wave) = 0;
    virtual double brightnessAtWavelength(double queryWavelength) = 0;
};

#endif /* defined(__cppxfel__Beam__) */
