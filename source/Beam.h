//
//  Beam.h
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


#ifndef __cppxfel__Beam__
#define __cppxfel__Beam__

#include <stdio.h>
#include "parameters.h"

class Beam
{
private:
    
public:
    Beam();
    
    virtual double integralBetweenEwaldWavelengths(double lowWavelength, double highWavelength) { return 0;};
    virtual double getNominalWavelength() { return 0;};
    
    virtual void addParameters(RefinementStrategyPtr map) {};
    virtual bool nonZeroPartialityExpected(double lowWavelength, double highWavelength) { return false;};
};

#endif /* defined(__cppxfel__Beam__) */
