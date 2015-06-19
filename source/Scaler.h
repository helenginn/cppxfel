//
//  Scaler.h
//  cppxfel
//
//  Created by Helen Ginn on 01/04/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Scaler__
#define __cppxfel__Scaler__

#include <stdio.h>
#include <vector>
#include "MtzManager.h"

typedef std::map<MtzPtr, std::vector<Holder *> > MtzDataMap;

class Scaler
{
private:
    MtzDataMap mtzData;
    MtzManager **groupedMtz;
    double evaluate();
    double evaluateForImage(MtzPtr mtz);
    double gradientForImageParameter(MtzPtr mtz, int paramNum);
    double hessianGradientForParamType(int paramNum);
    bool isRefiningParameter(int paramNum);
    double stepForParam(int paramNum);
    
    double xStepNorm(double step, int paramNum);
    double gNorm(int paramNum);
    double hessianGradientForParam(MtzPtr mtz, int paramNum);
    void calculateDiagonals();
    void calculateGradients();
    void loadParametersFromMtzs();
    void loadParametersIntoMtzs();
    
    int paramsPerImage();
    int parameterCount();
    scitbx::af::shared<double> g;
    scitbx::af::shared<double> diag;
    scitbx::af::shared<double> x;
public:
    Scaler(std::vector<MtzPtr> mtzs, MtzManager **grouped);
    void minimizeRMerge();
};

#endif /* defined(__cppxfel__Scaler__) */
