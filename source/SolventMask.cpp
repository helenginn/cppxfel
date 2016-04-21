//
//  SolventMask.cpp
//  cppxfel
//
//  Created by Helen Ginn on 21/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "SolventMask.h"

std::vector<SolventMaskPtr> SolventMask::solventMasks;

SolventMask::SolventMask(double lowRes, double highRes)
{
    lowResCutoff = 1 / lowRes;
    highResCutoff = 1 / highRes;
}

void SolventMask::addSolventMask(double lowRes, double highRes)
{
    SolventMaskPtr newSolventMask = SolventMaskPtr(new SolventMask(lowRes, highRes));
    
    solventMasks.push_back(newSolventMask);
}
