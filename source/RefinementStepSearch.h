//
//  RefinementStepSearch.h
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__RefinementStepSearch__
#define __cppxfel__RefinementStepSearch__

#include <stdio.h>
#include "parameters.h"
#include "RefinementStrategy.h"

class RefinementStepSearch : public RefinementStrategy
{
private:
    double minimizeParameter(int i, double *bestScore);

public:
    RefinementStepSearch() : RefinementStrategy()
    {
        
    };
    
    virtual void refine();
    
};

#endif /* defined(__cppxfel__RefinementStepSearch__) */
