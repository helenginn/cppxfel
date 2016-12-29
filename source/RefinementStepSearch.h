//
//  GetterSetterMap.h
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__GetterSetterMap__
#define __cppxfel__GetterSetterMap__

#include <stdio.h>
#include "parameters.h"
#include "RefinementStrategy.h"

class RefinementStepSearch : public RefinementStrategy, public LoggableObject
{
private:
    
public:
    RefinementStepSearch() : RefinementStrategy()
    {
        
    };
    
    void refine();
    
};

#endif /* defined(__cppxfel__GetterSetterMap__) */
