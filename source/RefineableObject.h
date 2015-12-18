//
//  RefineableObject.h
//  cppxfel
//
//  Created by Helen Ginn on 18/12/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__RefineableObject__
#define __cppxfel__RefineableObject__

#include <stdio.h>

class RefineableObject
{
private:
    
public:
    virtual std::vector<SetterFunction> getParamSetters() = 0;
    virtual std::vector<double> getParamSteps() = 0;
    virtual std::vector<GetterFunction> getParamGetters() = 0;
};

#endif /* defined(__cppxfel__RefineableObject__) */
