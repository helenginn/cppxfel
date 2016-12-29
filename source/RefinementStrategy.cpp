//
//  RefinementStrategy.cpp
//  cppxfel
//
//  Created by Helen Ginn on 29/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "RefinementStrategy.h"
#include "misc.h"

void RefinementStrategy::addParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag)
{
    objects.push_back(object);
    getters.push_back(getter);
    setters.push_back(setter);
    stepSizes.push_back(stepSize);
    stepConvergences.push_back(stepConvergence);
    
    if (!tag.length())
    {
        tag = "object" + i_to_str((int)objects.size());
    }
    
    tags.push_back(tag);
}
