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

void RefinementStrategy::refine()
{
    if (!jobName.length())
    {
        jobName = "Refinement procedure for " + i_to_str((int)objects.size()) + " objects";
    }
    
    logged << "--- " << jobName << " ---";
    
    if (tags.size() == 0)
    {
        logged << " No parameters to refine! Exiting." << std::endl;
        sendLog();
        return;
    }
    
    logged << " Refining ";
    
    for (int i = 0; i < tags.size() - 1; i++)
    {
        logged << tags[i] << ", ";
    }
    
    logged << tags[tags.size() - 1] << " --- " << std::endl;
}

void RefinementStrategy::reportProgress(double score)
{
    cycleNum++;
    logged << "Cycle " << cycleNum << "\t";
    
    for (int i = 0; i < objects.size(); i++)
    {
        logged << (*getters[i])(objects[i]) << "\t";
        logged << " - score: ";
        logged << score << std::endl;
    }
}

void RefinementStrategy::finish()
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.clear();
    logged.str("");
}