//
//  GetterSetterMap.cpp
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "misc.h"
#include "GetterSetterMap.h"
#include <float.h>
#include "MtzManager.h"

void GetterSetterMap::addParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag)
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

double GetterSetterMap::minimizeParameter(int whichParam)
{
    double param_trials[3];
    double param_scores[3];
    
    double step = stepSizes[whichParam];
    
    if (step < stepConvergences[whichParam])
        return 1;
    
    Getter getter = getters[whichParam];
    Setter setter = setters[whichParam];
    void *object = objects[whichParam];
    
    int j = 0;
    int param_min_num = 1;
    
    double bestParam = (*getter)(object);
    
    for (double i = bestParam - step; j < 3; i += step)
    {
        (*setter)(object, i);
        
        double aScore = (*evaluationFunction)(evaluateObject);
        
        if (aScore != aScore) aScore = FLT_MAX;
        
     //   logged << "(" << aScore << ") ";
        
        param_scores[j] = aScore;
        param_trials[j] = i;
        j++;
    }
    
    double param_min_score = param_scores[1];
    
    for (int i = 0; i < 3; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }
    
    (*setter)(object, param_trials[param_min_num]);
    
    if (param_min_num == 1)
        stepSizes[whichParam] /= 2;
    
    return 0;
}


void GetterSetterMap::refine(GetterSetterRefinementType type)
{
    if (!jobName.length())
    {
        jobName = "Refinement procedure for " + i_to_str((int)objects.size()) + " objects";
    }
    
    logged << "--- " << jobName << " ---" << std::endl;
    logged << "--- Refining ";
    
    for (int i = 0; i < tags.size() - 1; i++)
    {
        logged << tags[i] << ", ";
    }
    
    logged << tags[tags.size() - 1] << " --- " << std::endl;
    
    for (int i = 0; i < maxCycles; i++)
    {
        logged << "Cycle " << i << "\t";
        
        bool allFinished = true;
        
        for (int j = 0; j < objects.size(); j++)
        {
            allFinished *= minimizeParameter(j);
            
            logged << (*getters[j])(objects[j]) << "\t";
        }
        
        logged << " - score: ";
        
        double score = (*evaluationFunction)(evaluateObject);
        
        logged << score << std::endl;
        
        if (allFinished)
        {
            break;
        }
    }
    
    Logger::mainLogger->addStream(&logged, priority);
    logged.clear();
    logged.str("");
}