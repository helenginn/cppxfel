//
//  GetterSetterMap.cpp
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "GetterSetterMap.h"
#include <float.h>
#include "MtzManager.h"

void GetterSetterMap::addParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence)
{
    objects.push_back(object);
    getters.push_back(getter);
    setters.push_back(setter);
    stepSizes.push_back(stepSize);
    stepConvergences.push_back(stepConvergence);
}

double GetterSetterMap::minimizeParameter(int whichParam)
{
    double param_trials[3];
    double param_scores[3];
    
    double step = stepSizes[whichParam];
    
    if (step < stepConvergences[whichParam])
        return 0;
    
    Getter getter = getters[whichParam];
    Setter setter = setters[whichParam];
    void *object = objects[whichParam];
    
    int j = 0;
    int param_min_num = 1;
    
    double bestParam = (*getter)(object);
    
 //   logged << bestParam << "\t:\t";
    
    for (double i = bestParam - step; j < 3; i += step)
    {
        (*setter)(object, i);
        
        double aScore = (*evaluationFunction)(evaluateObject);
        
        if (aScore != aScore) aScore = FLT_MAX;
        
      //  logged << aScore << " (" << i << ")\t";
        
        param_scores[j] = aScore;
        param_trials[j] = i;
        j++;
    }
    
//    logged << std::endl;
//    sendLog();
    
    double param_min_score = param_scores[1];
    
    for (int i = 0; i < 3; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }
    
    (*setter)(object, param_trials[param_min_num]);
    
  //  std::cout << "Param_min_num trial is " << param_trials[param_min_num] << std::endl;
    
    if (param_min_num == 1)
        stepSizes[whichParam] /= 2;
    
    return param_scores[param_min_num];
}


void GetterSetterMap::refine(GetterSetterRefinementType type)
{
    for (int i = 0; i < maxCycles; i++)
    {
        std::ostringstream logged;

        for (int j = 0; j < objects.size(); j++)
        {
            minimizeParameter(j);
            
       //     logged << (*getters[j])(objects[j]) << "\t";
        }
        
        double score = (*evaluationFunction)(evaluateObject);
        
   //     logged << score << std::endl;
   //     Logger::log(logged);
    }
}