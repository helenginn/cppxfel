//
//  RefinementStepSearch.cpp
//  cppxfel
//
//  Created by Helen Ginn on 12/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "RefinementStepSearch.h"
#include "misc.h"
#include <float.h>
#include "MtzManager.h"

double RefinementStepSearch::minimizeParameter(int whichParam, double *bestScore)
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
    
    if (*bestScore == FLT_MAX)
    {
        param_scores[1] = *bestScore;
    }
    else
    {
        double aScore = (*evaluationFunction)(evaluateObject);
        if (aScore != aScore) aScore = FLT_MAX;
        param_scores[1] = aScore;
    }
    
    param_trials[1] = bestParam;
    
    for (double i = bestParam - step; j < 3; i += step * 2)
    {
        (*setter)(object, i);
        
        double aScore = (*evaluationFunction)(evaluateObject);
        
        if (aScore != aScore) aScore = FLT_MAX;
        
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

void RefinementStepSearch::refine()
{
    RefinementStrategy::refine();
    
    double bestScore = FLT_MAX;

    for (int i = 0; i < maxCycles; i++)
    {
        bool allFinished = true;
        
        for (int j = 0; j < objects.size(); j++)
        {
            allFinished *= minimizeParameter(j, &bestScore);
        }
        
        reportProgress(bestScore);
        
        if (allFinished)
        {
            break;
        }
    }
    
    finish();
}