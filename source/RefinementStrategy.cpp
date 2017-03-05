//
//  RefinementStrategy.cpp
//  cppxfel
//
//  Created by Helen Ginn on 29/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "RefinementStepSearch.h"
#include "NelderMead.h"
#include "RefinementStrategy.h"
#include "FileParser.h"
#include "misc.h"

RefinementStrategyPtr RefinementStrategy::userChosenStrategy()
{
    int miniMethod = FileParser::getKey("MINIMIZATION_METHOD", 1);
    MinimizationMethod method = (MinimizationMethod)miniMethod;
    RefinementStrategyPtr strategy;
    
    switch (method) {
        case MinimizationMethodStepSearch:
            strategy = boost::static_pointer_cast<RefinementStrategy>(RefinementStepSearchPtr(new RefinementStepSearch()));
            break;
        case MinimizationMethodNelderMead:
            strategy = boost::static_pointer_cast<RefinementStrategy>(NelderMeadPtr(new NelderMead()));
            break;
        default:
            break;
    }
    
    int cycles = FileParser::getKey("NELDER_MEAD_CYCLES", 30);
    strategy->setCycles(cycles);
    
    return strategy;
}

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
    couplings.push_back(1);
}

void RefinementStrategy::addCoupledParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag)
{
    couplings.at(couplings.size() - 1)++;
    addParameter(object, getter, setter, stepSize, stepConvergence, tag);
    couplings.at(couplings.size() - 1)++;
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
    
    reportProgress((*evaluationFunction)(evaluateObject));
}

void RefinementStrategy::reportProgress(double score)
{
    if (!(priority <= LogLevelNormal))
        return;
    
    logged << "Cycle " << cycleNum << "\t";
    
    for (int i = 0; i < objects.size(); i++)
    {
        logged << (*getters[i])(objects[i]) << "\t";
    }

    logged << " - score: ";
    logged << score << std::endl;
    cycleNum++;
}

void RefinementStrategy::finish()
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.clear();
    logged.str("");
    
    cycleNum = 0;
}