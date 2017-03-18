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
#include <iomanip>

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
    
    startingScore = (*evaluationFunction)(evaluateObject);
    
    for (int i = 0; i < objects.size(); i++)
    {
        double objectValue = (*getters[i])(objects[i]);
        startingValues.push_back(objectValue);
    }

    
    reportProgress(startingScore);
}

void RefinementStrategy::reportProgress(double score)
{
    if (!(priority <= LogLevelNormal))
        return;
    
    logged << "Cycle " << cycleNum << "\t";
    
    for (int i = 0; i < objects.size(); i++)
    {
        double objectValue = (*getters[i])(objects[i]);
        logged << std::setprecision(5) << objectValue << "\t";
    }

    logged << " - score: ";
    logged << score << std::endl;
    cycleNum++;
}

void RefinementStrategy::finish()
{
    double endScore = (*evaluationFunction)(evaluateObject);
    
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
    logged.clear();
    logged.str("");
    
    if (endScore > startingScore || endScore != endScore)
    {
        logged << "No change for " << jobName << " (" << startingScore << ")" << std::endl;
        sendLog();
        
        for (int i = 0; i < objects.size(); i++)
        {
            double objectValue = startingValues[i];
            (*setters[i])(objects[i], objectValue);
        }
    }
    else
    {
        double reduction = (startingScore - endScore) / startingScore;
        
        logged << "Reduction by " << std::fixed << std::setprecision(3) <<
        reduction * 100 << "% for " << jobName << ": ";
        
        for (int i = 0; i < objects.size(); i++)
        {
            double objectValue = (*getters[i])(objects[i]);
            logged << tags[i] << "=" << objectValue << ", ";
        }
        logged << std::endl;
        sendLog();
    }
    
    cycleNum = 0;
}