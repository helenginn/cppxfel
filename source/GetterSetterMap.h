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
#include "LoggableObject.h"

typedef enum
{
    GetterSetterStepSearch,
    GetterSetterNelderMead
} GetterSetterRefinementType;

class GetterSetterMap : public LoggableObject
{
private:
    Getter evaluationFunction;
    int maxCycles;
    void *evaluateObject;
    LogLevel priority;
    std::string jobName;
    
    std::vector<void *> objects;
    std::vector<Getter> getters;
    std::vector<Setter> setters;
    std::vector<double> stepSizes;
    std::vector<double> stepConvergences;
    std::vector<std::string> tags;
    double minimizeParameter(int i);
    
public:
    GetterSetterMap()
    {
        evaluationFunction = NULL;
        maxCycles = 30;
        priority = LogLevelDebug;
    };
    
    void refine(GetterSetterRefinementType type);
    
    void addParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag = "");
    void setEvaluationFunction(Getter function, void *evaluatedObject)
    {
        evaluationFunction = function;
        evaluateObject = evaluatedObject;
    }
    
    void setVerbose(bool verbose)
    {
        if (verbose)
        {
            priority = LogLevelNormal;
        }
    }
    
    void setCycles(int num)
    {
        maxCycles = num;
    }
    
    void setJobName(std::string job)
    {
        jobName = job;
    }
    
    void clearParameters()
    {
        getters.clear();
        setters.clear();
        objects.clear();
        stepSizes.clear();
        stepConvergences.clear();
        tags.clear();
    }
};

#endif /* defined(__cppxfel__GetterSetterMap__) */
