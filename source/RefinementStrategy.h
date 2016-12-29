//
//  RefinementStrategy.h
//  cppxfel
//
//  Created by Helen Ginn on 29/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__RefinementStrategy__
#define __cppxfel__RefinementStrategy__

#include <stdio.h>
#include "parameters.h"
#include "LoggableObject.h"

class RefinementStrategy : public LoggableObject
{
protected:
    Getter evaluationFunction;
    int maxCycles;
    void *evaluateObject;
    LogLevel priority;
    std::string jobName;
    int cycleNum;
    
<<<<<<< Updated upstream
=======
    std::vector<int> couplings;
>>>>>>> Stashed changes
    std::vector<void *> objects;
    std::vector<Getter> getters;
    std::vector<Setter> setters;
    std::vector<double> stepSizes;
    std::vector<double> stepConvergences;
    std::vector<std::string> tags;
    
    void reportProgress(double score);
    void finish();
public:
    RefinementStrategy()
    {
        evaluationFunction = NULL;
        maxCycles = 30;
        priority = LogLevelDebug;
        cycleNum = 0;
    };
    
    virtual void refine();
    
    void addParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag = "");
<<<<<<< Updated upstream
=======
    void addCoupledParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag = "");
>>>>>>> Stashed changes
    
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
    
    virtual void clearParameters()
    {
        getters.clear();
        setters.clear();
        objects.clear();
        stepSizes.clear();
        stepConvergences.clear();
        tags.clear();
    }
};

#endif /* defined(__cppxfel__RefinementStrategy__) */
