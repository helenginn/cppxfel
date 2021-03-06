//
//  RefinementStrategy.h
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __cppxfel__RefinementStrategy__
#define __cppxfel__RefinementStrategy__

#include <stdio.h>
#include "parameters.h"
#include "LoggableObject.h"

class RefinementStrategy : public LoggableObject
{
protected:
    Getter evaluationFunction;
        Getter finishFunction;
    int maxCycles;
    void *evaluateObject;
    LogLevel priority;
    std::string jobName;
    int cycleNum;
        int _changed;

    std::vector<int> couplings;
    std::vector<void *> objects;
    std::vector<Getter> getters;
    std::vector<Setter> setters;
    std::vector<double> stepSizes;
    std::vector<double> stepConvergences;
    std::vector<std::string> tags;
    std::vector<double> startingValues;
    double startingScore;

    void reportProgress(double score);
    void finish();
public:
    RefinementStrategy()
    {
        evaluationFunction = NULL;
        maxCycles = 30;
        priority = LogLevelDebug;
        cycleNum = 0;
        startingScore = 0;
                _changed = -1;
                finishFunction = NULL;
    };

    static RefinementStrategyPtr userChosenStrategy();

    virtual void refine();
        void resetToInitialParameters();

    void addParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag = "");
    void addCoupledParameter(void *object, Getter getter, Setter setter, double stepSize, double stepConvergence, std::string tag = "");

    void setEvaluationFunction(Getter function, void *evaluatedObject)
    {
        evaluationFunction = function;
        evaluateObject = evaluatedObject;
    }

        void setFinishFunction(Getter finishFunc)
        {
                finishFunction = finishFunc;
        }

    void setVerbose(bool verbose)
    {
        if (verbose)
        {
            priority = LogLevelNormal;
        }
        else
        {
            priority = LogLevelDebug;
        }
    }

        bool didChange()
        {
                return (_changed == 1);
        }

    void setCycles(int num)
    {
        maxCycles = num;
    }

    void setJobName(std::string job)
    {
        jobName = job;
    }

    void *getEvaluationObject()
    {
        return evaluateObject;
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
