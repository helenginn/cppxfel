//
//  RefinementGridSearch.cpp
//  cppxfel
//
//  Created by Helen Ginn on 25/01/2017.
//  Copyright (c) 2017 Division of Structural Biology Oxford. All rights reserved.
//

#include "RefinementGridSearch.h"
#include <float.h>
#include "CSV.h"

void RefinementGridSearch::recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results)
{
    size_t paramCount = objects.size();
    size_t workingCount = workingList.size();
    
    if (workingCount < paramCount)
    {
        if (workingCount == 1)
        {
            std::cout << "." << std::flush;
        }
        
        for (int i = -gridLength / 2; i <= (int)(gridLength / 2 + 0.5); i++)
        {
            double mean = referenceList[workingCount];
            double step = stepSizes[workingCount];
            double value = mean + i * step;
            
            ParamList extended = workingList;
            extended.push_back(value);
            recursiveEvaluation(referenceList, extended, results);
        }
        
        return;
    }
    
    for (int i = 0; i < workingList.size(); i++)
    {
        Setter setter = setters[i];
        (*setter)(objects[i], workingList[i]);
    }
    
    double result = (*evaluationFunction)(evaluateObject);
    (*results)[workingList] = result;
}

void RefinementGridSearch::refine()
{
    ParamList currentValues;
    CSVPtr csv = CSVPtr(new CSV());

    for (int i = 0; i < objects.size(); i++)
    {
        Getter getter = getters[i];
        currentValues.push_back((*getter)(objects[i]));
        csv->addHeader(tags[i]);
    }
    
    csv->addHeader("result");
    
    ResultMap results;
    recursiveEvaluation(currentValues, ParamList(), &results);
    
    std::cout << std::endl;
    
    double minResult = FLT_MAX;
    ParamList minParams;
    
    
    for (ResultMap::iterator it = results.begin(); it != results.end(); it++)
    {
        if (it->second < minResult)
        {
            minResult = it->second;
            minParams = it->first;
        }
        
        std::vector<double> result = it->first;
        result.push_back(it->second);
        
        csv->addEntry(result);
    }
    
    for (int i = 0; i < minParams.size(); i++)
    {
        Setter setter = setters[i];
        (*setter)(objects[i], minParams[i]);
    }
    
    csv->writeToFile(jobName + "_gridsearch.csv");
}