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
#include "FileParser.h"
#include "polyfit.hpp"

void RefinementGridSearch::recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results)
{
    size_t paramCount = objects.size();
    size_t workingCount = workingList.size();
    
    if (workingCount < paramCount)
    {
        if (workingCount == 1)
        {
       //     std::cout << "." << std::flush;
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
    
    orderedParams.push_back(workingList);
    orderedResults.push_back(result);
    
    reportProgress(result);
}

void RefinementGridSearch::refine()
{
    RefinementStrategy::refine();
    
    ParamList currentValues;
    CSVPtr csv = CSVPtr(new CSV());

    for (int i = 0; i < objects.size(); i++)
    {
        Getter getter = getters[i];
        currentValues.push_back((*getter)(objects[i]));
        csv->addHeader(tags[i]);
    }
    
    csv->addHeader("result");
    
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
    
    logged << "Setting params ";
    
    
    for (int i = 0; i < minParams.size(); i++)
    {
        Setter setter = setters[i];
        (*setter)(objects[i], minParams[i]);
        
        logged << tags[i] << " = " << minParams[i] << ", ";
    }
    
    double val = (*evaluationFunction)(evaluateObject);
    logged << "score = " << val << std::endl;
    
    logged << std::endl;
    
    csv->writeToFile(jobName + "_gridsearch.csv");
    
    finish();
}

void RefinementGridSearch::assignInterpanelMinimum()
{
    // assuming param0 = pokeY, param1 = pokeX
    
    double angle = FileParser::getKey("NUDGE_ROTATION", 0.0002);
    
    assert(stepSizes.size() == 2);
    
    double pokeXStep = stepSizes[1];
    const double pixRange = angle * 10;
    int gridJumps = pixRange / pokeXStep + 0.5;
    std::vector<std::pair<double, int> > lineDifferences;
    std::vector<std::pair<double, int> > lineSeparations;
    double coverPixels = angle * 5;
    double coverNum = coverPixels / pokeXStep;
    int coverPadding = (coverNum - 1) / 2;
    // y first, then x! Think about it!
    
    double maxSteepness = 0;
    int maxSteepnessValue = 0;
    
    for (int x = 0; x < gridLength; x++)
    {
        int lineMin = 0;
        double lineMinValue = FLT_MAX;
        
        for (int y = 0; y + gridJumps < gridLength; y++)
        {
            int gridPos = x * gridLength + y;
            
            for (int i = 0; i < gridJumps; i++)
            {
                double value = orderedResults[gridPos + i];
                if (value < lineMinValue)
                {
                    lineMinValue = value;
                    lineMin = y + i;
                }
            }
            
            std::vector<double> localRegion, increments;
            
            if (lineMin - coverPadding < 0)
                continue;
            
            if (lineMin + coverPadding >= gridLength)
                continue;
            
            for (int i = lineMin - coverPadding; i <= lineMin + coverPadding; i++)
            {
                int newGridPos = x * gridLength + i;
                localRegion.push_back(orderedResults[newGridPos]);
                increments.push_back(newGridPos);
            }
            
            std::vector<double> polynomial = polyfit(increments, localRegion, 2);
            
            double a = polynomial[2];
            
            if (a > maxSteepness)
            {
                maxSteepness = a;
                maxSteepnessValue = x * gridLength + lineMin;
            }
        }
    }
    
    ParamList minParams = orderedParams[maxSteepnessValue];
    
    for (int i = 0; i < minParams.size(); i++)
    {
        Setter setter = setters[i];
        (*setter)(objects[i], minParams[i]);
    }
    
    double val = (*evaluationFunction)(evaluateObject);
}