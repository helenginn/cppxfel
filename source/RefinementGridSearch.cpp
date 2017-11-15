//
//  RefinementGridSearch.cpp
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

#include "RefinementGridSearch.h"
#include <float.h>
#include "CSV.h"
#include "polyfit.hpp"
#include <iomanip>

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
    sendLog(LogLevelNormal);

    csv->writeToFile(jobName + "_gridsearch.csv");

    finish();
}

void RefinementGridSearch::assignInterpanelMinimum()
{
        // assuming param0 = pokeY, param1 = pokeX

        std::vector<std::pair<double, int> > lineDifferences;
        std::vector<std::pair<double, int> > lineSeparations;
        int coverPadding = gridJumps / 2;

        // y first, then x! Think about it!

        bool set = false;
        double maxSteepness = 0;
        int maxSteepnessValue = (int)(results.size() - 1) / 2;

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

                        if (coverPadding < 2) coverPadding = 2;

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

                        if (increments.size() <= 1)
                        {
                                continue;
                        }

                        std::vector<double> polynomial = polyfit(increments, localRegion, 2);

                        if (polynomial.size() == 0)
                        {
                                continue;
                        }

                        double a = polynomial[2];

                        if (a > maxSteepness)
                        {
                                maxSteepness = a;
                                maxSteepnessValue = x * gridLength + lineMin;
                                set = true;
                        }
                }
        }

        if (!set)
        {
                return;
        }

        ParamList minParams = orderedParams[maxSteepnessValue];

        logged << "Setting interpanel minimum (" << gridJumps << ", " << coverPadding << ") " << std::setprecision(5);

        for (int i = 0; i < minParams.size(); i++)
        {
                Setter setter = setters[i];
                (*setter)(objects[i], minParams[i]);
                logged << tags[i] << "=" << minParams[i] << ", ";
        }

        double val = (*evaluationFunction)(evaluateObject);

        logged << "score = " << val << std::endl;
        sendLog();

}
