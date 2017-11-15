//
//  RefinementStrategy.cpp
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
                case MinimizationMethodGridSearch:
                        strategy = boost::static_pointer_cast<RefinementStrategy>(RefinementGridSearchPtr(new RefinementGridSearch()));
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

        sendLog(priority);

    if (endScore >= startingScore || endScore != endScore)
    {
        logged << "No change for " << jobName << " (" << startingScore << ")" << std::endl;
        sendLog(LogLevelDetailed);

                resetToInitialParameters();
                _changed = 0;
    }
    else
    {
        double reduction = (startingScore - endScore) / startingScore;

        logged << "Reduction ";

        if (reduction == reduction)
        {
            logged << "by " << std::fixed << std::setprecision(4) <<
            reduction * 100 << "% ";
        }

        logged << "for " << jobName << ": ";

        for (int i = 0; i < objects.size(); i++)
        {
            double objectValue = (*getters[i])(objects[i]);
            logged << tags[i] << "=" << objectValue << ", ";
        }

                logged << "(" << startingScore << " to " << endScore << ")" << std::endl;
        sendLog(LogLevelDetailed);

                _changed = 1;
    }

    cycleNum = 0;

        if (finishFunction != NULL)
        {
                (*finishFunction)(evaluateObject);
        }
}

void RefinementStrategy::resetToInitialParameters()
{
        for (int i = 0; i < objects.size(); i++)
        {
                double objectValue = startingValues[i];
                (*setters[i])(objects[i], objectValue);
        }
}
