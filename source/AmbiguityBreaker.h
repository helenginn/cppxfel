//
//  AmbiguityBreaker.h
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


#ifndef __cppxfel__AmbiguityBreaker__
#define __cppxfel__AmbiguityBreaker__

#include "MtzManager.h"
#include "parameters.h"
#include <vector>

class StatisticsManager;
typedef std::map<MtzPtr, double> CrystalCorrelation;
typedef std::map<MtzPtr, CrystalCorrelation> CorrelationMap;

class AmbiguityBreaker : public LoggableObject
{
private:
    std::vector<CorrelationMap> correlationMaps;
    vector<MtzPtr> mtzs;
    int ambiguityCount;
    StatisticsManager *statsManager;
    double gridCorrelation(int imageNumI, int imageNumJ);
    double evaluation();
    double gradientForImage(int imageNum, int axis);
    
    double distance(int vectorNum, int centreNum);
    double toggleValue(int slowCloud, int fastCloud);
    double evaluationCloudCluster();
    double gradientCloudCluster(int centre, int axis);
    
    double dotProduct(int imageNumI, int imageNumJ);

	static void plotDifferenceThread(AmbiguityBreaker *me, int offset);
    void assignPartialities();
    void breakAmbiguity();
	static void calculateCorrelations(AmbiguityBreaker *me, int offset);
    void makeCorrelationGrid();
    void printResults();
    void merge();
    MtzPtr merged;
    
    
public:
    void setMtzs(vector<MtzPtr> newMtzs);
    
    AmbiguityBreaker(vector<MtzPtr> newMtzs);
    void run();
    void overrideAmbiguity(int newAmbiguity);
	void plotDifferences();

    MtzPtr getMergedMtz()
    {
        return merged;
    }
};

#endif /* defined(__cppxfel__AmbiguityBreaker__) */
