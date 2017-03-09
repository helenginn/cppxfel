//
//  AmbiguityBreaker.h
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

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
    
    MtzPtr getMergedMtz()
    {
        return merged;
    }
};

#endif /* defined(__cppxfel__AmbiguityBreaker__) */
