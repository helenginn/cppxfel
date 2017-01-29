//
//  RefinementGridSearch.h
//  cppxfel
//
//  Created by Helen Ginn on 25/01/2017.
//  Copyright (c) 2017 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__RefinementGridSearch__
#define __cppxfel__RefinementGridSearch__

#include <stdio.h>
#include "RefinementStrategy.h"

typedef std::vector<double> ParamList;
typedef std::map<ParamList, double> ResultMap;

class RefinementGridSearch : public RefinementStrategy
{
private:
    int gridLength;
    
public:
    RefinementGridSearch() : RefinementStrategy()
    {
        gridLength = 15;
        cycleNum = 1;
    };
    
    void recursiveEvaluation(ParamList referenceList, ParamList workingList, ResultMap *results);
    virtual void refine();
};

#endif /* defined(__cppxfel__RefinementGridSearch__) */