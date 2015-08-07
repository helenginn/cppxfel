//
//  CommonLine.h
//  cppxfel
//
//  Created by Helen Ginn on 27/07/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__CommonLine__
#define __cppxfel__CommonLine__
#include "parameters.h"

#include <stdio.h>

class Image;

class CommonLine
{
private:
    std::vector<SpotPtr> commonSpots;
    double angle;
    Image *parentImage;
    
public:
    CommonLine(std::vector<SpotPtr> newSpots, Image *parent);
    double angleWithCommonLine(CommonLinePtr otherLine);
    
    void signature(std::vector<double> *twoThetas);
    double similarityToCommonLine(CommonLinePtr otherLine);
    int makeSuccessful();
    void removeSuccess();
    
    std::vector<SpotPtr> spots()
    {
        return commonSpots;
    }
    
    int spotCount()
    {
        return (int)commonSpots.size();
    }
    
    SpotPtr spot(int i)
    {
        return commonSpots[i];
    }
    
    Image *getParentImage()
    {
        return parentImage;
    }
};

#endif /* defined(__cppxfel__CommonLine__) */
