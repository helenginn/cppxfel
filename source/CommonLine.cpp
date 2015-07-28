//
//  CommonLine.cpp
//  cppxfel
//
//  Created by Helen Ginn on 27/07/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "CommonLine.h"
#include "Spot.h"
#include <algorithm>

static bool higherResolutionThan(SpotPtr one, SpotPtr two)
{
    return (one->resolution() > two->resolution());
}

CommonLine::CommonLine(std::vector<SpotPtr> newSpots, Image *parent)
{
    commonSpots = newSpots;
    parentImage = parent;
    
    std::sort(commonSpots.begin(), commonSpots.end(), higherResolutionThan);
    
    if (commonSpots.size() > 0)
    {
        angle = commonSpots[0]->angleInPlaneOfDetector();
    }
    else
    {
        angle = 0;
    }
}

static bool isGreaterThan(double one, double two)
{
    return (one > two);
}



void CommonLine::signature(std::vector<double> *resolutions)
{
    for (int i = 0; i < commonSpots.size(); i++)
    {
        SpotPtr spot = commonSpots[i];
        double resolution = spot->resolution();
        resolutions->push_back(resolution);
    }
}

double CommonLine::similarityToCommonLine(CommonLinePtr otherLine)
{
    std::vector<double> thisSignature, otherSignature;
    this->signature(&thisSignature);
    otherLine->signature(&otherSignature);
    
    if (thisSignature.size() != otherSignature.size())
        return 10;
    
    double sumSquares = 0;
    
    for (int i = 0; i < thisSignature.size(); i++)
    {
        double thisResolution = thisSignature[i];
        double otherResolution = otherSignature[i];
        
        double squareDifference = pow(thisResolution - otherResolution, 2);
        sumSquares += squareDifference;
    }
    
    sumSquares /= thisSignature.size();
    
    return sqrt(sumSquares);
}

double CommonLine::angleWithCommonLine(CommonLine otherLine)
{
    double otherAngle = otherLine.angle;
    double thisAngle = this->angle;
    
    double biggerAngle = thisAngle > otherAngle ? thisAngle : otherAngle;
    double smallerAngle = thisAngle > otherAngle ? otherAngle : thisAngle;
    
    double difference = biggerAngle - smallerAngle;
    
    while (difference > 2 * M_PI)
    {
        difference -= 2 * M_PI;
    }
    
    if (difference > M_PI)
    {
        difference = 2 * M_PI - difference;
    }
    
    if (biggerAngle == otherAngle)
        difference = 0 - difference;
    
    return difference;
}