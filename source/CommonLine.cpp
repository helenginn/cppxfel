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

void CommonLine::signature(std::vector<double> *resolutions)
{
    for (int i = 0; i < commonSpots.size(); i++)
    {
        SpotPtr spot = commonSpots[i];
        double resolution = spot->resolution();
        if (angle < 0)
            resolution = -resolution;
        resolutions->push_back(resolution);
    }
}

double CommonLine::similarityToCommonLine(CommonLinePtr otherLine)
{
    std::vector<double> thisSignature, otherSignature;
    this->signature(&thisSignature);
    otherLine->signature(&otherSignature);
    bool thisGreaterThanOther = false;
    
    if (thisSignature.size() != otherSignature.size())
    {
        return 1000;
        
        thisGreaterThanOther = (thisSignature.size() > otherSignature.size());
        
        std::vector<double> bigger = thisGreaterThanOther ? thisSignature : otherSignature;
        std::vector<double> smaller = thisGreaterThanOther ? otherSignature : thisSignature;
        std::vector<double> reducedBigger;
        
        for (int i = 0; i < smaller.size(); i++)
        {
            double smallerOne = smaller[i];
            double minDifference = 1000;
            int minJ = -1;
            
            for (int j = 0; j < bigger.size(); j++)
            {
                double diff = (bigger[j] - smallerOne);
                if (minDifference > diff)
                {
                    minDifference = diff;
                    minJ = j;
                }
            }
            
            if (minJ != -1)
            {
                reducedBigger.push_back(bigger[minJ]);
                bigger.erase(bigger.begin() + minJ);
            }
        }
        
        thisSignature = thisGreaterThanOther ? reducedBigger : smaller;
        otherSignature = thisGreaterThanOther ? smaller : reducedBigger;
    }
    
    bool matchingPositive = true;
    bool matchingNegative = true;
    
    for (int i = 0; i < thisSignature.size(); i++)
    {
        if ((thisSignature[i] < 0 && otherSignature[i] < 0) || (thisSignature[i] > 0 && otherSignature[i] > 0))
            matchingNegative = false;
        
        if ((thisSignature[i] < 0 && otherSignature[i] > 0) || (thisSignature[i] > 0 && otherSignature[i] < 0))
            matchingPositive = false;
    }
    
    if (matchingNegative == false && matchingPositive == false)
    {
        return 1000;
    }
    
    double sumSquares = 0;
    
    for (int i = 0; i < thisSignature.size(); i++)
    {
        double thisResolution = thisSignature[i];
        double otherResolution = otherSignature[i];
        
        double squareDifference = pow(fabs(thisResolution) - fabs(otherResolution), 2);
        sumSquares += squareDifference;
    }
    
    sumSquares /= thisSignature.size();
    
    return sqrt(sumSquares);
}

double CommonLine::angleWithCommonLine(CommonLinePtr otherLine)
{
    double otherAngle = otherLine->angle;
    double thisAngle = this->angle;
    
    std::cout << "Common line for image " << getParentImage() << std::endl;
    
    for (int i = 0; i < commonSpots.size(); i++)
    {
        std::cout << commonSpots[i]->getX() << "\t" << commonSpots[i]->getY() << std::endl;
    }
    
    std::cout << "Common line for image " << otherLine->getParentImage() << std::endl;

    for (int i = 0; i < otherLine->commonSpots.size(); i++)
    {
        std::cout << otherLine->commonSpots[i]->getX() << "\t" << otherLine->commonSpots[i]->getY() << std::endl;
    }

    std::cout << std::endl;
    
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
    
    while (difference < -2 * M_PI)
    {
        difference += 2 * M_PI;
    }
    
    if (difference < -M_PI)
    {
        difference = -2 * M_PI + difference;
    }
    
    if (difference > M_PI / 2)
    {
        difference = M_PI - difference;
    }

 
    if (difference < -M_PI / 2)
    {
        difference = -M_PI - difference;
    }

    if (biggerAngle == otherAngle)
        difference = 0 - difference;
    
    return difference;
}

