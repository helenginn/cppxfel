//
//  SpotVector.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotVector.h"

SpotVector::SpotVector(SpotPtr first, SpotPtr second)
{
    firstSpot = first;
    secondSpot = second;
    
    vec firstVector = first->estimatedVector();
    vec secondVector = second->estimatedVector();
    
    spotDiff = copy_vector(secondVector);
    take_vector_away_from_vector(firstVector, &spotDiff);
}

double SpotVector::distance()
{
    return length_of_vector(spotDiff);
}

double SpotVector::angleWithVector(SpotVectorPtr spotVector2)
{
    vec otherDiff = spotVector2->spotDiff;
    
    return angleBetweenVectors(spotDiff, otherDiff);
}

double SpotVector::angleWithVertical()
{
    vec vertical = new_vector(0, 1, 0);
    
    double angle = angleBetweenVectors(spotDiff, vertical);
    
    if (spotDiff.h < 0)
        angle = -angle;
    
    return angle;
}

void SpotVector::projectedXYDisplacement(double *x, double *y)
{
    double length = distance();
    
    double angle = angleWithVertical();
    
    *x = length * sin(angle);
    *y = length * cos(angle);
}

bool SpotVector::isCloseToSpotVector(SpotVectorPtr spotVector2, double maxDistance)
{
    vec spotDiff2 = spotVector2->getSpotDiff();
    
    if (fabs(spotDiff2.h - spotDiff.h) > maxDistance || fabs(spotDiff2.k - spotDiff.k) > maxDistance
        || fabs(spotDiff2.l - spotDiff.l) > maxDistance)
        return false;
    
    return true;
}

double SpotVector::similarityToSpotVector(SpotVectorPtr spotVector2)
{
    vec displace = copy_vector(spotDiff);
    take_vector_away_from_vector(spotVector2->getSpotDiff(), &displace);
    
    return length_of_vector(displace);
}