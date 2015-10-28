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
    double x1, y1, x2, y2;
    projectedXYDisplacement(&x1, &y1);
    spotVector2->projectedXYDisplacement(&x2, &y2);
    
    if (fabs(x1 - x2) > maxDistance || fabs(y1 - y2) > maxDistance)
        return false;
    
    return true;
}

double SpotVector::similarityToSpotVector(SpotVectorPtr spotVector2)
{
    double x1, y1, x2, y2;
    projectedXYDisplacement(&x1, &y1);
    spotVector2->projectedXYDisplacement(&x2, &y2);
    
    vec displace1 = new_vector(x1, y1, 0);
    vec displace2 = new_vector(x2, y2, 0);
    take_vector_away_from_vector(displace1, &displace2);
    
    return length_of_vector(displace2);
}