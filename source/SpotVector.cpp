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