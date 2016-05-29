//
//  SpotFinderCorrelation.cpp
//  cppxfel
//
//  Created by Helen Ginn on 29/05/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotFinderCorrelation.h"
#include "Image.h"
/*
void SpotFinderCorrelation::focusOnMax(int *x, int *y)
{
    double bestX = *x;
    double bestY = *y;
    double maxValue = 0;
    
    for (int j = - jump; j < + jump; j++)
    {
        for (int i = - jump; i < +jump; i++)
        {
            double value = image->valueAt(*x + i, *y + j);
            
            if (value > maxValue)
            {
                maxValue = value;
                bestX = *x + i;
                bestY = *y + j;
            }
        }
    }
    
    *x = bestX;
    *y = bestY;
}
*/
void SpotFinderCorrelation::findSpecificSpots(std::vector<SpotPtr> *spots)
{
    int xDim = this->image->getXDim();
    int yDim = this->image->getYDim();
    
    SpotPtr spot = SpotPtr(new Spot(image));
    
    for (int y = jump; y < yDim - jump; y += jump * 2)
    {
        for (int x = jump; x < xDim - jump; x += jump * 2)
        {
            bool success = spot->focusOnNearbySpot(jump, x, y);
            
            if (success)
            {
                spots->push_back(spot);
                spot = SpotPtr(new Spot(image));
            }
            else
            {
                spot->setUpdate();
            }
        }
    }
}