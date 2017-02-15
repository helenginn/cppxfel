//
//  SpotFinderCorrelation.cpp
//  cppxfel
//
//  Created by Helen Ginn on 29/05/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotFinderCorrelation.h"
#include "Image.h"


void SpotFinderCorrelation::findSpecificSpots(std::vector<SpotPtr> *spots)
{
    bool verbose = (FileParser::getKey("VERBOSITY_LEVEL", 0) > 1);
    
    int xDim = this->image->getXDim();
    int yDim = this->image->getYDim();
    
    int *mask = new int[xDim * yDim];
    memset(mask, 0, sizeof(int) * xDim * yDim);
    
    SpotPtr spot = SpotPtr(new Spot(image));
    
    double minCorrelation = FileParser::getKey("IMAGE_MIN_CORRELATION", 0.7);
    
    for (int y = 0; y < yDim; y++)
    {
        for (int x = 0; x < xDim; x++)
        {
            if (mask[y * xDim + x])
            {
                continue;
            }
            
            double value = this->image->valueAt(x, y);
            
            if (value < threshold)
                continue;
            
            
            double correlation = spot->focusOnNearbySpot(0, x, y);
            bool success = (correlation > minCorrelation);
            
            if (success)
            {
                spot->recentreInWindow();
                spots->push_back(spot);
                spot->addToMask(mask, xDim, yDim);
            }
            
            if (verbose)
            {
                logged << "Pixel at (" << x << ", " << y << ") intensity " << value << ", correlation " << correlation << (success ? " - success" : " - abandoned") << std::endl;
                sendLog(LogLevelDetailed);
            }

            spot = SpotPtr(new Spot(image));
        }
    }
    
    delete [] mask;
}