//
//  SpotFinderCorrelation.cpp
//  cppxfel
//
//  Created by Helen Ginn on 29/05/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotFinderCorrelation.h"
#include "Image.h"
#include <float.h>

void SpotFinderCorrelation::findSpecificSpots(std::vector<SpotPtr> *spots)
{
    bool verbose = (FileParser::getKey("VERBOSITY_LEVEL", 0) > 1);

    int xDim = this->image->getXDim();
    int yDim = this->image->getYDim();

    int *mask = new int[xDim * yDim];
    memset(mask, 0, sizeof(int) * xDim * yDim);

    SpotPtr spot = SpotPtr(new Spot(image));

    double minCorrelation = FileParser::getKey("IMAGE_MIN_CORRELATION", 0.7);
    double minResolution = 1 / FileParser::getKey("MIN_INTEGRATED_RESOLUTION", 0.);
    double maxResolution = 1 / FileParser::getKey("MAX_INTEGRATED_RESOLUTION", 1.4);

    std::vector<double> acceptableValue = FileParser::getKey("IMAGE_ACCEPTABLE_COMBO", std::vector<double>());
    double gradient = 0; double offset = 0;

    if (acceptableValue.size())
    {
        gradient = (acceptableValue[1] - minCorrelation) / (acceptableValue[0] - threshold);
        offset = acceptableValue[1] - acceptableValue[0] * gradient;
    }

    if (minResolution != minResolution)
    {
        minResolution = 0;
    }

    if (maxResolution != maxResolution)
    {
        maxResolution = FLT_MAX;
    }

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

            bool success = false;

            if (!acceptableValue.size())
            {
                success = (correlation > minCorrelation);
            }
            else
            {
                double correlThresh = offset + value * gradient;
                success = (correlation > correlThresh);
            }

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
