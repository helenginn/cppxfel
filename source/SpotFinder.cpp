//
//  SpotFinder.cpp
//  cppxfel
//
//  Created by Helen Ginn on 29/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotFinder.h"
#include "Spot.h"

std::vector<SpotPtr> SpotFinder::findSpots()
{
    std::vector<SpotPtr> spots;

    peaks = NULL;
    totalPeaks = 0;

    findSpecificSpots(&spots);

    // now peaks is full of information

    for (int i = 0; i < totalPeaks; i++)
    {
        Peak *peak = &peaks[i];

        SpotPtr spot = SpotPtr(new Spot(image));
        spot->setXY(peak->centreX, peak->centreY);

        spots.push_back(spot);
    }

    free(peaks);
    peaks = NULL;

    return spots;
}
