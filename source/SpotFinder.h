//
//  SpotFinder.h
//  cppxfel
//
//  Created by Helen Ginn on 29/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__SpotFinder__
#define __cppxfel__SpotFinder__
#include "parameters.h"

#include <stdio.h>

class SpotFinder
{
protected:
    typedef struct
    {
        float centreX;
        float centreY;
        float totalSignal;
    } Peak;

    ImagePtr image;
    Peak *peaks;
    size_t totalPeaks;

public:
    SpotFinder(ImagePtr anImage)
    {
        image = anImage;
        peaks = NULL;
        totalPeaks = 0;
    }
    
    virtual void findSpecificSpots() {};
    std::vector<SpotPtr> findSpots();
};

#endif /* defined(__cppxfel__SpotFinder__) */
