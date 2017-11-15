//
//  SpotFinderQuick.h
//  cppxfel
//
//  Created by Helen Ginn on 29/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__SpotFinderQuick__
#define __cppxfel__SpotFinderQuick__

#include <stdio.h>
#include "SpotFinder.h"

class SpotFinderQuick : public SpotFinder
{
private:
    int minRadius;
    int maxRadius;
    int minPixels;
    int maxPixels;
    int maxHits;
    int minSeparation;
    float signalToNoiseThreshold;

    std::vector<int> backgroundShifts;

    template<class Value>
    void findSignalToNoise(Value *data, size_t position, int xDim, int yDim, float *signalToNoiseRatio, float *background, float *backgroundVariance);
    void calculateBackgroundShifts();
public:
    SpotFinderQuick(ImagePtr image) : SpotFinder(image)
    {
        // this has not been calibrated

        minRadius = 6;
        maxRadius = 8;
        minPixels = FileParser::getKey("SPOT_FINDING_MIN_PIXELS", 2);
        minSeparation = 3;
        maxPixels = FileParser::getKey("SPOT_FINDING_MAX_PIXELS", 40);
        maxHits = 500;

        // For cheetah, Takanori Nakane says:
        // "In LCLS, 4. For SACLA, 4 leads to many false positives"
        signalToNoiseThreshold = FileParser::getKey("SPOT_FINDING_SIGNAL_TO_NOISE", 5.);

        calculateBackgroundShifts();
    }

    virtual void findSpecificSpots(std::vector<SpotPtr> *spots);
};

#endif /* defined(__cppxfel__SpotFinderQuick__) */
