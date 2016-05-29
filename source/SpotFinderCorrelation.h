//
//  SpotFinderCorrelation.h
//  cppxfel
//
//  Created by Helen Ginn on 29/05/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__SpotFinderCorrelation__
#define __cppxfel__SpotFinderCorrelation__

#include <stdio.h>
#include "SpotFinder.h"
#include "FileParser.h"

class SpotFinderCorrelation : public SpotFinder
{
private:
    double jump;
    
    void focusOnMax(int *x, int *y);

public:
    SpotFinderCorrelation(ImagePtr image) : SpotFinder(image)
    {
        jump = FileParser::getKey("IMAGE_PIXEL_JUMP", 10);
        
        
    }
    
    virtual void findSpecificSpots(std::vector<SpotPtr> *spots);
};

#endif /* defined(__cppxfel__SpotFinderCorrelation__) */
