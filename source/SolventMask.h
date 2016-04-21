//
//  SolventMask.h
//  cppxfel
//
//  Created by Helen Ginn on 21/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__SolventMask__
#define __cppxfel__SolventMask__

#include <stdio.h>
#include <vector>
#include "parameters.h"

class SolventMask
{
private:
    // to be expressed as 1 / resolution.
    double lowResCutoff;
    double highResCutoff;
    
    template <class ObjectPtr>
    bool thisSolventMasks(ObjectPtr objectPtr)
    {
        double objectResolution = objectPtr->resolution();
        
        return (objectResolution > lowResCutoff && objectResolution < highResCutoff);
    }
    
    static std::vector<SolventMaskPtr> solventMasks;
public:
    SolventMask(double lowRes, double highRes);
   
    // plan to take either Miller or Spot.
    template <class ObjectPtr>
    static bool isMasked(ObjectPtr objectPtr)
    {
        for (int i = 0; i < solventMaskCount(); i++)
        {
            if (solventMasks[i]->thisSolventMasks(objectPtr))
                return true;
        }
        
        return false;
    }
    
    static void addSolventMask(double lowRes, double highRes);

    
    static int solventMaskCount()
    {
        return (int)solventMasks.size();
    }
    
};

#endif /* defined(__cppxfel__SolventMask__) */
