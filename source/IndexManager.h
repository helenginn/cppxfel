//
//  IndexManager.h
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__IndexManager__
#define __cppxfel__IndexManager__

#include <stdio.h>
#include <vector>
#include "Image.h"

class IndexManager
{
protected:
    std::vector<Image *> images;
    
public:
    Image *getImage(int i)
    {
        return images[i];
    }
    
    IndexManager(std::vector<Image *>images);
};

#endif /* defined(__cppxfel__IndexManager__) */
