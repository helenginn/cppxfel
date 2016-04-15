//
//  Hdf5Image.h
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5Image__
#define __cppxfel__Hdf5Image__

#include <stdio.h>
#include "Image.h"

class Hdf5Image : public Image
{
private:
    virtual void loadImage();
public:
    
};

#endif /* defined(__cppxfel__Hdf5Image__) */
