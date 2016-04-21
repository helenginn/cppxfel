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
#include "Hdf5ManagerCheetahSacla.h"

class Hdf5Image : public Image
{
private:
    void failureMessage();
    virtual void loadImage();
    std::string imageAddress;
    
public:
    Hdf5Image(std::string filename = "", double wavelength = 0,
              double distance = 0) : Image(filename, wavelength, distance)
    {
        imageAddress = std::string();
    };
  
    virtual ~Hdf5Image() {};
};

#endif /* defined(__cppxfel__Hdf5Image__) */
