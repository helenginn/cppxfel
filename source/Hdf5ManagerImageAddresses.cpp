//
//  Hdf5ManagerImageAddresses.cpp
//  cppxfel
//
//  Created by Helen Ginn on 21/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5ManagerImageAddresses.h"
#include "misc.h"

std::string Hdf5ManagerImageAddresses::addressForImage(std::string imageName)
{
    // maybe imageName has .img extension, so let's get rid of it
    std::string baseName = getBaseFilename(imageName);
    
    for (int i = 0; i < imagePaths.size(); i++)
    {
        std::string lastName = lastComponent(imagePaths[i]);
        
        if (lastName == baseName)
        {
            return imagePaths[i];
        }
    }
    
    return "";
}