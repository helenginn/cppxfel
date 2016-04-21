//
//  Hdf5ManagerImageAddresses.h
//  cppxfel
//
//  Created by Helen Ginn on 21/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5ManagerImageAddresses__
#define __cppxfel__Hdf5ManagerImageAddresses__

#include <stdio.h>
#include "Hdf5Manager.h"

class Hdf5ManagerImageAddresses : public Hdf5Manager
{
protected:
    std::vector<std::string> imagePaths;
    
    Hdf5ManagerImageAddresses(std::string newName, Hdf5AccessType access = Hdf5AccessTypeReadOnly) : Hdf5Manager(newName, access)
    {
        groupsWithPrefix(&imagePaths, "tag");
    }

public:
    std::string addressForImage(std::string imageName);
    
};

#endif /* defined(__cppxfel__Hdf5ManagerImageAddresses__) */

