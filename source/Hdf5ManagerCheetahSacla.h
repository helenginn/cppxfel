//
//  Hdf5ManagerCheetahSacla.h
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5ManagerCheetahSacla__
#define __cppxfel__Hdf5ManagerCheetahSacla__

#include <stdio.h>
#include "Hdf5Manager.h"
#include <iostream>

class Hdf5ManagerCheetahSacla : public Hdf5Manager
{
private:
    std::vector<std::string> imagePaths;
    static std::vector<Hdf5ManagerCheetahSaclaPtr> cheetahManagers;

public:
    static void initialiseSaclaManagers();
    
    static Hdf5ManagerCheetahSaclaPtr hdf5ManagerForImage(std::string imageName);
    
    std::string addressForImage(std::string imageName);
    bool dataForImage(std::string address, void **buffer);
    size_t bytesPerTypeForImageAddress(std::string address);
    
    virtual ~Hdf5ManagerCheetahSacla() {};
    
    Hdf5ManagerCheetahSacla(std::string newName) : Hdf5Manager(newName)
    {
        groupsWithPrefix(&imagePaths, "tag");
    }
    
    int hdf5MallocBytesForImage(std::string address, void **buffer);
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahSacla__) */
