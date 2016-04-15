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
    static std::vector<Hdf5ManagerCheetahSacla *> cheetahManagers;

public:
    static Hdf5ManagerCheetahSacla *hdf5ManagerForImage(std::string imageName);
    
    std::string addressForImage(std::string imageName);
    bool dataForImage(std::string address, void **buffer);
    size_t bytesPerTypeForImageAddress(std::string address);
    bool hdf5MallocBytesForImage(std::string address, void **buffer);
    
    virtual ~Hdf5ManagerCheetahSacla() {};
    
    Hdf5ManagerCheetahSacla(std::string newName) : Hdf5Manager(newName)
    {
        groupsWithPrefix(&imagePaths, "tag");
        
        cheetahManagers.push_back(this);
        
        /*
        std::string test = imagePaths[0];
        std::string imageName = lastComponent(test);
        
        std::cout << "Examining image name " << imageName << std::endl;
        
        short *buffer;
        hdf5MallocBytesForImage(imageName, (void **)&buffer);
        dataForImage(imageName, (void **)&buffer);
        
        for (int i = 0; i < 10; i++)
        {
            std::cout << buffer[i] << " ";
        }
        
        std::cout << std::endl;*/
    }
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahSacla__) */
