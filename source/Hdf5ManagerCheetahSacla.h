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
#include "Hdf5ManagerCheetah.h"
#include <iostream>
#include <mutex>

class Hdf5ManagerCheetahSacla : public Hdf5ManagerCheetah
{
private:    
public:
    static Hdf5ManagerCheetahPtr makeManager(std::string filename);
    
    virtual bool dataForImage(std::string address, void **buffer);
    virtual ~Hdf5ManagerCheetahSacla() {};
    virtual double wavelengthForImage(std::string address, void **buffer);
    virtual int hdf5MallocBytesForImage(std::string address, void **buffer);
    virtual bool getImageSize(std::string address, int *dims);
   
    size_t bytesPerTypeForImageAddress(std::string address);
    
    Hdf5ManagerCheetahSacla(std::string newName) : Hdf5ManagerCheetah(newName)
    {
        groupsWithPrefix(&imagePaths, "tag");
    }
    
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahSacla__) */
