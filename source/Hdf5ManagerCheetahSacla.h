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
#include "Hdf5ManagerImageAddresses.h"
#include <iostream>
#include <mutex>

class Hdf5ManagerCheetahSacla : public Hdf5ManagerImageAddresses
{
private:
    static std::vector<Hdf5ManagerCheetahSaclaPtr> cheetahManagers;

    std::mutex readingHdf5;
    
public:
    static void initialiseSaclaManagers();
    
    static Hdf5ManagerCheetahSaclaPtr hdf5ManagerForImage(std::string imageName);
    
    bool dataForImage(std::string address, void **buffer);
    size_t bytesPerTypeForImageAddress(std::string address);
    
    static void closeHdf5Files();
    virtual ~Hdf5ManagerCheetahSacla() {};
    
    Hdf5ManagerCheetahSacla(std::string newName) : Hdf5ManagerImageAddresses(newName)
    {

    }
    
    int hdf5MallocBytesForImage(std::string address, void **buffer);
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahSacla__) */
