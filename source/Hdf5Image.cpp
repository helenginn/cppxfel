//
//  Hdf5Image.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5Image.h"
#include "Hdf5ManagerCheetahSacla.h"

void Hdf5Image::loadImage()
{
    Hdf5ManagerCheetahSacla *manager = Hdf5ManagerCheetahSacla::hdf5ManagerForImage(getFilename());
    std::string address = manager->addressForImage(getFilename());
    
    useShortData = (manager->bytesPerTypeForImageAddress(address) == 2);
    
    void *buffer = &data[0];
    
    if (useShortData)
    {
        buffer = &shortData[0];
    }
    
    bool success = manager->hdf5MallocBytesForImage(address, &buffer);
    
    if (success)
    {
        success = manager->dataForImage(address, &buffer);
    }
    
    if (!success)
    {
        Logger::mainLogger->addString("Unable to get data from any HDF5 file");
        sendLog();
    }
}