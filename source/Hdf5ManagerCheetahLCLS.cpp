//
//  Hdf5ManagerCheetahLCLS.cpp
//  cppxfel
//
//  Created by Helen Ginn on 09/10/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5ManagerCheetahLCLS.h"
#include <iterator>
#include <algorithm>

Hdf5ManagerCheetahPtr Hdf5ManagerCheetahLCLS::makeManager(std::string filename)
{
    Hdf5ManagerCheetahLCLSPtr cheetahPtr = Hdf5ManagerCheetahLCLSPtr(new Hdf5ManagerCheetahLCLS(filename));
    
    return std::static_pointer_cast<Hdf5ManagerCheetah>(cheetahPtr);
}

bool Hdf5ManagerCheetahLCLS::dataForImage(std::string address, void **buffer)
{
    // need to find the index of the correct image

    // lazy as.
    
    std::vector<std::string>::iterator it = std::find(imagePaths.begin(),
                                                imagePaths.end(),
                                                address);
    
    if (it != imagePaths.end() && *it == address)
    {
        // am a terrible caster. sorry.
        int index = (int)(it - imagePaths.begin());
        
        return Hdf5Manager::dataForAddress(dataAddress, buffer, index);
    }
    
    return false;
}

int Hdf5ManagerCheetahLCLS::hdf5MallocBytesForImage(std::string address, void **buffer)
{
    return Hdf5Manager::hdf5MallocBytesForDataset(dataAddress, buffer);
}

size_t Hdf5ManagerCheetahLCLS::bytesPerTypeForImageAddress(std::string address)
{
    return Hdf5Manager::bytesPerTypeForDatasetAddress(dataAddress);
}