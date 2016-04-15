//
//  Hdf5ManagerCheetahSacla.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5ManagerCheetahSacla.h"
#include "misc.h"

std::vector<Hdf5ManagerCheetahSacla *> Hdf5ManagerCheetahSacla::cheetahManagers;

Hdf5ManagerCheetahSacla *Hdf5ManagerCheetahSacla::hdf5ManagerForImage(std::string imageName)
{
    for (int i = 0; i < cheetahManagers.size(); i++)
    {
        Hdf5ManagerCheetahSacla *cheetahManager = cheetahManagers[i];
        
        std::string address = cheetahManager->addressForImage(imageName);
        
        if (address.length() > 0)
        {
            return cheetahManager;
        }
    }
    
    return NULL;
}

std::string Hdf5ManagerCheetahSacla::addressForImage(std::string imageName)
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

size_t Hdf5ManagerCheetahSacla::bytesPerTypeForImageAddress(std::string address)
{
    hid_t dataset = H5Dopen1(handle, address.c_str());
    hid_t type = H5Dget_type(dataset);
    
    size_t sizeType = H5Tget_size(type);
    
    return sizeType;
}

bool Hdf5ManagerCheetahSacla::hdf5MallocBytesForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    
    std::cout << dataAddress << std::endl;
    
    size_t sizeSet = 0;
    size_t sizeType = 0;
    
    try
    {
        // h5sget_simple_extent_npoints_f(space_id, npoints, hdferr)
        hid_t dataset = H5Dopen1(handle, dataAddress.c_str());
        hid_t type = H5Dget_type(dataset);
        hid_t space = H5Dget_space(dataset);
        sizeType = H5Tget_size(type);
        sizeSet = H5Sget_simple_extent_npoints(space);
        
        std::cout << "Size: " << sizeSet << " of a " << sizeType << " byte type." << std::endl;
    }
    catch (std::exception e)
    {
        return false; // failure
    }
    
    *buffer = malloc(sizeSet * sizeType);
    
    return true;
}

bool Hdf5ManagerCheetahSacla::dataForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    
    try
    {
        hid_t dataset = H5Dopen1(handle, dataAddress.c_str());
        hid_t type = H5Dget_type(dataset);
        H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
        
        return true; // success
    }
    catch (std::exception e)
    {
        return false; // failure
    }
}