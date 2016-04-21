//
//  Hdf5ManagerCheetahSacla.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5ManagerCheetahSacla.h"
#include "misc.h"
#include "FileParser.h"

std::vector<Hdf5ManagerCheetahSaclaPtr> Hdf5ManagerCheetahSacla::cheetahManagers;

void Hdf5ManagerCheetahSacla::initialiseSaclaManagers()
{
    if (cheetahManagers.size() > 0)
        return;
    
    std::vector<std::string> hdf5FileGlobs = FileParser::getKey("HDF5_SOURCE_FILES", std::vector<std::string>());
    std::ostringstream logged;
    
    for (int i = 0; i < hdf5FileGlobs.size(); i++)
    {
        std::vector<std::string> hdf5Files = glob(hdf5FileGlobs[i]);
        
        logged << "HDF5_SOURCE_FILES entry " << i << " matches: " << std::endl;
        
        for (int j = 0; j < hdf5Files.size(); j++)
        {
            std::string aFilename = hdf5Files[i];
            
            Hdf5ManagerCheetahSaclaPtr cheetahPtr = Hdf5ManagerCheetahSaclaPtr(new Hdf5ManagerCheetahSacla(aFilename));
            cheetahManagers.push_back(cheetahPtr);
            
            logged << hdf5Files[j] << ", ";
        }
        
        logged << std::endl;
    }
    
    logged << "... now managing " << cheetahManagers.size() << " hdf5 image source files." << std::endl;
    
    Logger::mainLogger->addStream(&logged);
}

Hdf5ManagerCheetahSaclaPtr Hdf5ManagerCheetahSacla::hdf5ManagerForImage(std::string imageName)
{
    for (int i = 0; i < cheetahManagers.size(); i++)
    {
        Hdf5ManagerCheetahSaclaPtr cheetahManager = cheetahManagers[i];
        
        std::string address = cheetahManager->addressForImage(imageName);
        
        if (address.length() > 0)
        {
            return cheetahManager;
        }
    }
    
    return Hdf5ManagerCheetahSaclaPtr();
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

int Hdf5ManagerCheetahSacla::hdf5MallocBytesForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    
    std::cout << dataAddress << std::endl;
    
    size_t sizeSet = 0;
    size_t sizeType = 0;
    
    readingHdf5.lock();
    
    try
    {
        // h5sget_simple_extent_npoints_f(space_id, npoints, hdferr)
        hid_t dataset = H5Dopen1(handle, dataAddress.c_str());
        hid_t type = H5Dget_type(dataset);
        hid_t space = H5Dget_space(dataset);
        sizeType = H5Tget_size(type);
        sizeSet = H5Sget_simple_extent_npoints(space);
        
        if (space >= 0)
            H5Sclose(space);
        
        if (type >= 0)
            H5Tclose(type);
        
        if (dataset >= 0)
            H5Dclose(dataset);
        
        logged << "Size: " << sizeSet << " of a " << sizeType << " byte type." << std::endl;
        sendLog();
    }
    catch (std::exception e)
    {
        return 0; // failure
    }
    
    readingHdf5.unlock();
    
    *buffer = malloc(sizeSet * sizeType);
    
    return (int)(sizeSet * sizeType);
}

size_t Hdf5ManagerCheetahSacla::bytesPerTypeForImageAddress(std::string address)
{
    std::string dataAddress = concatenatePaths(address, "data");
    
    hid_t dataset, type;
    
    readingHdf5.lock();
    
    try
    {
        dataset = H5Dopen1(handle, dataAddress.c_str());
        type = H5Dget_type(dataset);
        H5Dclose(dataset);
        H5Tclose(type);
    }
    catch (std::exception e)
    {
        return 0; // failure
    }
    
    readingHdf5.unlock();

    size_t sizeType = H5Tget_size(type);
    
    return sizeType;
}



bool Hdf5ManagerCheetahSacla::dataForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    
    readingHdf5.lock();
    
    try
    {
        hid_t dataset = H5Dopen1(handle, dataAddress.c_str());
        hid_t type = H5Dget_type(dataset);
        H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
        H5Dclose(dataset);
        H5Tclose(type);

        readingHdf5.unlock();

        return true; // success
    }
    catch (std::exception e)
    {
        readingHdf5.unlock();

        return false; // failure
    }
}


void Hdf5ManagerCheetahSacla::closeHdf5Files()
{
    for (int i = 0; i < cheetahManagers.size(); i++)
    {
        cheetahManagers[i]->closeHdf5();
    }
}