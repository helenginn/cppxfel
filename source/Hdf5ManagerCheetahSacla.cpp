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
        
        logged << "HDF5_SOURCE_FILES entry (no. " << i << ") matches: " << std::endl;
        
        for (int j = 0; j < hdf5Files.size(); j++)
        {
            std::string aFilename = hdf5Files[j];
            
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

int Hdf5ManagerCheetahSacla::hdf5MallocBytesForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    return Hdf5Manager::hdf5MallocBytesForDataset(dataAddress, buffer);
}

size_t Hdf5ManagerCheetahSacla::bytesPerTypeForImageAddress(std::string address)
{
    std::string dataAddress = concatenatePaths(address, "data");
    return Hdf5Manager::bytesPerTypeForDatasetAddress(dataAddress);
}

bool Hdf5ManagerCheetahSacla::dataForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    return Hdf5Manager::dataForAddress(dataAddress, buffer);
}

void Hdf5ManagerCheetahSacla::closeHdf5Files()
{
    for (int i = 0; i < cheetahManagers.size(); i++)
    {
        cheetahManagers[i]->closeHdf5();
    }
}

