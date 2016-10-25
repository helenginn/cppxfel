//
//  Hdf5ManagerCheetahSacla.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5ManagerCheetahSacla.h"
#include "FileParser.h"

Hdf5ManagerCheetahPtr Hdf5ManagerCheetahSacla::makeManager(std::string filename)
{
    Hdf5ManagerCheetahSaclaPtr cheetahPtr = Hdf5ManagerCheetahSaclaPtr(new Hdf5ManagerCheetahSacla(filename));
    
    return std::static_pointer_cast<Hdf5ManagerCheetah>(cheetahPtr);
}


bool Hdf5ManagerCheetahSacla::dataForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    return Hdf5Manager::dataForAddress(dataAddress, buffer);
}

double Hdf5ManagerCheetahSacla::wavelengthForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "photon_wavelength_A");
    return Hdf5Manager::dataForAddress(dataAddress, buffer);
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