//
//  Hdf5Image.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5Image.h"
#include "FileParser.h"
#include "Hdf5ManagerCheetahSacla.h"
#include <fstream>

void Hdf5Image::failureMessage()
{
    logged << "(" << getBasename() << ") Unable to open image from HDF5 file." << std::endl;
    sendLog();
}

std::string Hdf5Image::findAddress(Hdf5ManagerCheetahSaclaPtr manager)
{
    std::string address = imageAddress;
    
    if (!address.length())
    {
        if (!manager)
        {
            manager = Hdf5ManagerCheetahSacla::hdf5ManagerForImage(getFilename());
        }
        
        address = manager->addressForImage(getFilename());
        imageAddress = address;
    }
    
    if (!address.length())
    {
        failureMessage();
    }
    
    return address;
}

void Hdf5Image::loadImage()
{
    Hdf5ManagerCheetahSaclaPtr manager = Hdf5ManagerCheetahSacla::hdf5ManagerForImage(getFilename());
    
    if (!manager)
        return;
    
    std::string address = imageAddress;
    
    if (!address.length())
    {
        address = manager->addressForImage(getFilename());
    }
    
    if (!address.length())
    {
        failureMessage();
    }
    
    useShortData = (manager->bytesPerTypeForImageAddress(address) == 2);
    
    char *buffer;
    
    int size = manager->hdf5MallocBytesForImage(address, (void **)&buffer);
    
    if (size > 0)
    {
        bool success = manager->dataForImage(address, (void **)&buffer);
        
        if (!success)
        {
            failureMessage();
        }
        
        if (!useShortData)
        {
            data.resize(size / sizeof(int));
            memcpy(&data[0], &buffer[0], size);
        }
        else
        {
            shortData.resize(size / sizeof(short));
            memcpy(&shortData[0], &buffer[0], size);
        }

    }
    else
    {
        Logger::mainLogger->addString("Unable to get data from any HDF5 file");
        sendLog();
    }
    
    bool dumpImages = FileParser::getKey("DUMP_IMAGES", false);
    
    if (dumpImages)
    {
        std::ofstream imgStream;
        imgStream.open(getFilename().c_str(), std::ios::binary);
        
        long int size = useShortData ? shortData.size() : data.size();
        size *= useShortData ? sizeof(short) : sizeof(int);
        char *start = useShortData ? (char *)&shortData[0] : (char *)&data[0];
        
        imgStream.write(start, size);
        
        imgStream.close();
    }
}
