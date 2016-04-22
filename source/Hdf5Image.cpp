//
//  Hdf5Image.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5Table.h"
#include "Hdf5Image.h"
#include "FileParser.h"
#include "Hdf5ManagerCheetahSacla.h"
#include "Hdf5ManagerProcessing.h"
#include <fstream>

typedef struct
{
    double x;
    double y;
    double intensity;
    double sigma;
} Hdf5Spot;

#define HDF5SPOT_FIELD_COUNT 4

void Hdf5Image::failureMessage()
{
    logged << "(" << getBasename() << ") Unable to open image from HDF5 file." << std::endl;
    sendLog();
}

std::string Hdf5Image::findAddress()
{
    std::string address = imageAddress;
    
    if (!address.length())
    {
        Hdf5ManagerCheetahSaclaPtr manager = getManager();
        
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
    Hdf5ManagerCheetahSaclaPtr manager = getManager();
    
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

Hdf5ManagerCheetahSaclaPtr Hdf5Image::getManager()
{
    if (!chManager)
    {
        chManager = Hdf5ManagerCheetahSacla::hdf5ManagerForImage(getFilename());
    }

    return chManager;
}

void Hdf5Image::writeSpotsList(std::string spotFile)
{
    logged << "Writing spot list to HDF5." << std::endl;
    
    Hdf5ManagerProcessingPtr processingManager = Hdf5ManagerProcessing::getProcessingManager();
    
    std::string imageAddress = getAddress();
    std::string spotsName = "spots_v1";
    
    std::string concat = Hdf5Manager::concatenatePaths(imageAddress, spotsName);
    
    int spotCount = this->spotCount();
    
    Hdf5Spot *spotData = (Hdf5Spot *)malloc(sizeof(Hdf5Spot) * spotCount);
    
    for (int i = 0; i < spotCount; i++)
    {
        spotData[i].intensity = 0; // write later if necessary
        spotData[i].sigma = 0; // write later if necessary
        spotData[i].x = spot(i)->getRawXY().first;
        spotData[i].y = spot(i)->getRawXY().second;
    }

    size_t *fieldOffsets = (size_t *)malloc(sizeof(size_t) * HDF5SPOT_FIELD_COUNT);
    fieldOffsets[0] = HOFFSET(Hdf5Spot, x);
    fieldOffsets[1] = HOFFSET(Hdf5Spot, y);
    fieldOffsets[2] = HOFFSET(Hdf5Spot, intensity);
    fieldOffsets[3] = HOFFSET(Hdf5Spot, sigma);
    

    hid_t *fieldTypes = (hid_t *)malloc(sizeof(hid_t) * HDF5SPOT_FIELD_COUNT);
    fieldTypes[0] = H5T_NATIVE_DOUBLE;
    fieldTypes[1] = H5T_NATIVE_DOUBLE;
    fieldTypes[2] = H5T_NATIVE_DOUBLE;
    fieldTypes[3] = H5T_NATIVE_DOUBLE;
    
    size_t *fieldSizes = (size_t *)malloc(sizeof(size_t) * HDF5SPOT_FIELD_COUNT);
    fieldSizes[0] = member_size(Hdf5Spot, x);
    fieldSizes[1] = member_size(Hdf5Spot, y);
    fieldSizes[2] = member_size(Hdf5Spot, intensity);
    fieldSizes[3] = member_size(Hdf5Spot, sigma);
    
    char **fieldNames = (char **)malloc(sizeof(char **) * HDF5SPOT_FIELD_COUNT);
    
    const char *toCopy[HDF5SPOT_FIELD_COUNT] =
    {
        "x",
        "y",
        "intensity",
        "sigma"
    };
    
    for (int i = 0; i < HDF5SPOT_FIELD_COUNT; i++)
    {
        fieldNames[i] = (char *)malloc(sizeof(char*) * strlen(toCopy[i]));
        strcpy(fieldNames[i], toCopy[i]);
    }
    
    Hdf5Table table;
    table.setTableTitle(getBasename() + "_spots_v1");
    table.setTableName(getBasename() + "_spots_v1");
    table.setHeaders((const char**)fieldNames);
    table.setRecordSize(sizeof(Hdf5Spot));
    table.setOffsets(fieldOffsets);
    table.setNumberOfRecords(Image::spotCount());
    table.setNumberOfFields(HDF5SPOT_FIELD_COUNT);
    table.setFieldSizes(fieldSizes);
    table.setTypes(fieldTypes);
    table.setData(spotData);
    table.setCompress(false);
    
    bool success = table.writeToManager(processingManager, imageAddress);
    
    free(fieldOffsets);
    free(fieldSizes);
    free(fieldTypes);
    free(spotData);
    
    logged << "Writing spots to HDF5 file was " << (success ? "successful." : "a failure.") << std::endl;
    
    sendLog();
}
