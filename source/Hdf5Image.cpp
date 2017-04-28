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
#include "Hdf5Crystal.h"

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
        Hdf5ManagerCheetahPtr manager = getManager();
        
        if (!manager)
        {
            failureMessage();
            return "";
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

void Hdf5Image::loadCrystals()
{
    Hdf5ManagerProcessingPtr processingManager = Hdf5ManagerProcessing::getProcessingManager();
    
    if (!processingManager)
        return;
    
    std::string address = findAddress();
    
    if (!processingManager->groupExists(address))
    {
        return;
    }
    
    std::vector<std::string> addresses = processingManager->getSubGroupNames(address);
    
    for (int i = 0; i < addresses.size(); i++)
    {
        std::string lastComponent = Hdf5Manager::lastComponent(addresses[i]);
        
        if (lastComponent.substr(0, 3) == "img")
        {
            Hdf5CrystalPtr crystal = Hdf5CrystalPtr(new Hdf5Crystal(lastComponent));
            crystal->setImage(shared_from_this());
            crystal->setAddress(addresses[i]);
            
            MtzPtr mtz = boost::static_pointer_cast<MtzManager>(crystal);
            mtz->loadReflections();
            addMtz(mtz);
        }
    }
}

void Hdf5Image::loadImage()
{
    if (isLoaded())
    {
        return;
    }
    
    bool asFloat = FileParser::getKey("HDF5_AS_FLOAT", false);
    Hdf5ManagerCheetahPtr manager = getManager();
    
    if (!manager)
    {
        failureMessage();
    }
    
    std::string address = imageAddress;

    if (!address.length())
    {
		address = manager->addressForImage(getFilename());
    }
    
    if (!address.length())
    {
        failureMessage();
    }

	logged << "Chosen manager " << manager->getFilename() << std::endl;
	sendLog();

/*
	if (!manager->datasetExists(address) && !manager->groupExists(address))
	{
		logged << "Cannot find an image with data address " << getBasename() << ", nevermind." << std::endl;
		sendLog();
		return;
	}*/

    useShortData = (manager->bytesPerTypeForImageAddress(address) == 2);
    
    char *buffer;
    
    int size = manager->hdf5MallocBytesForImage(address, (void **)&buffer);
    
    if (size > 0)
    {
        int dims[2];
        bool success = manager->getImageSize(address, dims);
        
        xDim = dims[1];
        yDim = dims[0];
        
        if (!success)
        {
            failureMessage();
        }
        
        success = manager->dataForImage(address, (void **)&buffer, _isMask);
        
        if (!success)
        {
            failureMessage();
        }
        
        if (asFloat)
        {
            useShortData = false;
            data.reserve(size / sizeof(int));
            for (int i = 0; i < size; i += sizeof(float))
            {
                float value = 0;
                memcpy(&value, &buffer[i], sizeof(float));
                data.push_back(value);
            }
        }
        else if (!useShortData)
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
    
    free(buffer);
    
    logged << "Loaded data for " << getFilename() << std::endl;
    sendLog();
    
    bool dumpImages = FileParser::getKey("DUMP_IMAGES", false);
    
    if (dumpImages)
    {
        dumpImage();
    }
    
    int totalPixels = xDim * yDim;
    
    overlapMask = vector<signed char>(totalPixels, 0);
    
    checkAndSetupLookupTable();
}

Hdf5ManagerCheetahPtr Hdf5Image::getManager()
{
    if (!chManager)
    {
        chManager = Hdf5ManagerCheetah::hdf5ManagerForImage(getFilename());
    }
    
    if (!chManager && Hdf5ManagerCheetah::cheetahManagerCount() > 0)
    {
        chManager = Hdf5ManagerCheetah::cheetahManager(0);
    }

    return chManager;
}

void Hdf5Image::writeSpotsList(std::string spotFile)
{
    Hdf5ManagerProcessingPtr processingManager = Hdf5ManagerProcessing::getProcessingManager();
    
    if (!processingManager)
    {
        Image::writeSpotsList(spotFile);
        return;
    }
    
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
    
    spotTable.setData(spotData);
    spotTable.setNumberOfRecords(Image::spotCount());

    bool success = spotTable.writeToManager(processingManager, imageAddress);
    
    logged << "Writing spots to HDF5 file was " << (success ? "successful." : "a failure.") << std::endl;
    
    free(spotData);
    
    sendLog();
}

void Hdf5Image::processSpotList()
{
    if (loadedSpots)
    {
        return;
    }
    
    Hdf5ManagerProcessingPtr processingManager = Hdf5ManagerProcessing::getProcessingManager();
    
    if (!processingManager)
    {
        Image::processSpotList();
        return;
    }
    
    std::string imageAddress = getAddress();
    std::string spotsName = getBasename() + "_spots_v1";
    std::string concat = Hdf5Manager::concatenatePaths(imageAddress, spotsName);

    Hdf5Spot *spotData;
    
    bool forceSpotFinding = FileParser::getKey("FORCE_SPOT_FINDING", false);
    
    loadedSpots = true;

    if (!forceSpotFinding && processingManager->datasetExists(concat))
    {
        spots.clear();
        
        int size = spotTable.readSizeFromManager(processingManager, concat);
        spotData = (Hdf5Spot *)malloc(size);

        bool success = spotTable.readFromManager(processingManager, concat, (void *)spotData);
        int count = 0;
        
        logged << "Reading spots from HDF5 file was " << (success ? "successful." : "a failure.") << std::endl;
        sendLog();
        
        for (int i = 0; i < size; i += sizeof(Hdf5Spot))
        {
            SpotPtr newSpot = SpotPtr(new Spot(shared_from_this()));
            newSpot->setXY(spotData[count].x, spotData[count].y);
            
            addSpotIfNotMasked(newSpot);
            
            count++;
        }
        
        if (size > 0)
        {
            free(spotData);
        }
        
        return;
    }

    findSpots();
}

void Hdf5Image::createSpotTable()
{
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
    
    spotTable.setTableTitle(getBasename() + "_spots_v1");
    spotTable.setTableName(getBasename() + "_spots_v1");
    spotTable.setHeaders((const char**)fieldNames);
    spotTable.setRecordSize(sizeof(Hdf5Spot));
    spotTable.setOffsets(fieldOffsets);
    spotTable.setNumberOfFields(HDF5SPOT_FIELD_COUNT);
    spotTable.setFieldSizes(fieldSizes);
    spotTable.setTypes(fieldTypes);
    spotTable.setCompress(false);
}

void Hdf5Image::getWavelengthFromHdf5()
{
    bool useHdf5Wavelength = FileParser::getKey("USE_HDF5_WAVELENGTH", true);
    
    if (!useHdf5Wavelength)
    {
        return;
    }
    
    std::string address = getAddress();

    double wavelength = 0;
    double *wavePtr = &wavelength;

    Hdf5ManagerCheetahPtr manager = this->getManager();
    
    if (!manager)
    {
        return;
    }
/*
	if (!manager->groupExists(address) && !manager->datasetExists(address))
	{
		return;
	}
*/
    manager->wavelengthForImage(address, (void **)&wavePtr);
    
    if (wavelength > 0)
    {
        this->setWavelength(wavelength);
    }
}

int Hdf5Image::getFiducial()
{
	std::string address = "LCLS/fiducial";
	int fiducial = 0;
	int *fidPtr = &fiducial;

	Hdf5ManagerCheetahPtr manager = this->getManager();

	if (!manager)
	{
		return 0;
	}

	int index = manager->numberForAddress(imageAddress);

	if (index < 0)
	{
		return -1;
	}
	
	manager->dataForAddress(address, (void **)&fidPtr, index);

	return fiducial;
}