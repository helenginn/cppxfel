//
//  Hdf5ManagerCheetah.h
//  cppxfel
//
//  Created by Helen Ginn on 09/10/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5ManagerCheetah__
#define __cppxfel__Hdf5ManagerCheetah__

#include <stdio.h>
#include "parameters.h"
#include "Hdf5Manager.h"
#include "FileParser.h"

typedef enum
{
    FreeElectronLaserTypeLCLS = 0,
    FreeElectronLaserTypeSACLA = 1,
    FreeElectronLaserTypeEuropeanXFEL = 2,
    FreeElectronLaserTypeSwissFEL = 3
} FreeElectronLaserType;

class Hdf5ManagerCheetah : public Hdf5Manager
{
protected:
    static std::vector<Hdf5ManagerCheetahPtr> cheetahManagers;
    std::vector<std::string> imagePaths;
    std::map<std::string, int> imagePathMap;
    static std::mutex readingPaths;
    
    static std::string maskAddress;
    
public:
    Hdf5ManagerCheetah(std::string newName,
                       Hdf5AccessType accessType = Hdf5AccessTypeReadOnly) : Hdf5Manager(newName, accessType)
    {
        maskAddress = FileParser::getKey("HDF5_MASK_ADDRESS", std::string("/entry_1/instrument_1/detector_1/mask_shared"));
    };
    
    static Hdf5ManagerCheetahPtr hdf5ManagerForImage(std::string imageName);

    static void initialiseCheetahManagers();
    static void closeHdf5Files();
    
    virtual double wavelengthForImage(std::string address, void **buffer) { return 0; };
    virtual bool dataForImage(std::string address, void **buffer, bool rawAddress = false) { return false; };
    virtual int hdf5MallocBytesForImage(std::string address, void **buffer) { return 0; };
    virtual size_t bytesPerTypeForImageAddress(std::string address) { return 0; };
    
    static std::string getMaskAddress()
    {
        return maskAddress;
    }
    
    static int cheetahManagerCount()
    {
        return (int)cheetahManagers.size();
    }
    
    static Hdf5ManagerCheetahPtr cheetahManager(int i)
    {
        return cheetahManagers[i];
    }
    
    // stuff from imageaddresses
    
    std::string addressForImage(std::string imageName);
    
    int imageAddressCount()
    {
        return (int)imagePaths.size();
    }
    
	int numberForAddress(std::string address)
    {
        if (imagePathMap.count(address) == 0)
		{
			return -1;
		}

		return imagePathMap[address];
    }

	std::string imageAddress(int i)
	{
		return imagePaths[i];
	}
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetah__) */
