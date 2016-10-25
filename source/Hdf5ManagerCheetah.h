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

typedef enum
{
    FreeElectronLaserTypeLCLS,
    FreeElectronLaserTypeSACLA,
    FreeElectronLaserTypeEuropeanXFEL,
    FreeElectronLaserTypeSwissFEL
} FreeElectronLaserType;

class Hdf5ManagerCheetah : public Hdf5Manager
{
protected:
    static std::vector<Hdf5ManagerCheetahPtr> cheetahManagers;
    std::vector<std::string> imagePaths;
    static std::mutex readingPaths;
    
public:
    Hdf5ManagerCheetah(std::string newName,
                       Hdf5AccessType accessType = Hdf5AccessTypeReadOnly) : Hdf5Manager(newName, accessType)
    {
    
    };
    
    static Hdf5ManagerCheetahPtr hdf5ManagerForImage(std::string imageName);

    static void initialiseCheetahManagers();
    static void closeHdf5Files();
    
    virtual double wavelengthForImage(std::string address, void **buffer) { return 0; };
    virtual bool dataForImage(std::string address, void **buffer) { return false; };
    virtual int hdf5MallocBytesForImage(std::string address, void **buffer) { return 0; };
    virtual size_t bytesPerTypeForImageAddress(std::string address) { return 0; };
    
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
    
    std::string imageAddress(int i)
    {
        return imagePaths[i];
    }
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetah__) */
