//
//  Hdf5ManagerCheetahLCLS.h
//  cppxfel
//
//  Created by Helen Ginn on 09/10/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5ManagerCheetahLCLS__
#define __cppxfel__Hdf5ManagerCheetahLCLS__

#include <stdio.h>
#include "Hdf5ManagerCheetah.h"
#include "FileParser.h"

class Hdf5ManagerCheetahLCLS : public Hdf5ManagerCheetah
{
private:
    // place where all the stuff is
    // stupid comments, sorry helen.
    // I'm tired.
    std::string idAddress;
    std::string dataAddress;
    std::string wavelengthAddress;
    std::vector<double> wavelengths;
    
public:
    Hdf5ManagerCheetahLCLS(std::string newName) : Hdf5ManagerCheetah(newName)
    {
        idAddress = FileParser::getKey("CHEETAH_ID_ADDRESSES",
                                                   std::string("entry_1/data_1/experiment_identifier"));
        dataAddress = FileParser::getKey("CHEETAH_DATA_ADDRESSES", std::string("entry_1/data_1/data"));
        wavelengthAddress = FileParser::getKey("CHEETAH_WAVELENGTH_ADDRESSES", std::string("LCLS/photon_wavelength_A"));
        identifiersFromAddress(&imagePathMap, &imagePaths, idAddress);
        prepareWavelengths();
    }
    
    static Hdf5ManagerCheetahPtr makeManager(std::string filename);
    
    void prepareWavelengths();
    virtual double wavelengthForImage(std::string address, void **buffer);
    virtual bool dataForImage(std::string address, void **buffer);
    virtual int hdf5MallocBytesForImage(std::string address, void **buffer);
    virtual size_t bytesPerTypeForImageAddress(std::string address);
    virtual bool getImageSize(std::string address, int *dims);
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahLCLS__) */
