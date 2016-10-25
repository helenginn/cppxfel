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
    
public:
    Hdf5ManagerCheetahLCLS(std::string newName) : Hdf5ManagerCheetah(newName)
    {
        idAddress = FileParser::getKey("CHEETAH_ID_ADDRESSES",
                                                   std::string("entry_1/data_1/experiment_identifier"));
        dataAddress = FileParser::getKey("CHEETAH_DATA_ADDRESSES", std::string("entry_1/data_1/data"));
        identifiersFromAddress(&imagePaths, idAddress);
    }
    
    static Hdf5ManagerCheetahPtr makeManager(std::string filename);
    
    virtual bool dataForImage(std::string address, void **buffer);
    virtual int hdf5MallocBytesForImage(std::string address, void **buffer);
    virtual size_t bytesPerTypeForImageAddress(std::string address);
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahLCLS__) */
