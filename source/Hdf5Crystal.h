//
//  Hdf5Crystal.h
//  cppxfel
//
//  Created by Helen Ginn on 24/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5Crystal__
#define __cppxfel__Hdf5Crystal__

#include <stdio.h>
#include "MtzManager.h"
#include "Hdf5Table.h"

class Hdf5Crystal : public MtzManager
{
private:
    Hdf5Table millerTable;
    std::string address;
    
    void createMillerTable();
    void writeCrystalData(std::string address);
    void writeReflectionData(std::string address);
public:
    Hdf5Crystal(std::string _filename) : MtzManager()
    {
        setFilename(_filename);
        createMillerTable();
    };
    
    void setAddress(std::string newAddress)
    {
        address = newAddress;
    }
    
    void writeToFile(std::string newFilename = "", bool announce = false);
    
    void loadReflections(PartialityModel model, bool special = false);
};

#endif /* defined(__cppxfel__Hdf5Crystal__) */
