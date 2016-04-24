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
    
    void createMillerTable();
    void writeCrystalData(std::string address);
    void writeReflectionData(std::string address);
public:
    Hdf5Crystal(std::string filename) : MtzManager()
    {
        this->filename = filename;
        createMillerTable();
    };
    
    void writeToFile(std::string newFilename = "", bool announce = false, bool shifts = false, bool includeAmbiguity = false, bool useCountingSigma = false);
};

#endif /* defined(__cppxfel__Hdf5Crystal__) */
