//
//  Hdf5Table.cpp
//  cppxfel
//
//  Created by Helen Ginn on 22/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5Table.h"
#include "Hdf5ManagerProcessing.h"

bool Hdf5Table::writeToManager(Hdf5ManagerProcessingPtr manager, std::string address)
{
    return manager->writeTable(*this, address);
}

int Hdf5Table::readSizeFromManager(Hdf5ManagerProcessingPtr manager, std::string address)
{
    return manager->readSizeForTable(*this, address);
}

int Hdf5Table::readFromManager(Hdf5ManagerProcessingPtr manager, std::string address, void *data)
{
    return manager->recordsForTable(*this, address, data);
}

Hdf5Table::~Hdf5Table()
{
    if (fieldOffsets != NULL)
        free(fieldOffsets);
    fieldOffsets = NULL;
    
    if (fieldSizes != NULL)
        free(fieldSizes);
    fieldSizes = NULL;
    
    if (fieldTypes != NULL)
        free(fieldTypes);
    fieldTypes = NULL;

}