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