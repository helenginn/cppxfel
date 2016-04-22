//
//  Hdf5ManagerProcessing.cpp
//  cppxfel
//
//  Created by Helen Ginn on 21/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5ManagerProcessing.h"

Hdf5ManagerProcessingPtr Hdf5ManagerProcessing::processingManager = Hdf5ManagerProcessingPtr();


void Hdf5ManagerProcessing::setupProcessingManager()
{
    std::string filename = FileParser::getKey("HDF5_OUTPUT_FILE", std::string("processing.h5"));
    
    processingManager = Hdf5ManagerProcessingPtr(new Hdf5ManagerProcessing(filename));
}

Hdf5ManagerProcessingPtr Hdf5ManagerProcessing::getProcessingManager()
{
    if (!processingManager)
    {
        setupProcessingManager();
    }
    
    return processingManager;
}
