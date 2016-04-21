//
//  Hdf5ManagerProcessing.h
//  cppxfel
//
//  Created by Helen Ginn on 21/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5ManagerProcessing__
#define __cppxfel__Hdf5ManagerProcessing__

#include <stdio.h>
#include "Hdf5ManagerImageAddresses.h"

class Hdf5ManagerProcessing  : public Hdf5ManagerImageAddresses
{
private:
    
    
public:
    Hdf5ManagerProcessing(std::string filename) : Hdf5ManagerImageAddresses(filename, Hdf5AccessTypeReadWrite)
    {
        
    }
    
    
};

#endif /* defined(__cppxfel__Hdf5ManagerProcessing__) */
