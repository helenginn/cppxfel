//
//  Hdf5Manager.h
//  cppxfel
//
//  Created by Helen Ginn on 14/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5Manager__
#define __cppxfel__Hdf5Manager__

#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string>

class Hdf5Manager
{
private:
    std::string filename;
    hid_t handle;
    
public:
    Hdf5Manager(std::string newName);
    
    void setFilename(std::string newName)
    {
        filename = newName;
    }
    
    std::string getFilename()
    {
        return filename;
    }
};

#endif /* defined(__cppxfel__Hdf5Manager__) */
