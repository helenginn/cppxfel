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
#include <vector>
#include "LoggableObject.h"

class Hdf5Manager : public LoggableObject
{
private:
    std::string filename;
    
protected:
    hid_t handle;
    std::string lastComponent(std::string path);
    std::string concatenatePaths(std::string path1, std::string path2);
    
public:
    Hdf5Manager(std::string newName);
    virtual ~Hdf5Manager();
    
    int getSubTypeForIndex(std::string address, int objIdx, H5G_obj_t *type);
    std::vector<H5G_obj_t> getSubGroupTypes(std::string address);
    std::vector<std::string> getSubGroupNames(std::string address);
    
    void groupsWithPrefix(std::vector<std::string> *list, std::string prefix, std::string startAddress = "/");

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
