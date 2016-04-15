//
//  Hdf5Manager.cpp
//  cppxfel
//
//  Created by Helen Ginn on 14/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5Manager.h"


Hdf5Manager::Hdf5Manager(std::string newName)
{
    filename = newName;
    
    // test code - does this work?
    
    handle = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    hsize_t objectNum = 0;
    char objectName[50];
    
    hid_t groupAll = H5Gopen1(handle, "/");
    H5Gget_num_objs(groupAll, &objectNum);
    H5Gget_objname_by_idx(groupAll, 0, objectName, 50);
    
    groupAll = groupAll;
}

