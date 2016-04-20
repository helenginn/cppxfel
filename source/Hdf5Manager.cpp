//
//  Hdf5Manager.cpp
//  cppxfel
//
//  Created by Helen Ginn on 14/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5Manager.h"
#include "FileReader.h"
#include "Logger.h"

#define MAX_NAME 100

void Hdf5Manager::groupsWithPrefix(std::vector<std::string> *list, std::string prefix, std::string startAddress)
{
    std::vector<std::string> names = getSubGroupNames(startAddress);
    
    for (int i = 0; i < names.size(); i++)
    {
        std::string last = lastComponent(names[i]);
        
        if (last.substr(0, prefix.length()) == prefix)
        {
            list->push_back(names[i]);
        }
        else
        {
            H5G_obj_t type;
            getSubTypeForIndex(startAddress, i, &type);
            
            if (type == H5G_GROUP)
            {
                std::cout << "Searching " << names[i] << std::endl;
                groupsWithPrefix(list, prefix, names[i]);
            }
        }
    }
    
    
}
                
int Hdf5Manager::getSubTypeForIndex(std::string address, int objIdx, H5G_obj_t *type)
{
    try
    {
        hid_t group = H5Gopen1(handle, address.c_str());
        *type = H5Gget_objtype_by_idx(group, objIdx);
        H5Gclose(group);
        return 1;
    }
    catch (std::exception e)
    {
        return 0;
    }
}

std::vector<H5G_obj_t> Hdf5Manager::getSubGroupTypes(std::string address)
{
    std::vector<H5G_obj_t> types;
    
    try
    {
        hid_t groupAll = H5Gopen1(handle, address.c_str());
        hsize_t objectNum = 0;
        H5Gget_num_objs(groupAll, &objectNum);
        
        for (int i = 0; i < objectNum; i++)
        {
            H5G_obj_t obj = H5Gget_objtype_by_idx(groupAll, i);
            types.push_back(obj);
        }
        
        H5Gclose(groupAll);
    }
    catch (std::exception e)
    {
        // just return something empty
    }
    
    return types;
}

std::vector<std::string> Hdf5Manager::getSubGroupNames(std::string address)
{
    hid_t groupAll = H5Gopen1(handle, address.c_str());
    hsize_t objectNum = 0;
    H5Gget_num_objs(groupAll, &objectNum);
    
    std::vector<std::string> names;
    
    for (int i = 0; i < objectNum; i++)
    {
        char objectName[MAX_NAME];
        H5Gget_objname_by_idx(groupAll, i, objectName, MAX_NAME);
        std::string concatenated = concatenatePaths(address, objectName);
        names.push_back(concatenated);
    }
    
    H5Gclose(groupAll);
    
    return names;
}

std::string Hdf5Manager::lastComponent(std::string path)
{
    std::vector<std::string> components = FileReader::split(path, "/");
    
    if (components[components.size() - 1] == "")
    {
        components.erase(components.end());
    }
    
    return components[components.size() - 1];
}

std::string Hdf5Manager::concatenatePaths(std::string path1, std::string path2)
{
    if (path1[path1.length() - 1] == '/')
    {
        path1 += path2;
    }
    else
    {
        path1 += "/";
        path1 += path2;
    }
    
    return path1;
}


Hdf5Manager::Hdf5Manager(std::string newName)
{
    filename = newName;
    
    // test code - does this work?
    
    handle = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
}


void Hdf5Manager::closeHdf5()
{
    logged << "Closing HDF5 file " << getFilename() << std::endl;
    sendLog();
    H5Fclose(handle);
}

Hdf5Manager::~Hdf5Manager()
{
    
}