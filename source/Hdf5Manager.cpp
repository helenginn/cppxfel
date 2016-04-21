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

void Hdf5Manager::turnOffErrors()
{
    /* Save old error handler */
//    hid_t error_stack = H5Eget_current_stack();
    
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    
    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    
//    H5Eset_current_stack(error_stack);
}

void Hdf5Manager::turnOnErrors()
{
    hid_t error_stack = H5Eget_current_stack();
    H5Eset_auto(error_stack, old_func, old_client_data);
}

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
        if (group >= 0)
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
        
        if (groupAll >= 0)
            H5Gclose(groupAll);
    }
    catch (std::exception e)
    {
        // just return something empty
    }
    
    return types;
}

bool Hdf5Manager::createLastComponentAsGroup(std::string address)
{
    int components = pathComponentCount(address);
    std::string parentAddress = truncatePath(address, components - 1);
    std::string newGroupName = lastComponent(address);
 
    hid_t groupAll = H5Gopen1(handle, address.c_str());
    
    if (groupAll >= 0)
    {
        H5Gcreate1(groupAll, newGroupName.c_str(), 0);
        
        H5Gclose(groupAll);
        
        return true;
    }
    else
    {
        return false;
    }
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
    
    if (groupAll >= 0)
        H5Gclose(groupAll);
    
    return names;
}

bool Hdf5Manager::groupExists(std::string address)
{
    turnOffErrors();
    
    bool success = false;
    
    hid_t group = H5Gopen1(handle, address.c_str());
    
    if (group >= 0)
    {
        H5Gclose(group);
        success = true;
    }
    
    turnOnErrors();
    
    return success;
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

int Hdf5Manager::pathComponentCount(std::string path)
{
    std::vector<std::string> components = FileReader::split(path, "/");
    
    return (int)components.size();
}

std::string Hdf5Manager::truncatePath(std::string path, int numToTruncate)
{
    std::vector<std::string> components = FileReader::split(path, "/");
    std::string added;
    
    for (int i = 0; i < numToTruncate && i < components.size(); i++)
    {
        added += components[i] + "/";
    }
    
    added.erase(added.end() - 1);
    
    return added;
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

bool Hdf5Manager::createGroupsFromAddress(std::string address)
{
    int componentCount = pathComponentCount(address);
    
    for (int i = 1; i < componentCount; i++)
    {
        std::string pathToCheck = truncatePath(address, i);
        
        bool exists = groupExists(pathToCheck);
        
        if (exists)
        {
            logged << "Address " << pathToCheck << " exists." << std::endl;
        }
        else
        {
            bool success = createLastComponentAsGroup(pathToCheck);
            
            if (success)
            {
                logged << "Address " << pathToCheck << " successfully created." << std::endl;
            }
            else
            {
                logged << "Address " << pathToCheck << " creation failed." << std::endl;
            }
        }
    }
 
    sendLog();

    return true;
}

Hdf5Manager::Hdf5Manager(std::string newName, Hdf5AccessType hdf5AccessType)
{
    filename = newName;
    accessType = hdf5AccessType;
    
    unsigned int accessFlag = H5F_ACC_RDONLY;

    if (hdf5AccessType == Hdf5AccessTypeReadWrite)
    {
        accessFlag = H5F_ACC_RDWR;
    }
    
    turnOffErrors();
    
    handle = H5Fopen(filename.c_str(), accessFlag, H5P_DEFAULT);
    
    if (handle < 0 && accessFlag != H5F_ACC_RDONLY)
    {
        handle = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    
    turnOnErrors();
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