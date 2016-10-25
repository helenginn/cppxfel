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

class Hdf5Table;

typedef enum
{
    Hdf5AccessTypeReadOnly,
    Hdf5AccessTypeReadWrite,
} Hdf5AccessType;

class Hdf5Manager : public LoggableObject
{
private:
    // for disabling errors momentarily
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;
    
    std::string filename;
    Hdf5AccessType accessType;
    void turnOffErrors();
    void turnOnErrors();

protected:
    hid_t handle;
    static std::mutex readingHdf5;
    
public:
    static int pathComponentCount(std::string path);
    static std::string concatenatePaths(std::string path1, std::string path2);
    static std::string truncatePath(std::string path, int numToTruncate);
    
    Hdf5Manager(std::string newName, Hdf5AccessType hdf5AccessType = Hdf5AccessTypeReadOnly);
    void closeHdf5();
    virtual ~Hdf5Manager();
    
    static std::string truncateLastComponent(std::string path);
    static std::string lastComponent(std::string path);
    int getSubTypeForIndex(std::string address, int objIdx, H5G_obj_t *type);
    std::vector<H5G_obj_t> getSubGroupTypes(std::string address);
    std::vector<std::string> getSubGroupNames(std::string address);
    bool createGroupsFromAddress(std::string address);
    bool createLastComponentAsGroup(std::string address);
    bool groupExists(std::string address);

    bool datasetExists(std::string address);
    int readSizeForTable(Hdf5Table &table, std::string address);
    bool recordsForTable(Hdf5Table &table, std::string address, void *data);
    bool writeTable(Hdf5Table &table, std::string address);
    
    size_t bytesPerTypeForDatasetAddress(std::string dataAddress);
    int hdf5MallocBytesForDataset(std::string dataAddress, void **buffer);
    bool createDataset(std::string address, int nDimensions, hsize_t *dims, hid_t type);
    bool writeDataset(std::string address, void **buffer, hid_t type);
    bool dataForAddress(std::string address, void **buffer, int offset = -1);
    void identifiersFromAddress(std::vector<std::string> *list, std::string idAddress);
    
    template <class Value>
    bool writeSingleValueDataset(std::string address, Value value, hid_t type)
    {
        hsize_t dim = 1;
        
        Value *valuePtr = &value;
        
        bool success = createDataset(address, 1, &dim, type);
        
        if (success)
        {
            success = writeDataset(address, (void **)valuePtr, type);
        }
        
        return success;
    }
    
    template <class Value>
    bool readDatasetValue(std::string address, Value *value)
    {
        bool success = dataForAddress(address, (void **)&value);
        
        return success;
    }

    
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
