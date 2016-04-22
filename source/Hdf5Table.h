//
//  Hdf5Table.h
//  cppxfel
//
//  Created by Helen Ginn on 22/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5Table__
#define __cppxfel__Hdf5Table__

#include <stdio.h>
#include "parameters.h"
#include "Hdf5Manager.h"
#include <cstdarg>

class Hdf5Table
{
private:
    std::string tableTitle;
    std::string tableName;
    int nfields;
    int nrecords;
    size_t type_size;
    const char **fieldNames;
    size_t *fieldOffsets;
    size_t *fieldSizes;
    hid_t *fieldTypes;
    int chunk_size;
    int *fill_data;
    int compress;
    void *data;
public:
    ~Hdf5Table();
    Hdf5Table()
    {
        fill_data = NULL;
        chunk_size = 10;
        compress = false;
        data = NULL;
        fieldTypes = NULL;
        fieldOffsets = NULL;
        fieldNames = NULL;
        type_size = 0;
        nrecords = 0;
        nfields = 0;
    };
    
    void setTableTitle(std::string title)
    {
        tableTitle = title;
    }
    
    const char *getTableTitle()
    {
        return tableTitle.c_str();
    }
    
    void setTableName(std::string name)
    {
        tableName = name;
    }
    
    const char *getTableName()
    {
        return tableName.c_str();
    }
    
    void setNumberOfFields(int fieldNum)
    {
        nfields = fieldNum;
    }
    
    int getNumberOfFields()
    {
        return nfields;
    }
    
    void setNumberOfRecords(int recordNum)
    {
        nrecords = recordNum;
    }
    
    int getNumberOfRecords()
    {
        return nrecords;
    }
    
    void setRecordSize(int newSize)
    {
        type_size = newSize;
    }
    
    int getRecordSize()
    {
        return (int)type_size;
    }
    
    void setOffsets(size_t *offsets)
    {
        fieldOffsets = offsets;
    }
    
    size_t *getOffsets()
    {
        return fieldOffsets;
    }
    
    void setHeaders(const char **pointerToCharArrayArray)
    {
        fieldNames = pointerToCharArrayArray;
    }
    
    const char **getHeaders()
    {
        return fieldNames;
    }
    
    void setTypes(hid_t *hidTTypes)
    {
        fieldTypes = hidTTypes;
    }
    
    hid_t *getTypes()
    {
        return fieldTypes;
    }
    
    void setCompress(bool shouldCompress)
    {
        compress = shouldCompress;
    }
    
    int getCompress()
    {
        return compress;
    }
    
    int getChunkSize()
    {
        return chunk_size;
    }
    
    void setData(void *buffer)
    {
        data = buffer;
    }
    
    void *getData()
    {
        return data;
    }
    
    int *getFillData()
    {
        return fill_data;
    }
    
    void setFieldSizes(size_t *sizes)
    {
        fieldSizes = sizes;
    }
    
    size_t *getFieldSizes()
    {
        return fieldSizes;
    }

    bool writeToManager(Hdf5ManagerProcessingPtr manager, std::string address);
    int readFromManager(Hdf5ManagerProcessingPtr manager, std::string address, void *data);
    int readSizeFromManager(Hdf5ManagerProcessingPtr manager, std::string address);
};

#endif /* defined(__cppxfel__Hdf5Table__) */
