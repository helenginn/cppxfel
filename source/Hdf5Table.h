//
//  Hdf5Table.h
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
