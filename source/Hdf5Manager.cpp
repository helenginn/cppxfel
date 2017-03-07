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
#include "Hdf5Table.h"

#define MAX_NAME 100
std::mutex Hdf5Manager::readingHdf5;

int Hdf5Manager::readSizeForTable(Hdf5Table &table, std::string address)
{
    // herr_t H5TBread_table( hid_t loc_id, const char *table_name, size_t dst_size,  const size_t *dst_offset, const size_t *dst_sizes, void *dst_buf )
    
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    std::string groupOnly = truncateLastComponent(address);
    std::string tableName = lastComponent(address);
    hid_t group = H5Gopen2(handle, groupOnly.c_str(), H5P_DEFAULT);
    
    if (group < 0)
    {
        return 0;
    }
    
    hsize_t nFields = 0;
    hsize_t nRecords = 0;
    size_t recordSize = 0;
    
    // herr_t H5TBget_field_info ( hid_t loc_id, const char *table_name, char *field_names[], size_t *field_sizes, size_t *field_offsets, size_t *type_size  )
    herr_t error = H5TBget_table_info(group, table.getTableName(), &nFields, &nRecords);
    
    if (error < 0)
    {
        H5Gclose(group);
        return 0;
    }
    
    error = H5TBget_field_info(group, table.getTableName(), NULL, NULL, NULL, &recordSize);
    
    H5Gclose(group);

    if (error < 0)
    {
        return 0;
    }
    
    int totalToAllocate = (int)(recordSize * nRecords);
    
//    logged << "Allocating " << totalToAllocate << " bytes to read in table." << std::endl;
//    sendLog();
    
    return totalToAllocate;
}


bool Hdf5Manager::recordsForTable(Hdf5Table &table, std::string address, void *data)
{
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    std::string groupOnly = truncateLastComponent(address);
    hid_t group = H5Gopen2(handle, groupOnly.c_str(), H5P_DEFAULT);
    
    if (group < 0)
    {
        return 0;
    }
    
    herr_t error = H5TBread_table(group, table.getTableName(), table.getRecordSize(), table.getOffsets(), table.getFieldSizes(), data);
    
    H5Gclose(group);
    
    if (error < 0)
    {
        return 0;
    }
    
    return true;
}

bool Hdf5Manager::datasetExists(std::string address)
{
    hid_t table;
    
    std::lock_guard<std::mutex> lg(readingHdf5);
    
 //   logged << "Checking if table " << address << " exists." << std::endl;
 //   sendLog();
    
    turnOffErrors();
    table = H5Dopen2(handle, address.c_str(), H5P_DEFAULT);
    turnOnErrors();
    
    if (table < 0)
    {
        return false;
    }
    
    H5Dclose(table);
    
    return true;
}

bool Hdf5Manager::writeTable(Hdf5Table &table, std::string address)
{
    std::string groupsOnly = address;
    std::string tableAddress = concatenatePaths(address, table.getTableName());
    
    createGroupsFromAddress(address);
    
    readingHdf5.lock();

    hid_t group = H5Gopen1(handle, groupsOnly.c_str());

    readingHdf5.unlock();

    if (group < 0)
    {
        return false;
    }
    
    bool exists = datasetExists(tableAddress);
    
    if (exists)
    {
     //   logged << "Overwriting existing table." << std::endl;
     //   sendLog();
        
        std::lock_guard<std::mutex> lg(readingHdf5);
        
        hsize_t nFields = 0;
        hsize_t nRecords = 0;
        
        herr_t error = H5TBget_table_info(group, table.getTableName(), &nFields, &nRecords);
        hsize_t nRecordsNew = table.getNumberOfRecords();
        
        hsize_t maxRecords = std::min(nRecords, nRecordsNew);
        
        error = H5TBwrite_records(group,
                                      table.getTableName(),
                                         0,
                                      maxRecords,
                                      table.getRecordSize(),
                                      table.getOffsets(),
                                      table.getFieldSizes(),
                                      table.getData());

        if (error < 0)
        {
            logged << "Writing error (maxRecords = " << maxRecords << ") " << nRecords << ", " << nRecordsNew << std::endl;
            sendLog();
        }

        if (maxRecords < nRecordsNew)
        {
            hsize_t extraRecords = nRecordsNew - maxRecords;
            
            error = H5TBappend_records(group,
                                       table.getTableName(),
                                       extraRecords,
                                       table.getRecordSize(),
                                       table.getOffsets(),
                                       table.getFieldSizes(),
                                       table.getData());
            
            error = H5TBget_table_info(group, table.getTableName(), &nFields, &nRecords);
            
            if (error < 0)
            {
                logged << "Appending error" << std::endl;
                sendLog();
            }
        }
        
        if (error < 0)
        {
            H5Gclose(group);
            return false;
        }
        
        hsize_t records_to_delete = nRecords - table.getNumberOfRecords();
        
        if (records_to_delete <= 0)
        {
            H5Gclose(group);
            return true;
        }
        
        error = H5TBdelete_record(group, table.getTableName(), table.getNumberOfRecords(), records_to_delete);
        
        if (error < 0)
        {
            logged << "Deleting error" << std::endl;
            sendLog();
        }

        if (error < 0)
        {
            H5Gclose(group);
            return false;
        }
        
        H5Gclose(group);
        return true;
    }
    
//    logged << "Making new table " << table.getTableTitle() << " in group " << address << std::endl;
//    sendLog();
    
    const char *title = table.getTableTitle();  // ok
    const char *name = table.getTableName();    // ok
    int nrecords = table.getNumberOfRecords();  // ok
    int nfields = table.getNumberOfFields();    // ok
    int recordSize = table.getRecordSize();     // ok
    const char **headers = table.getHeaders();  // ok
    size_t *offsets = table.getOffsets();       // ok
    hid_t *types = table.getTypes();            // ok
    int chunksize = table.getChunkSize();
    void *fill_data = table.getFillData();
    int compress = table.getCompress();
    void *data = table.getData();
    
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    herr_t error = H5TBmake_table(title,
                                  group,
                                  name,
                                  nfields,
                                  nrecords,
                                  recordSize,
                                  headers,
                                  offsets,
                                  types,
                                  chunksize,
                                  fill_data,
                                  compress,
                                  data);
    
    if (error >= 0)
    {
  //      logged << "Written " << nfields << " fields, " << nrecords << " records of " << recordSize << " bytes." << std::endl;
  //      sendLog();
        
        H5Gclose(group);
        return true;
    }
    
    H5Gclose(group);
    return false;
}

bool Hdf5Manager::createDataset(std::string address, int nDimensions, hsize_t *dims, hid_t type)
{
    bool existsAlready = datasetExists(address);
    
    if (existsAlready)
        return true;
    
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    hid_t dataspace_id = dataspace_id = H5Screate_simple(nDimensions, dims, NULL);

    if (dataspace_id < 0)
    {
        return false;
    }
    
    hid_t dataset_id = H5Dcreate2(handle, address.c_str(), type, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    if (dataset_id < 0)
    {
        return false;
    }
    
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    
    return true;
}

bool Hdf5Manager::writeDataset(std::string address, void **buffer, hid_t type)
{
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    hid_t dataset_id = H5Dopen2(handle, address.c_str(), H5P_DEFAULT);
    
    if (dataset_id < 0)
    {
        return false;
    }
    
    /* Write the dataset. */
    herr_t status = H5Dwrite(dataset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      buffer);

    return (status >= 0);
}

int Hdf5Manager::hdf5MallocBytesForDataset(std::string dataAddress, void **buffer)
{
    size_t sizeSet = 0;
    size_t sizeType = 0;
    
    std::lock_guard<std::mutex> lg(readingHdf5);
    size_t numPixels = 0;
    
    try
    {
        // h5sget_simple_extent_npoints_f(space_id, npoints, hdferr)
        hid_t dataset = H5Dopen1(handle, dataAddress.c_str());
        hid_t type = H5Dget_type(dataset);
        hid_t space = H5Dget_space(dataset);
        sizeType = H5Tget_size(type);
        sizeSet = H5Sget_simple_extent_npoints(space);
        
        int numDims = H5Sget_simple_extent_ndims(space);
        
        if (numDims < 3)
        {
            // This is a single image, like a SACLA image
            numPixels = sizeSet;
        }
        else
        {
            // This is probably - hopefully - a pack of images, like at LCLS
            // int H5Sget_simple_extent_dims(hid_t space_id, hsize_t *dims, hsize_t *maxdims )
            hsize_t dims[numDims];
            H5Sget_simple_extent_dims(space, dims, NULL);
            numPixels = dims[1] * dims[2];
        }
        
        if (space >= 0)
            H5Sclose(space);
        
        if (type >= 0)
            H5Tclose(type);
        
        if (dataset >= 0)
            H5Dclose(dataset);
        
    }
    catch (std::exception e)
    {
        return 0; // failure
    }
    
    *buffer = malloc(numPixels * sizeType);
    
    return (int)(numPixels * sizeType);
}

size_t Hdf5Manager::bytesPerTypeForDatasetAddress(std::string dataAddress)
{
    hid_t dataset, type;
    
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    try
    {
        dataset = H5Dopen1(handle, dataAddress.c_str());
        type = H5Dget_type(dataset);
        
        if (dataset >= 0)
        {
            H5Dclose(dataset);
        }
        else
        {
            return 0;
        }
    }
    catch (std::exception e)
    {
        return 0; // failure
    }
    
    size_t sizeType = 0;
    
    if (type >= 0)
    {
        sizeType = H5Tget_size(type);
        H5Tclose(type);
    }
    
    return sizeType;
}

bool Hdf5Manager::getImageSize(std::string dataAddress, int *finalDims)
{
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    try
    {
        hid_t dataset = H5Dopen1(handle, dataAddress.c_str());
        hid_t space = H5Dget_space(dataset);
        int numDims = H5Sget_simple_extent_ndims(space);
        hsize_t dims[numDims];
        H5Sget_simple_extent_dims(space, dims, NULL);
        
        int start = (numDims == 3) ? 1 : 0;
        
        finalDims[0] = (int)dims[start];
        finalDims[1] = (int)dims[start + 1];
    }
    catch (std::exception e)
    {
        return false;
    }
    
    return true;
}

bool Hdf5Manager::dataForAddress(std::string dataAddress, void **buffer, int offset)
{
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    try
    {
        hid_t dataset = H5Dopen1(handle, dataAddress.c_str());
        hid_t type = H5Dget_type(dataset);
        
        if (offset == -1) // numDims = 2
        {
            H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
        }
        else
        {
            hid_t space = H5Dget_space(dataset);
            hsize_t unsigned_offset[3];
            hsize_t count[3];
            int numDims = H5Sget_simple_extent_ndims(space);
            
            if (numDims < 2)
            {
                logged << "I have no idea what this error message should say!" << std::endl;
                sendLog();
            }
            else if (numDims < 3)
            {
                H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
            }
            else
            {
                /* make data space thing */
                hsize_t dims[numDims];
                H5Sget_simple_extent_dims(space, dims, NULL);
                
                count[0] = 1;
                count[1] = dims[1];
                count[2] = dims[2];
                
                unsigned_offset[0] = offset;
                unsigned_offset[1] = 0;
                unsigned_offset[2] = 0;
                
                H5Sselect_hyperslab(space, H5S_SELECT_SET, unsigned_offset, NULL, count, NULL);
                
                /* make memory space thing */
                
                hid_t memspace_id = H5Screate_simple (3, count, NULL);
                
                /* read */
                
                H5Dread(dataset, type, memspace_id, space, H5P_DEFAULT, *buffer);
                
                H5Sclose(space);
                H5Sclose(memspace_id);
            }
        }
        
        H5Dclose(dataset);
        H5Tclose(type);
        
        return true; // success
    }
    catch (std::exception e)
    {
        return false; // failure
    }
}

void Hdf5Manager::turnOffErrors()
{
 //   return;
    
    /* Save old error handler */
//   hid_t error_stack = H5Eget_current_stack();
    
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    
    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    
//    H5Eset_current_stack(error_stack);
}

void Hdf5Manager::turnOnErrors()
{
 //   return;
    
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
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

void Hdf5Manager::identifiersFromAddress(std::map<std::string, int> *map, std::vector<std::string> *list, std::string idAddress)
{
    char *buffer;
    int idSizes = hdf5MallocBytesForDataset(idAddress, (void **)&buffer);
    dataForAddress(idAddress, (void **)&buffer);
    // nah you need to make this memspace thingy
    int lastStart = 0;
    int count = 0;
    
    for (int i = 0; i < idSizes; i++)
    {
        if (buffer[i] == '\0')
        {
            std::string anId = std::string(&buffer[lastStart]);
            lastStart = i + 1;
            
            if (anId.size() > 0)
            {
                list->push_back(anId);
                (*map)[lastComponent(anId)] = count;
                count++;
            }
        }
    }
}
    
int Hdf5Manager::getSubTypeForIndex(std::string address, int objIdx, H5G_obj_t *type)
{
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    try
    {
        hid_t group = H5Gopen1(handle, address.c_str());
        *type = H5Gget_objtype_by_idx(group, objIdx);
        if (group >= 0)
        {
            H5Gclose(group);
        }
        else
        {
            logged << "Could not get sub group types for " << address << std::endl;
            sendLog();
        }

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
        std::lock_guard<std::mutex> lg(readingHdf5);
        hid_t groupAll = H5Gopen1(handle, address.c_str());
        hsize_t objectNum = 0;
        herr_t error = H5Gget_num_objs(groupAll, &objectNum);
        
        if (error < 0)
        {
            logged << "Could not get sub group types for " << address << std::endl;
            sendLog();
            return types;
        }
        
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
 
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    hid_t groupAll = H5Gopen1(handle, parentAddress.c_str());
    
    if (groupAll >= 0)
    {
        hid_t newGroup = H5Gcreate1(groupAll, newGroupName.c_str(), 0);
        
        if (newGroup < 0)
        {
            return false;
        }
        else
        {
            H5Gclose(newGroup);
        }
        
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
    std::lock_guard<std::mutex> lg(readingHdf5);
    
    std::vector<std::string> names;
    
    hid_t groupAll = H5Gopen(handle, address.c_str(), H5P_DEFAULT);
    hsize_t objectNum = 0;
    herr_t error = H5Gget_num_objs(groupAll, &objectNum);
    
    if (error < 0)
    {
        logged << "Could not get sub group names for " << address << std::endl;
        sendLog();

        return names;
    }
    
    for (int i = 0; i < objectNum; i++)
    {
        char objectName[MAX_NAME] = "";
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
    std::lock_guard<std::mutex> lg(readingHdf5);
    
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
    std::vector<std::string> components = FileReader::split(path, '/');
    
    if (components.size() == 0)
    {
        return std::string("");
    }
    
    if (components[components.size() - 1] == "")
    {
        components.pop_back();
    }
    
    return components[components.size() - 1];
}

int Hdf5Manager::pathComponentCount(std::string path)
{
    std::vector<std::string> components = FileReader::split(path, '/');
    
    return (int)components.size();
}

std::string Hdf5Manager::truncateLastComponent(std::string path)
{
    int componentCount = pathComponentCount(path);
    
    std::string truncated = truncatePath(path, componentCount - 1);
    
    return truncated;
}

std::string Hdf5Manager::truncatePath(std::string path, int numToTruncate)
{
    std::vector<std::string> components = FileReader::split(path, '/');
    std::string added = "/";
    
    for (int i = 0; i < numToTruncate && i < components.size(); i++)
    {
        if (components[i].length() > 0)
        {
            added += components[i] + "/";
        }
    }
    
    if (added.length() > 1 && added[added.size() - 1] == '/')
    {
        added.erase(added.size() - 1);

    }
    
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
    
    for (int i = 0; i <= componentCount; i++)
    {
        std::string pathToCheck = truncatePath(address, i);
        
        if (pathToCheck.length() == 0)
        {
            continue;
        }
        
        bool exists = groupExists(pathToCheck);
        
        if (exists)
        {
       //     logged << "Address " << pathToCheck << " exists." << std::endl;
        }
        else
        {
            bool success = createLastComponentAsGroup(pathToCheck);
            
            if (success)
            {
          //      logged << "Address " << pathToCheck << " successfully created." << std::endl;
            }
            else
            {
          //      logged << "Address " << pathToCheck << " creation failed." << std::endl;
            }
        }
    }
 
  //  sendLog();

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
    
    if (handle < 0)
    {
        logged << "Warning: could not open file handle " << newName << std::endl;
        Logger::setShouldExit();
        sendLog();
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