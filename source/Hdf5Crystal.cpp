//
//  Hdf5Crystal.cpp
//  cppxfel
//
//  Created by Helen Ginn on 24/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Hdf5Crystal.h"
#include "Hdf5Image.h"
#include "Hdf5ManagerProcessing.h"
#include "Miller.h"
#include "misc.h"
#include "Matrix.h"

#define HDF5MILLER_FIELD_COUNT 18

typedef struct
{
    int h;
    int k;
    int l;
    bool free;
    double rawIntensity;
    double sigma; // error used during refinement
    double countingSigma;
    double partiality;
    double wavelength;
    float resolution;
    float phase;
    float bFactor; // make scale part of crystal metadata
    float partialCutoff;
    float lastX; // lab coordinates, no shifts
    float lastY; // lab coordinates, no shifts
    float shiftX; // shift from panel-corrected lab coordinates to highest pixel
    float shiftY; // shift from panel-corrected lab coordinates to highest pixel
    RejectReason rejectReason;
} Hdf5Miller;

void Hdf5Crystal::writeReflectionData(std::string address)
{
    Hdf5ManagerProcessingPtr processingManager = Hdf5ManagerProcessing::getProcessingManager();
    
    if (!processingManager)
    {
        return;
    }
    
    int count = 0;
    int millerCount = this->millerCount();
    Hdf5Miller *millerData = (Hdf5Miller *)malloc(sizeof(Hdf5Miller) * millerCount);
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        for (int j = 0; j < reflection(i)->millerCount(); j++)
        {
            MillerPtr miller = reflection(i)->miller(j);
            
            millerData[i].h = miller->getH();
            millerData[i].k = miller->getK();
            millerData[i].l = miller->getL();
            millerData[i].free = miller->isFree();
            millerData[i].rawIntensity = miller->getRawestIntensity();
            millerData[i].sigma = miller->getSigma();
            millerData[i].countingSigma = miller->getRawCountingSigma();
            millerData[i].partiality = miller->getPartiality();
            millerData[i].wavelength = miller->getWavelength();
            millerData[i].resolution = miller->resolution();
            millerData[i].phase = miller->getPhase();
            millerData[i].bFactor = miller->getBFactor();
            millerData[i].partialCutoff = miller->getPartialCutoff();
            millerData[i].lastX = miller->getLastX();
            millerData[i].lastY = miller->getLastY();
            millerData[i].shiftX = miller->getShift().first;
            millerData[i].shiftY = miller->getShift().second;
            millerData[i].rejectReason = miller->getRejectedReason();
        }
        
        count++;
    }
    
    millerTable.setData(millerData);
    millerTable.setNumberOfRecords(millerCount);
    
    bool success = millerTable.writeToManager(processingManager, address);
    
    
    logged << "Writing reflection data to HDF5 file was " << (success ? "successful." : "a failure.") << std::endl;
    
    free(millerData);
    
    sendLog();
}

void Hdf5Crystal::writeCrystalData(std::string address)
{
    Hdf5ManagerProcessingPtr manager = Hdf5ManagerProcessing::getProcessingManager();
    
    // space group
    
    std::string spaceGroupAddress = Hdf5Manager::concatenatePaths(address, "spaceGroup");
    manager->writeSingleValueDataset(spaceGroupAddress, low_group->spg_num, H5T_NATIVE_INT);
    
    // rotation matrix
    
    std::string unitCellMatrixAddress = Hdf5Manager::concatenatePaths(address, "unitcell_matrix");
    std::string rotationAddress = Hdf5Manager::concatenatePaths(address, "rotation_matrix");
    
    this->getMatrix()->getUnitCell()->writeToHdf5(unitCellMatrixAddress);
    this->getMatrix()->getRotation()->writeToHdf5(rotationAddress);
    
    // unit cell data
    
    double *unitCell = new double[6];
    memcpy(&unitCell[0], &cellDim[0], sizeof(double) * 3);
    memcpy(&unitCell[3], &cellAngles[0], sizeof(double) * 3);
    
    std::string unitCellParamsAddress = Hdf5Manager::concatenatePaths(address, "unitcell_params");
    
    hsize_t dim = 6;
    manager->createDataset(unitCellParamsAddress, 1, &dim, H5T_NATIVE_DOUBLE);
    manager->writeDataset(unitCellParamsAddress, (void **)&unitCell, H5T_NATIVE_DOUBLE);
    
    delete [] unitCell;
    
    // general scaling/indexing ambiguity
    
    std::string scaleAddress = Hdf5Manager::concatenatePaths(address, "scale");
    std::string ambiguityAddress = Hdf5Manager::concatenatePaths(address, "ambiguity");
    
    manager->writeSingleValueDataset(scaleAddress, scale, H5T_NATIVE_DOUBLE);
    manager->writeSingleValueDataset(ambiguityAddress, activeAmbiguity, H5T_NATIVE_INT32);
   
    // partiality-specific parameters
    
    std::string hRotAddress = Hdf5Manager::concatenatePaths(address, "hRot");
    std::string kRotAddress = Hdf5Manager::concatenatePaths(address, "kRot");
    std::string mosaicityAddress = Hdf5Manager::concatenatePaths(address, "mosaicity");
    std::string rlpSizeAddress = Hdf5Manager::concatenatePaths(address, "rlpSize");
    
    manager->writeSingleValueDataset(hRotAddress, hRot, H5T_NATIVE_DOUBLE);
    manager->writeSingleValueDataset(kRotAddress, kRot, H5T_NATIVE_DOUBLE);
    manager->writeSingleValueDataset(mosaicityAddress, mosaicity, H5T_NATIVE_DOUBLE);
    manager->writeSingleValueDataset(rlpSizeAddress, spotSize, H5T_NATIVE_DOUBLE);
    
    // these will need turning into a Beam object in time.
    
    std::string wavelengthAddress = Hdf5Manager::concatenatePaths(address, "wavelength");
    std::string bandwidthAddress = Hdf5Manager::concatenatePaths(address, "bandwidth");
    std::string exponentAddress = Hdf5Manager::concatenatePaths(address, "exponent");
    
    manager->writeSingleValueDataset(wavelengthAddress, wavelength, H5T_NATIVE_DOUBLE);
    manager->writeSingleValueDataset(bandwidthAddress, bandwidth, H5T_NATIVE_DOUBLE);
    manager->writeSingleValueDataset(exponentAddress, exponent, H5T_NATIVE_DOUBLE);
}

// This has a lot of extra commands for back-compatibility to storing MTZs
// but no attention will be paid to the values of these commands
void Hdf5Crystal::writeToFile(std::string newFilename, bool announce, bool shifts, bool includeAmbiguity, bool useCountingSigma)
{
    if (getImagePtr()->getClass() != ImageClassHdf5)
    {
        return;
    }
    
    ImagePtr superParent = getImagePtr();
    
    Hdf5ImagePtr parent = boost::static_pointer_cast<Hdf5Image>(superParent);
    std::string address = parent->getAddress();
    std::string crystalAddress = Hdf5Manager::concatenatePaths(address, newFilename);
    
    Hdf5ManagerProcessingPtr manager = Hdf5ManagerProcessing::getProcessingManager();
    manager->createGroupsFromAddress(crystalAddress);
    
    writeReflectionData(crystalAddress);
    writeCrystalData(crystalAddress);
    
    MtzManager::writeToFile(filename, announce, shifts, includeAmbiguity, useCountingSigma);
}

void Hdf5Crystal::createMillerTable()
{
    size_t *fieldOffsets = (size_t *)malloc(sizeof(size_t) * HDF5MILLER_FIELD_COUNT);
    fieldOffsets[0] = HOFFSET(Hdf5Miller, h);
    fieldOffsets[1] = HOFFSET(Hdf5Miller, k);
    fieldOffsets[2] = HOFFSET(Hdf5Miller, l);
    fieldOffsets[3] = HOFFSET(Hdf5Miller, free);
    fieldOffsets[4] = HOFFSET(Hdf5Miller, rawIntensity);
    fieldOffsets[5] = HOFFSET(Hdf5Miller, sigma);
    fieldOffsets[6] = HOFFSET(Hdf5Miller, countingSigma);
    fieldOffsets[7] = HOFFSET(Hdf5Miller, partiality);
    fieldOffsets[8] = HOFFSET(Hdf5Miller, wavelength);
    fieldOffsets[9] = HOFFSET(Hdf5Miller, resolution);
    fieldOffsets[10] = HOFFSET(Hdf5Miller, phase);
    fieldOffsets[11] = HOFFSET(Hdf5Miller, bFactor);
    fieldOffsets[12] = HOFFSET(Hdf5Miller, partialCutoff);
    fieldOffsets[13] = HOFFSET(Hdf5Miller, lastX);
    fieldOffsets[14] = HOFFSET(Hdf5Miller, lastY);
    fieldOffsets[15] = HOFFSET(Hdf5Miller, shiftX);
    fieldOffsets[16] = HOFFSET(Hdf5Miller, shiftY);
    fieldOffsets[17] = HOFFSET(Hdf5Miller, rejectReason);
    
    
    hid_t *fieldTypes = (hid_t *)malloc(sizeof(hid_t) * HDF5MILLER_FIELD_COUNT);
    fieldTypes[0] = H5T_NATIVE_INT32;
    fieldTypes[1] = H5T_NATIVE_INT32;
    fieldTypes[2] = H5T_NATIVE_INT32;
    fieldTypes[3] = H5T_NATIVE_CHAR;
    fieldTypes[4] = H5T_NATIVE_DOUBLE;
    fieldTypes[5] = H5T_NATIVE_DOUBLE;
    fieldTypes[6] = H5T_NATIVE_DOUBLE;
    fieldTypes[7] = H5T_NATIVE_DOUBLE;
    fieldTypes[8] = H5T_NATIVE_DOUBLE;
    fieldTypes[9] = H5T_NATIVE_FLOAT;
    fieldTypes[10] = H5T_NATIVE_FLOAT;
    fieldTypes[11] = H5T_NATIVE_FLOAT;
    fieldTypes[12] = H5T_NATIVE_FLOAT;
    fieldTypes[13] = H5T_NATIVE_FLOAT;
    fieldTypes[14] = H5T_NATIVE_FLOAT;
    fieldTypes[15] = H5T_NATIVE_FLOAT;
    fieldTypes[16] = H5T_NATIVE_FLOAT;
    fieldTypes[17] = H5T_NATIVE_INT32;
    
    size_t *fieldSizes = (size_t *)malloc(sizeof(size_t) * HDF5MILLER_FIELD_COUNT);
    fieldSizes[0] = member_size(Hdf5Miller, h);
    fieldSizes[1] = member_size(Hdf5Miller, k);
    fieldSizes[2] = member_size(Hdf5Miller, l);
    fieldSizes[3] = member_size(Hdf5Miller, free);
    fieldSizes[4] = member_size(Hdf5Miller, rawIntensity);
    fieldSizes[5] = member_size(Hdf5Miller, sigma);
    fieldSizes[6] = member_size(Hdf5Miller, countingSigma);
    fieldSizes[7] = member_size(Hdf5Miller, partiality);
    fieldSizes[8] = member_size(Hdf5Miller, wavelength);
    fieldSizes[9] = member_size(Hdf5Miller, resolution);
    fieldSizes[10] = member_size(Hdf5Miller, phase);
    fieldSizes[11] = member_size(Hdf5Miller, bFactor);
    fieldSizes[12] = member_size(Hdf5Miller, partialCutoff);
    fieldSizes[13] = member_size(Hdf5Miller, lastX);
    fieldSizes[14] = member_size(Hdf5Miller, lastY);
    fieldSizes[15] = member_size(Hdf5Miller, shiftX);
    fieldSizes[16] = member_size(Hdf5Miller, shiftY);
    fieldSizes[17] = member_size(Hdf5Miller, rejectReason);
    
    char **fieldNames = (char **)malloc(sizeof(char **) * HDF5MILLER_FIELD_COUNT);
    
    const char *toCopy[HDF5MILLER_FIELD_COUNT] =
    {
        "h",
        "k",
        "l",
        "free",
        "rawIntensity",
        "sigma",
        "countingSigma",
        "partiality",
        "wavelength",
        "resolution",
        "phase",
        "bFactor",
        "partialCutoff",
        "lastX",
        "lastY",
        "shiftX",
        "shiftY",
        "rejectReason",
    };
    
    for (int i = 0; i < HDF5MILLER_FIELD_COUNT; i++)
    {
        fieldNames[i] = (char *)malloc(sizeof(char*) * strlen(toCopy[i]));
        strcpy(fieldNames[i], toCopy[i]);
    }
    
    millerTable.setTableTitle("refls");
    millerTable.setTableName("refls");
    millerTable.setHeaders((const char**)fieldNames);
    millerTable.setRecordSize(sizeof(Hdf5Miller));
    millerTable.setOffsets(fieldOffsets);
    millerTable.setNumberOfFields(HDF5MILLER_FIELD_COUNT);
    millerTable.setFieldSizes(fieldSizes);
    millerTable.setTypes(fieldTypes);
    millerTable.setCompress(false);
}

void Hdf5Crystal::loadReflections(PartialityModel model, bool special)
{
    MtzManager::loadReflections(model, special);
    
    if (address.length() == 0)
    {
        std::cerr
        << "Cannot load reflections as no address has been specified."
        << std::endl;
        return;
    }
    
    // get space group
    
    std::string scaleAddress = Hdf5Manager::concatenatePaths(address, "scale");
    
    Hdf5ManagerProcessingPtr manager = Hdf5ManagerProcessing::getProcessingManager();
    
    scale = 1;
    double newScale = 1;
    manager->readDatasetValue(scaleAddress, &newScale);
    applyScaleFactor(newScale);

    int spgnum = 0;
    manager->readDatasetValue(scaleAddress, &spgnum);
    setSpaceGroup(spgnum);
    
    if (low_group == NULL)
    {
        setRejected(true);
        return;
    }
    
    double *unitCell = new double[6];
    std::string unitCellParamsAddress = Hdf5Manager::concatenatePaths(address, "unitcell_params");
    
    manager->readDatasetValue(unitCellParamsAddress, &unitCell);
    
    memcpy(&cellDim[0], &unitCell[0], sizeof(double) * 3);
    memcpy(&cellAngles[0], &unitCell[3], sizeof(double) * 3);
    
    delete [] unitCell;
    
    double *unitCellMatrix = new double[16];
    double *rotationMatrix = new double[16];
    
    manager->readDatasetValue(unitCellParamsAddress, &unitCellMatrix);
    manager->readDatasetValue(unitCellParamsAddress, &rotationMatrix);
    
    MatrixPtr unitMat = MatrixPtr(new Matrix());
    unitMat->setComponents(unitCellMatrix);

    MatrixPtr rotMat = MatrixPtr(new Matrix());
    rotMat->setComponents(unitCellMatrix);
    
    delete [] unitCellMatrix;
    delete [] rotationMatrix;

    MatrixPtr newMat = MatrixPtr(new Matrix());
    newMat->setComplexMatrix(unitMat, rotMat);
    
    setMatrix(newMat);
    
    std::string ambiguityAddress = Hdf5Manager::concatenatePaths(address, "ambiguity");
    manager->readDatasetValue(scaleAddress, &activeAmbiguity);
    
    std::string hRotAddress = Hdf5Manager::concatenatePaths(address, "hRot");
    manager->readDatasetValue(hRotAddress, &hRot);
    
    std::string kRotAddress = Hdf5Manager::concatenatePaths(address, "kRot");
    manager->readDatasetValue(kRotAddress, &kRot);
    
    std::string mosaicityAddress = Hdf5Manager::concatenatePaths(address, "mosaicity");
    manager->readDatasetValue(mosaicityAddress, &mosaicity);
    
    std::string rlpSizeAddress = Hdf5Manager::concatenatePaths(address, "rlpSize");
    manager->readDatasetValue(rlpSizeAddress, &spotSize);
    
    std::string wavelengthAddress = Hdf5Manager::concatenatePaths(address, "wavelength");
    manager->readDatasetValue(wavelengthAddress, &wavelength);
    
    std::string bandwidthAddress = Hdf5Manager::concatenatePaths(address, "bandwidth");
    manager->readDatasetValue(bandwidthAddress, &bandwidth);
    
    std::string exponentAddress = Hdf5Manager::concatenatePaths(address, "exponent");
    manager->readDatasetValue(exponentAddress, &exponent);
    
    
    /*
     if (mtz == NULL)
        return;
    
    float *refldata = (float *) CCP4::ccp4_utils_malloc(
                                                        (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    memset(refldata, '\0',
           (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    
    float *adata = (float *) CCP4::ccp4_utils_malloc(
                                                     (mtz->ncol_read) * sizeof(float));
    
    MtzRrefl(mtz->filein, mtz->ncol_read * mtz->nref_filein, refldata);
    
    MTZCOL *col_f = MtzColLookup(mtz, "I");
    
    if (col_f == NULL)
    {
        col_f = MtzColLookup(mtz, "IMEAN");
    }
    
    if (col_f == NULL)
    {
        col_f = MtzColLookup(mtz, "FC");
        
        if (col_f == NULL)
        {
            col_f = MtzColLookup(mtz, "FP");
            
            if (col_f == NULL)
                col_f = MtzColLookup(mtz, "F");
            
            if (col_f != NULL)
            {
                std::cout << "Warning: using observed amplitude instead of intensity" << std::endl;
            }
            else
            {
                std::cout << "Warning: could not find any amplitude columns" << std::endl;
            }
        }
        else
        {
            std::cout << "Warning: using calculated amplitude instead of intensity" << std::endl;
        }
    }
    
    if (col_f == NULL)
    {
        std::cout << "Warning: intensity column not labelled I or IMEAN or amplitude" << std::endl;
        exit(1);
    }
    
    MTZCOL *col_sigf = MtzColLookup(mtz, "SIGI");
    
    if (col_sigf == NULL)
    {
        col_sigf = MtzColLookup(mtz, "SIGIMEAN");
    }
    
    if (col_sigf == NULL)
    {
        col_sigf = MtzColLookup(mtz, "SIGF");
        
        if (col_sigf != NULL)
        {
            std::cout << "Warning: using SIGF column" << std::endl;
        }
    }
    
    if (col_sigf == NULL)
    {
        std::cout << "Warning: sigma column not labelled SIGI or SIGIMEAN" << std::endl;
    }
    
    MTZCOL *col_wave = MtzColLookup(mtz, "WAVE");
    MTZCOL *col_partials = MtzColLookup(mtz, "PART");
    MTZCOL *col_phase = MtzColLookup(mtz, "PHIC");
    MTZCOL *col_shiftx = MtzColLookup(mtz, "SHIFTX");
    MTZCOL *col_shifty = MtzColLookup(mtz, "SHIFTY");
    
    int multiplier = MULTIPLIER;
    
    MTZXTAL **xtals = MtzXtals(mtz);
    float cell[6];
    double coefhkl[6];
    ccp4_lrcell(xtals[0], cell);
    MtzHklcoeffs(cell, coefhkl);
    
    setUnitCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    
    
    int num = 0;
    
    for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
    {
        memcpy(adata, &refldata[i], mtz->ncol_read * sizeof(float));
        
        int h;
        int k;
        int l;
        int offset;
        
        hkls_for_reflection(mtz, adata, &h, &k, &l, &multiplier, &offset);
        
        int reflection_reflection_index = index_for_reflection(h, k, l, false);
        
        float intensity = adata[col_f->source - 1];
        
        float sigma = 1;
        if (col_sigf != NULL)
            sigma = adata[col_sigf->source - 1];
        float wavelength = 1;
        float partiality = 1;
        float phase = 0;
        float shiftX = 0;
        float shiftY = 0;
        
        if (col_wave != NULL && col_partials != NULL)
        {
            wavelength = adata[col_wave->source - 1];
            partiality = adata[col_partials->source - 1];
        }
        
        if (col_shiftx != NULL && col_shifty != NULL)
        {
            shiftX = adata[col_shiftx->source - 1];
            shiftY = adata[col_shifty->source - 1];
        }
        
        if (col_phase != NULL)
        {
            phase = adata[col_phase->source - 1];
        }
        
        if (col_partials == NULL && intensity != intensity)
            continue;
        
        MillerPtr miller = MillerPtr(new Miller(this, h, k, l));
        miller->setData(intensity, sigma, partiality, wavelength);
        miller->setCountingSigma(sigma);
        miller->setFilename(filename);
        //    miller->setPartialityModel(model);
        miller->setPhase(phase);
        miller->setLastX(shiftX);
        miller->setLastY(shiftY);
        //        miller->setShift(std::make_pair(shiftX, shiftY));
        miller->matrix = this->matrix;
        
        Reflection *prevReflection;
        
        this->findReflectionWithId(reflection_reflection_index, &prevReflection);
        
        if (prevReflection != NULL)
        {
            // Exclude unobserved reflections by checking for nan
            if (adata[col_f->source - 1] == adata[col_f->source - 1])
            {
                prevReflection->addMiller(miller); // TODO
                miller->setParent(prevReflection);
            }
            
            // reflection is a repeat so set flag.
        }
        
        if (prevReflection == NULL)
        {
            Reflection *newReflection = new Reflection();
            reflections.push_back(newReflection);
            newReflection->setUnitCell(cell);
            newReflection->setSpaceGroup(low_group->spg_num);
            
            reflections[num]->addMiller(miller);
            miller->setParent(reflections[num]);
            reflections[num]->calculateResolution(this);
            
            this->sortLastReflection();
            
            num++;
        }
    }
    
    free(refldata);
    free(adata);
    
    //	insertionSortReflections();
    
    std::ostringstream log;
    
    log << "Loaded " << mtz->nref_filein << " reflections (" << accepted()
    << " accepted)." << std::endl;
    
    Logger::mainLogger->addStream(&log);
    
    MtzFree(mtz);*/
}