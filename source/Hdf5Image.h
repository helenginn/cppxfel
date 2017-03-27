//
//  Hdf5Image.h
//  cppxfel
//
//  Created by Helen Ginn on 15/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Hdf5Image__
#define __cppxfel__Hdf5Image__

#include <stdio.h>
#include "Image.h"
#include "Hdf5ManagerCheetahSacla.h"
#include "Hdf5Table.h"

class Hdf5Image : public Image
{
private:
    void failureMessage();
    virtual void loadImage();
    std::string imageAddress;
    std::string findAddress();
    Hdf5ManagerCheetahPtr chManager;
    Hdf5Table spotTable;
    void createSpotTable();
    void getWavelengthFromHdf5();
protected:
    
public:
    Hdf5Image(std::string filename = "", double wavelength = 0,
              double distance = 0) : Image(filename, wavelength, distance)
    {
        _isMask = false;
        imageAddress = std::string();
        createSpotTable();
        getWavelengthFromHdf5();

        if (!imageAddress.length())
        {
            return;
        }
    
        getWavelengthFromHdf5();
    };
    
    void setMask()
    {
        _isMask = true;
        imageAddress = Hdf5ManagerCheetah::getMaskAddress();
        setFilename("mask.img");
        loadImage();

        if (shortData.size() || data.size())
        {
            setImageMask(shared_from_this());
        }
        else
        {
            std::ostringstream logged;
            logged << "Cannot find mask from CHEETAH_MASK_ADDRESS." << std::endl;
            Logger::log(logged);
        }
    }
    
    Hdf5ManagerCheetahPtr getManager();
    
    void writeSpotsList(std::string spotFile);
    void processSpotList();
    void loadCrystals();
	virtual int getFiducial();

    virtual ImageClass getClass()
    {
        return ImageClassHdf5;
    }

    std::string getAddress()
    {
        return findAddress();
    }
    
    virtual ~Hdf5Image() {};
};

#endif /* defined(__cppxfel__Hdf5Image__) */
