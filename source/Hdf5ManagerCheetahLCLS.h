//
//  Hdf5ManagerCheetahLCLS.h
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

#ifndef __cppxfel__Hdf5ManagerCheetahLCLS__
#define __cppxfel__Hdf5ManagerCheetahLCLS__

#include <stdio.h>
#include "Hdf5ManagerCheetah.h"
#include "FileParser.h"

class Hdf5ManagerCheetahLCLS : public Hdf5ManagerCheetah
{
private:
    // place where all the stuff is
    // stupid comments, sorry helen.
    // I'm tired.
    std::string idAddress;
    std::string dataAddress;
    std::string wavelengthAddress;
    std::vector<double> wavelengths;


public:
    Hdf5ManagerCheetahLCLS(std::string newName) : Hdf5ManagerCheetah(newName)
    {
        idAddress = FileParser::getKey("CHEETAH_ID_ADDRESSES",
                                                   std::string("entry_1/data_1/experiment_identifier"));
        dataAddress = FileParser::getKey("CHEETAH_DATA_ADDRESSES", std::string("entry_1/data_1/data"));
        wavelengthAddress = FileParser::getKey("CHEETAH_WAVELENGTH_ADDRESSES", std::string("LCLS/photon_wavelength_A"));
        identifiersFromAddress(&imagePathMap, &imagePaths, idAddress);
        prepareWavelengths();
    }

    static Hdf5ManagerCheetahPtr makeManager(std::string filename);

    void prepareWavelengths();
    virtual double wavelengthForImage(std::string address, void **buffer);
    virtual bool dataForImage(std::string address, void **buffer, bool rawAddress = false);

    virtual int hdf5MallocBytesForImage(std::string address, void **buffer);
    virtual size_t bytesPerTypeForImageAddress(std::string address);
    virtual bool getImageSize(std::string address, int *dims);
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahLCLS__) */
