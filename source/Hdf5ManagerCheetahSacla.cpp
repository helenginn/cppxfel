//
//  Hdf5ManagerCheetahSacla.cpp
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

#include "Hdf5ManagerCheetahSacla.h"
#include "FileParser.h"

Hdf5ManagerCheetahPtr Hdf5ManagerCheetahSacla::makeManager(std::string filename)
{
    Hdf5ManagerCheetahSaclaPtr cheetahPtr = Hdf5ManagerCheetahSaclaPtr(new Hdf5ManagerCheetahSacla(filename));

    return std::static_pointer_cast<Hdf5ManagerCheetah>(cheetahPtr);
}

bool Hdf5ManagerCheetahSacla::getImageSize(std::string address, int *dims)
{
    std::string dataAddress = concatenatePaths(address, "data");
    return Hdf5Manager::getImageSize(dataAddress, dims);
}

bool Hdf5ManagerCheetahSacla::dataForImage(std::string address, void **buffer, bool rawAddress)
{
    std::string dataAddress = concatenatePaths(address, "data");

    if (rawAddress)
    {
    //    dataAddress = address;
    }

        return Hdf5Manager::dataForAddress(dataAddress, buffer);
}

double Hdf5ManagerCheetahSacla::wavelengthForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "photon_wavelength_A");
    return Hdf5Manager::dataForAddress(dataAddress, buffer);
}

int Hdf5ManagerCheetahSacla::hdf5MallocBytesForImage(std::string address, void **buffer)
{
    std::string dataAddress = concatenatePaths(address, "data");
    return Hdf5Manager::hdf5MallocBytesForDataset(dataAddress, buffer);
}

size_t Hdf5ManagerCheetahSacla::bytesPerTypeForImageAddress(std::string address)
{
    std::string dataAddress = concatenatePaths(address, "data");
    return Hdf5Manager::bytesPerTypeForDatasetAddress(dataAddress);
}
