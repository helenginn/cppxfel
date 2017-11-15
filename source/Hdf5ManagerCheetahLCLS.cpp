//
//  Hdf5ManagerCheetahLCLS.cpp
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

#include "Hdf5ManagerCheetahLCLS.h"
#include <iterator>
#include <algorithm>

Hdf5ManagerCheetahPtr Hdf5ManagerCheetahLCLS::makeManager(std::string filename)
{
    Hdf5ManagerCheetahLCLSPtr cheetahPtr = Hdf5ManagerCheetahLCLSPtr(new Hdf5ManagerCheetahLCLS(filename));

    return std::static_pointer_cast<Hdf5ManagerCheetah>(cheetahPtr);
}

bool Hdf5ManagerCheetahLCLS::getImageSize(std::string address, int *dims)
{
    return Hdf5Manager::getImageSize(dataAddress, dims);
}

bool Hdf5ManagerCheetahLCLS::dataForImage(std::string address, void **buffer, bool rawAddress)
{
    // need to find the index of the correct image

    // lazy as.

    if (rawAddress)
    {
        return Hdf5Manager::dataForAddress(address, buffer, true);
    }

    std::vector<std::string>::iterator it = std::find(imagePaths.begin(),
                                                imagePaths.end(),
                                                address);

    if (it != imagePaths.end() && *it == address)
    {
        // am a terrible caster. sorry.
        int index = (int)(it - imagePaths.begin());

        return Hdf5Manager::dataForAddress(dataAddress, buffer, index);
    }

    return false;
}

void Hdf5ManagerCheetahLCLS::prepareWavelengths()
{
    double *tempWaves;
    size_t bytes = Hdf5Manager::hdf5MallocBytesForDataset(wavelengthAddress, (void **)&tempWaves);
    size_t waveNum = bytes / sizeof(double);
    Hdf5Manager::dataForAddress(wavelengthAddress, (void **)&tempWaves);
    wavelengths.resize(waveNum);
    memcpy(&wavelengths[0], &tempWaves[0], bytes);
}

double Hdf5ManagerCheetahLCLS::wavelengthForImage(std::string address, void **buffer)
{
    std::vector<std::string>::iterator it = std::find(imagePaths.begin(),
                                                      imagePaths.end(),
                                                      address);

    if (it != imagePaths.end() && *it == address)
    {
        size_t index = it - imagePaths.begin();
        memcpy(*buffer, &wavelengths[index], sizeof(double));
        return wavelengths[index];
    }

    return 0;
}

int Hdf5ManagerCheetahLCLS::hdf5MallocBytesForImage(std::string address, void **buffer)
{
    return Hdf5Manager::hdf5MallocBytesForDataset(dataAddress, buffer);
}

size_t Hdf5ManagerCheetahLCLS::bytesPerTypeForImageAddress(std::string address)
{
    return Hdf5Manager::bytesPerTypeForDatasetAddress(dataAddress);
}
