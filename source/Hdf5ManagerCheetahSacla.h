//
//  Hdf5ManagerCheetahSacla.h
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

#ifndef __cppxfel__Hdf5ManagerCheetahSacla__
#define __cppxfel__Hdf5ManagerCheetahSacla__

#include <stdio.h>
#include "Hdf5ManagerCheetah.h"
#include <iostream>
#include <mutex>

class Hdf5ManagerCheetahSacla : public Hdf5ManagerCheetah
{
private:    
public:
    static Hdf5ManagerCheetahPtr makeManager(std::string filename);
    
    virtual bool dataForImage(std::string address, void **buffer, bool rawAddress = false);
    virtual ~Hdf5ManagerCheetahSacla() {};
    virtual double wavelengthForImage(std::string address, void **buffer);
    virtual int hdf5MallocBytesForImage(std::string address, void **buffer);
    virtual bool getImageSize(std::string address, int *dims);
   
    size_t bytesPerTypeForImageAddress(std::string address);
    
    Hdf5ManagerCheetahSacla(std::string newName) : Hdf5ManagerCheetah(newName)
    {
        groupsWithPrefix(&imagePaths, "tag");

		for (int i = 0; i < imagePaths.size(); i++)
        {
            std::string last = lastComponent(imagePaths[i]);
            imagePathMap[last] = i;
        }
    }
    
};

#endif /* defined(__cppxfel__Hdf5ManagerCheetahSacla__) */
