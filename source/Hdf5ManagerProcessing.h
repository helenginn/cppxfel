//
//  Hdf5ManagerProcessing.h
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

#ifndef __cppxfel__Hdf5ManagerProcessing__
#define __cppxfel__Hdf5ManagerProcessing__

#include <stdio.h>
#include "FileParser.h"
#include "Hdf5ManagerCheetah.h"

class Hdf5ManagerProcessing  : public Hdf5ManagerCheetah
{
private:
    static Hdf5ManagerProcessingPtr processingManager;
    
public:
    Hdf5ManagerProcessing(std::string filename) : Hdf5ManagerCheetah(filename, Hdf5AccessTypeReadWrite)
    {
        
    }
    
    static void setupProcessingManager();
    
    static Hdf5ManagerProcessingPtr getProcessingManager();
    
};

#endif /* defined(__cppxfel__Hdf5ManagerProcessing__) */
