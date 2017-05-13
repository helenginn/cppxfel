//
//  Hdf5ManagerProcessing.cpp
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

#include "Hdf5ManagerProcessing.h"

Hdf5ManagerProcessingPtr Hdf5ManagerProcessing::processingManager = Hdf5ManagerProcessingPtr();


void Hdf5ManagerProcessing::setupProcessingManager()
{
    std::string filename = FileParser::getKey("HDF5_OUTPUT_FILE", std::string(""));
    
    if (filename == "")
    {
        processingManager = Hdf5ManagerProcessingPtr();
        return;
    }
    
    processingManager = Hdf5ManagerProcessingPtr(new Hdf5ManagerProcessing(filename));
}

Hdf5ManagerProcessingPtr Hdf5ManagerProcessing::getProcessingManager()
{
    return processingManager;
}
