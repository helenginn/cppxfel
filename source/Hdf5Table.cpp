//
//  Hdf5Table.cpp
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

#include "Hdf5Table.h"
#include "Hdf5ManagerProcessing.h"

bool Hdf5Table::writeToManager(Hdf5ManagerProcessingPtr manager, std::string address)
{
    return manager->writeTable(*this, address);
}

int Hdf5Table::readSizeFromManager(Hdf5ManagerProcessingPtr manager, std::string address)
{
    return manager->readSizeForTable(*this, address);
}

int Hdf5Table::readFromManager(Hdf5ManagerProcessingPtr manager, std::string address, void *data)
{
    return manager->recordsForTable(*this, address, data);
}

Hdf5Table::~Hdf5Table()
{
    if (fieldOffsets != NULL)
        free(fieldOffsets);
    fieldOffsets = NULL;

    if (fieldSizes != NULL)
        free(fieldSizes);
    fieldSizes = NULL;

    if (fieldTypes != NULL)
        free(fieldTypes);
    fieldTypes = NULL;

}
