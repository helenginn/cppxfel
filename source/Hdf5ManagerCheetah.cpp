//
//  Hdf5ManagerCheetah.cpp
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

#include <stdio.h>
#include "Hdf5ManagerCheetah.h"
#include "Hdf5ManagerCheetahLCLS.h"
#include "Hdf5ManagerCheetahSacla.h"
#include "FileParser.h"
#include "misc.h"

std::string Hdf5ManagerCheetah::maskAddress;
std::vector<Hdf5ManagerCheetahPtr> Hdf5ManagerCheetah::cheetahManagers;
std::mutex Hdf5ManagerCheetah::readingPaths;

void Hdf5ManagerCheetah::initialiseCheetahManagers()
{
    if (cheetahManagers.size() > 0)
        return;

    std::vector<std::string> hdf5FileGlobs = FileParser::getKey("HDF5_SOURCE_FILES", std::vector<std::string>());
    std::ostringstream logged;

    for (int i = 0; i < hdf5FileGlobs.size(); i++)
    {
        std::vector<std::string> hdf5Files = glob(hdf5FileGlobs[i]);

        logged << "HDF5_SOURCE_FILES entry (no. " << i << ") matches: " << std::endl;

        for (int j = 0; j < hdf5Files.size(); j++)
        {
            std::string aFilename = hdf5Files[j];
            Hdf5ManagerCheetahPtr cheetahPtr;

            int guessLCLS = (aFilename.find("cxi") != std::string::npos);
            int guessLaser = (guessLCLS) ? 0 : 1;

            int laserInt = FileParser::getKey("FREE_ELECTRON_LASER", guessLaser);
            FreeElectronLaserType laser = (FreeElectronLaserType)laserInt;

            switch (laser) {
                case FreeElectronLaserTypeLCLS:
                    cheetahPtr = Hdf5ManagerCheetahLCLS::makeManager(aFilename);
                    break;
                                case FreeElectronLaserTypeEuropeanXFEL:
                                        cheetahPtr = Hdf5ManagerCheetahLCLS::makeManager(aFilename);
                                        break;
                case FreeElectronLaserTypeSACLA:
                    cheetahPtr = Hdf5ManagerCheetahSacla::makeManager(aFilename);
                    break;
                default:
                    cheetahPtr = Hdf5ManagerCheetahSacla::makeManager(aFilename);
                    break;
            }

            cheetahManagers.push_back(cheetahPtr);

            logged << hdf5Files[j] << ", ";
        }

        logged << std::endl;
    }

    if (cheetahManagers.size())
    {
        logged << "... now managing " << cheetahManagers.size() << " hdf5 image source files." << std::endl;
    }
    if (!cheetahManagers.size())
    {
        if (hdf5FileGlobs.size())
        {
            logged << "You specified some hdf5 image source files but none could be found. Please check existence and location of files." << std::endl;
            LoggableObject::staticLogAndExit(logged);
        }
    }

    Logger::mainLogger->addStream(&logged);
}

Hdf5ManagerCheetahPtr Hdf5ManagerCheetah::hdf5ManagerForImage(std::string imageName)
{
    // make a map!
    for (int i = 0; i < cheetahManagers.size(); i++)
    {
        Hdf5ManagerCheetahPtr cheetahManager = cheetahManagers[i];

        std::string address = cheetahManager->addressForImage(imageName);

        if (address.length() > 0)
        {
            return cheetahManager;
        }
    }

    return Hdf5ManagerCheetahSaclaPtr();
}

void Hdf5ManagerCheetah::closeHdf5Files()
{
    for (int i = 0; i < cheetahManagers.size(); i++)
    {
        cheetahManagers[i]->closeHdf5();
    }
}

std::string Hdf5ManagerCheetah::addressForImage(std::string imageName)
{
    // maybe imageName has .img extension, so let's get rid of it
        unsigned long h5Pos = imageName.find(".h5");
        if (h5Pos != std::string::npos)
        {
                imageName[h5Pos] = '_';
        }
        std::string baseName = getBaseFilename(imageName);

    if (!imagePathMap.count(baseName))
    {
        return "";
    }

    int index = imagePathMap[baseName];
    return imagePaths[index];
}
