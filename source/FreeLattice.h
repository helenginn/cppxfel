//
//  FreeLattice.h
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


#ifndef __cppxfel__FreeLattice__
#define __cppxfel__FreeLattice__

#include <stdio.h>
#include <vector>
#include "parameters.h"

class FreeLattice
{
protected:
    std::vector<SpotVectorPtr> spotVectors;
    std::vector<SpotVectorPtr> expandedSpotVectors;
    std::vector<double> orderedDistances;

    void calculateExpandedVectors(bool originOnly);
    void startingAngles(double a, double b, double c, double alpha, double beta, double gamma);
public:
    FreeLattice();
    FreeLattice(double a, double b, double c, double alpha, double beta, double gamma);

    void addExpanded();
    void powderPattern(bool originOnly = false, std::string filename = "latticePowder.csv");
    void anglePattern(bool originOnly = false);

    size_t orderedDistanceCount()
    {
        return orderedDistances.size();
    }

    double orderedDistance(int i)
    {
        return orderedDistances[i];
    }
};

#endif /* defined(__cppxfel__FreeLattice__) */
