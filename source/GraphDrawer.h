/*
 * GraphDrawer.h
 */
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

#ifndef GRAPHDRAWER_H_
#define GRAPHDRAWER_H_

#include "parameters.h"
#include "StatisticsManager.h"
#include "MtzManager.h"
#include <boost/variant.hpp>
#include <map>

typedef std::map<std::string, boost::variant<double, std::string> > GraphMap;

class GraphDrawer
{
private:
	MtzManager *mtz;

public:
	GraphDrawer(MtzManager *mtz);
	virtual ~GraphDrawer();

    void resolutionStatsCSV(std::vector<MtzManager *>& managers);
	
    void partialityPNG(MtzPtr mtz, double maxRes = 0);

    void plotSingleMillerFromMtzs(std::vector<MtzPtr> mtzs, int h, int k, int l);
    void plotReflectionFromMtzs(std::vector<MtzPtr> mtzs, int h = 0, int k = 0, int l = 0);
    void plotOrientationStats(vector<MtzPtr> mtzs);
    void cutoutIntegrationAreas(std::vector<MtzPtr> mtzs, int h = 0, int k = 0, int l = 0);
    
    void partialityPNGResolutionShell(std::string filename, double meanWavelength,
                                      std::vector<ReflectionPtr> refRefls,
                                      std::vector<ReflectionPtr> imageRefls,
                                      double minRes, double maxRes);
    
    MtzManager*& getMtz()
	{
		return mtz;
	}

	void setMtz(MtzManager*& mtz)
	{
		this->mtz = mtz;
	}
};

#endif /* GRAPHDRAWER_H_ */
