/*
 * GraphDrawer.h
 *
 *  Created on: 22 Oct 2014
 *      Author: helenginn
 */

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
