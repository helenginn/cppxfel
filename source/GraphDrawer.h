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

	static std::string generateFilename(std::string stem);
	static std::string generateFilename(std::string stem, std::string ext);


	std::string plot(std::string filename, GraphMap properties,
			vector<double> x, vector<double> y,
			vector<double> x2, vector<double> y2);

	std::string plot(std::string filename, GraphMap properties,
			vector<vector<double> > x, vector<vector<double> > y,
			vector<double> x2, vector<double> y2);

	std::string plot(std::string filename,
			GraphMap properties,
			vector<double> x, vector<double> y);

	std::string plot(std::string filename, GraphMap properties,
			vector<vector<double> > xs, vector<vector<double> > ys);

    void resolutionStatsCSV(std::vector<MtzManager *>& managers);
	
    void plotPolarisation(vector<MtzPtr> mtzs);
	void correlationPlot(std::string filename, double xMax = 0, double yMax = 0);
	void partialityPNG(MtzPtr mtz, double maxRes = 0);

    void plotSingleMillerFromMtzs(std::vector<MtzPtr> mtzs, int h, int k, int l);
    void plotReflectionFromMtzs(std::vector<MtzPtr> mtzs, int h = 0, int k = 0, int l = 0);
    void plotOrientationStats(vector<MtzPtr> mtzs);
    void plotPartialityStats(int h = 0, int k = 0, int l = 0);
    void cutoutIntegrationAreas(std::vector<MtzPtr> mtzs, int h = 0, int k = 0, int l = 0);
    
    void partialityPNGResolutionShell(std::string filename,
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
