/*
 * GraphDrawer.h
 *
 *  Created on: 22 Oct 2014
 *      Author: helenginn
 */

#ifndef GRAPHDRAWER_H_
#define GRAPHDRAWER_H_

#include "StatisticsManager.h"
#include "MtzManager.h"
#include <boost/variant.hpp>
#include <map>

typedef std::map<string, boost::variant<double, string> > GraphMap;

class GraphDrawer
{
private:
	MtzManager *mtz;
//	vector<MtzManager *>mtzs;

public:
	GraphDrawer(MtzManager *mtz);
	virtual ~GraphDrawer();

	static std::string generateFilename(string stem);
	static std::string generateFilename(std::string stem, std::string ext);

#ifdef MAC
	std::string plot(string filename, GraphMap properties,
			std::vector<double> x, std::vector<double> y,
			std::vector<double> x2, std::vector<double> y2);

	std::string plot(string filename, GraphMap properties,
			std::vector<std::vector<double> > x, std::vector<std::vector<double> > y,
			std::vector<double> x2, std::vector<double> y2);

	std::string plot(string filename,
			GraphMap properties,
			std::vector<double> x, std::vector<double> y);

	std::string plot(string filename, GraphMap properties,
			std::vector<std::vector<double> > xs, std::vector<std::vector<double> > ys);

	void resolutionStatsPlot(vector<MtzManager *>& managers, string filename = "resolution_stats",
			GraphMap properties = GraphMap(), bool intensityBins = false, bool image = false);

    void plotPartialityStats();
    void plotPolarisation(std::vector<MtzPtr> mtzs);
	void correlationPlot(string filename, double xMax = 0, double yMax = 0);
	void partialityPlot(string filename, GraphMap properties = GraphMap());
	void bFactorPlot(vector<MtzManager *>& managers,
			string filename = "all_gradients", GraphMap properties = GraphMap());
#endif

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
