/*
 * Panel.h
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#ifndef PANEL_H_
#define PANEL_H_

#include "parameters.h"
#include <sstream>
#include <iostream>
#include <vector>
#include "LoggableObject.h"
#include <boost/thread/thread.hpp>

typedef enum
{
    PlotTypeAbsolute,
    PlotTypeRelative
} PlotType;

class Panel : public LoggableObject
{
private:
	Coord topLeft;
	Coord bottomRight;
	Coord bestShift;
    Coord originalShift;
	Coord tilt;
    double tiltWindowSize;
    double averageIntensity;
    bool tiltHorizontalAxis;
	bool beingWrittenTo;
	int allowedSearchSpace;
    PanelTag tag;
    int panelNum;
    double swivel;
    double sinSwivel;
    double cosSwivel;
    double gainScale;
    double height();
    double width();
    static double distanceMultiplier;
    static Coord beamCentre;
    Coord midPoint();

    Coord getSwivelShift(Coord millerCoord, bool isSpot = false);
    Coord getTiltShift(Coord millerCoord);
    Coord millerToSpotCoordShift(Coord millerCoord);
    Coord millerToSpotCoord(Coord millerCoord);
    
    void fractionalCoordinates(Coord coord, Coord *frac);
    void fractionalCoordinates(Miller *miller, Coord *frac);

	static bool usePanelInfo;
    static vector<PanelPtr> panels;
    static vector<PanelPtr> badPanels;
    bool addMiller(MillerPtr miller);

    void centreWindowShift();
    void findAllParameters();
    void findShift(double windowSize, double step, double x = 0, double y = 0);

    std::map<boost::thread::id, vector<MillerPtr> > tempMillers;
    
    static void calculateMetrologyThread(int offset);
    static std::string printAllThreaded();
    double detectorGain(double *error);

    Coord relativeToMidPointForMiller(Coord coord, bool isSpot = false);
    double angleForMiller(Miller *miller);
    double distanceFromMidPointForMiller(Miller *miller);
    void refreshMillerPositions();

    void score();
    static double globalScoreWrapper();
    vector<MillerPtr> millers;
	int defaultShift;
    
    void setPanelNum(int new_value)
    {
        panelNum = new_value;
    }

public:
	Panel(vector<double> dimensions, PanelTag tag = PanelTagNormal);
    Panel(double x1, double y1, double x2, double y2, PanelTag newTag);
    void init(vector<double> dimensions, PanelTag newTag);
    virtual ~Panel();

    bool isMillerCoordInPanel(Coord coord, Coord *topLeft = NULL, Coord *bottomRight = NULL);
    bool isMillerInPanel(Miller *miller);
	static void addMillerToPanelArray(MillerPtr miller);
	static PanelPtr panelForMiller(Miller *miller);
    static PanelPtr panelForSpot(Spot *spot);
    static PanelPtr panelForSpotCoord(Coord coord, PanelPtr *anyBadPanel = NULL);
    static PanelPtr panelForCoord(Coord coord);
    static PanelPtr spotCoordFallsInMask(Coord shifted);
    Coord spotToMillerCoord(Coord xy);
    static void setupPanel(PanelPtr panel);
    static void removePanel(PanelPtr panel);
	void plotVectors(int i, PlotType plotType);
	static void plotAll(PlotType plotType);
    void stepSearch();
    static void detectorStepSearch();
    double stepScore();
    static double scoreWrapper(void *object);

    void finaliseMillerArray();
    static void finaliseMillerArrays();
    static Coord shiftForMiller(Miller *miller);
    Coord realCoordsForImageCoord(Coord xy);
    static Coord translationShiftForSpot(Spot *spot);
    static Coord swivelShiftForSpot(Spot *spot);
    static double scaleForMiller(Miller *miller);
	void print(std::ostringstream *stream);
	static void printToFile(std::string filename);
	static std::string printAll();
    
    static bool hasMillers();
    static int panelCount();
    static void clearAllMillers();
    void clearMillers();
    
    
    static void setBestShiftX(void *object, double newValue)
    {
        static_cast<Panel *>(object)->bestShift.first = newValue;
    }
    
    static void setBestShiftY(void *object, double newValue)
    {
        static_cast<Panel *>(object)->bestShift.second = newValue;
    }
    
    static double getBestShiftX(void *object)
    {
        return static_cast<Panel *>(object)->bestShift.first;
    }
    
    static double getBestShiftY(void *object)
    {
        return static_cast<Panel *>(object)->bestShift.second;
    }
    
    static double getSwivel(void *object)
    {
        return static_cast<Panel *>(object)->swivel;
    }
    
    static void setSwivel(void *object, double newValue)
    {
        Panel *panel = static_cast<Panel *>(object);
        
        panel->swivel = newValue;
        
        panel->sinSwivel = sin(panel->swivel * M_PI / 180);
        panel->cosSwivel = cos(panel->swivel * M_PI / 180);
        
    }
    
    
	static bool shouldUsePanelInfo()
	{
		return usePanelInfo;
	}

	static void setUsePanelInfo(bool useInfo)
	{
		usePanelInfo = useInfo;
	}

	const Coord& getBestShift() const
	{
		return bestShift;
	}
    
    double *pointerToBestShiftX()
    {
        return &(bestShift.first);
    }
    
    double *pointerToBestShiftY()
    {
        return &(bestShift.second);
    }

	void setBestShift(const Coord& bestShift)
	{
		this->bestShift = bestShift;
	}
    
    double getGainScale()
    {
        return gainScale;
    }
    
    static PanelPtr getPanel(int i)
    {
        return panels[i];
    }
    
    PanelTag getTag()
    {
        return tag;
    }
    
    void setTag(PanelTag newBad)
    {
        tag = newBad;
    }
    
};

#endif /* PANEL_H_ */
