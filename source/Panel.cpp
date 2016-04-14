/*
 * Panel.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#include "Panel.h"
#include <iostream>
#include "Miller.h"
#include "FileParser.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include "GraphDrawer.h"
#include "Vector.h"
#include "Spot.h"
#include "CSV.h"
#include <complex>

vector<PanelPtr> Panel::panels;
vector<PanelPtr> Panel::badPanels;
bool Panel::usePanelInfo;
double Panel::distanceMultiplier;
Coord Panel::beamCentre;

Panel::Panel(double x1, double y1, double x2, double y2, PanelTag newTag)
{
    std::vector<double> dimensions;
    dimensions.push_back(x1);
    dimensions.push_back(y1);
    dimensions.push_back(x2);
    dimensions.push_back(y2);
    
    init(dimensions, newTag);
}

Panel::Panel(vector<double> dimensions, PanelTag newTag)
{
    init(dimensions, newTag);
}

void Panel::init(vector<double> dimensions, PanelTag newTag)
{
    if (panels.size() == 0)
        usePanelInfo = true;
    
    tag = newTag;
    
    beingWrittenTo = false;
    // TODO Auto-generated constructor stub
    defaultShift = FileParser::getKey("METROLOGY_SEARCH_SIZE",
                                      METROLOGY_SEARCH_SIZE);
    bestShift = std::make_pair(0, 0);
    allowedSearchSpace = 0;
    swivel = 0;
    tilt = std::make_pair(0, 0);
    gainScale = 1;
    
    topLeft = std::make_pair(dimensions[0], dimensions[1]);
    bottomRight = std::make_pair(dimensions[2], dimensions[3]);
    
    if (dimensions.size() >= 6)
    {
        originalShift = std::make_pair(dimensions[4], dimensions[5]);
        bestShift = std::make_pair(dimensions[4], dimensions[5]);
    }
    else
    {
        usePanelInfo = false;
    }
    if (dimensions.size() >= 8)
    {
        tilt = std::make_pair(dimensions[6], dimensions[7]);
    }
    if (dimensions.size() >= 9)
    {
        swivel = dimensions[8];
    }
    
    sinSwivel = sin(swivel * M_PI / 180);
    cosSwivel = cos(swivel * M_PI / 180);
    
    if (dimensions.size() >= 10)
    {
        gainScale = dimensions[9];
    }
}

void Panel::setupPanel(PanelPtr panel)
{
    if (panel->getTag() == PanelTagBad)
    {
        badPanels.push_back(panel);
        return;
    }
    
    panels.push_back(panel);
}


void Panel::removePanel(PanelPtr panel)
{
    std::vector<PanelPtr> *panelList = &panels;
    
    if (panel->getTag() == PanelTagBad)
    {
        panelList = &badPanels;
    }
    
    for (int i = 0; i < panelList->size(); i++)
    {
        if ((*panelList)[i] == panel)
        {
            panelList->erase(panelList->begin() + i);
            return;
        }
    }

}

Panel::~Panel()
{
    millers.clear();
}

double Panel::width()
{
    return bottomRight.first - topLeft.first;
}

double Panel::height()
{
    return bottomRight.second - topLeft.second;
}

bool Panel::isCoordInPanel(Coord coord, Coord *topLeft, Coord *bottomRight)
{
    double x = coord.first;
    double y = coord.second;
    
    if (topLeft == NULL)
        topLeft = &this->topLeft;
    
    if (bottomRight == NULL)
        bottomRight = &this->bottomRight;
    
    if (x < topLeft->first || x > bottomRight->first)
        return false;
    
    if (y < topLeft->second || y > bottomRight->second)
        return false;
    /*
    logged << "Success!" << std::endl;
    sendLog();*/
    
    return true;
}

bool Panel::isMillerInPanel(Miller *miller)
{
    return isCoordInPanel(miller->getLastXY());
}

Coord Panel::midPoint()
{
    double x = topLeft.first + width() / 2;
    double y = topLeft.second + height() / 2;
    
    return std::make_pair(x, y);
}

bool Panel::addMiller(MillerPtr miller)
{
    Coord coord = std::make_pair(miller->getLastX(), miller->getLastY());
    
    if (!isCoordInPanel(coord))
        return false;
    
    boost::thread::id thread_id = boost::this_thread::get_id();
    
    if (tempMillers.count(thread_id) == 0)
        tempMillers[thread_id] = vector<MillerPtr>();
    
    tempMillers[thread_id].push_back(miller);
    
    return true;
}

void Panel::addMillerToPanelArray(MillerPtr miller)
{
    for (int i = 0; i < panels.size(); i++)
    {
        if (panels[i]->addMiller(miller))
            break;
    }
}

void Panel::finaliseMillerArrays()
{
    for (int i = 0; i < panels.size(); i++)
    {
        panels[i]->finaliseMillerArray();
    }
}

void Panel::finaliseMillerArray()
{
    unsigned int total = 0;
    
    for (std::map<boost::thread::id, vector<MillerPtr> >::iterator it = tempMillers.begin();
         it != tempMillers.end(); ++it)
    {
        total += tempMillers[it->first].size();
    }
    
    millers.reserve(total);
    
    int offset = 0;
    
    for (std::map<boost::thread::id, vector<MillerPtr> >::iterator it = tempMillers.begin();
         it != tempMillers.end(); ++it)
    {
        millers.insert(millers.begin() + offset, tempMillers[it->first].begin(), tempMillers[it->first].end());
        offset += tempMillers[it->first].size();
    }
    
}

void Panel::calculateMetrologyThread(int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < panels.size(); i += maxThreads)
    {
        panels[i]->findAllParameters();
    }
}

void Panel::print(std::ostringstream *stream)
{
    *stream << "PANEL ";
    *stream << topLeft.first << " " << topLeft.second << " " << bottomRight.first
    << " " << bottomRight.second << " ";
    *stream << bestShift.first << " " << bestShift.second << " ";
    *stream << tilt.first << " " << tilt.second << " ";
    *stream << swivel << " " << gainScale << std::endl;
    
    Logger::mainLogger->addStream(&logged);
    logged.str("");
}

void Panel::printToFile(std::string filename)
{
    std::string info = printAllThreaded();
    
    std::ofstream file;
    file.open(filename);
    file << info;
    file.close();
    
    usePanelInfo = true;
    
    Logger::mainLogger->addString("Written to file new_panels.txt");
}

std::string Panel::printAllThreaded()
{
    int totalThreads = FileParser::getMaxThreads();
    
    boost::thread_group threads;
    
    Logger::mainLogger->addString("Calculating panel shifts");
    
    for (int i = 0; i < totalThreads; i++)
    {
        boost::thread *thr = new boost::thread(calculateMetrologyThread, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    return printAll();
}

std::string Panel::printAll()
{
    std::ostringstream stream;
    
    for (int i = 0; i < panels.size(); i++)
    {
        panels[i]->print(&stream);
    }
    
    return stream.str();
}

void Panel::clearAllMillers()
{
    for (int i = 0; i < panels.size(); i++)
    {
        panels[i]->clearMillers();
    }
}

void Panel::clearMillers()
{
    millers.clear();
    vector<MillerPtr>().swap(millers);
}

int Panel::panelCount()
{
    return (int)panels.size();
}

bool scoreComparison(std::pair<Coord, double> score1,
                     std::pair<Coord, double> score2)
{
    return (score1.second > score2.second);
}

bool scoreComparisonDescending(std::pair<double, double> score1,
                               std::pair<double, double> score2)
{
    return (score1.second > score2.second);
}

PanelPtr Panel::panelForMiller(Miller *miller)
{
    for (int i = 0; i < badPanels.size(); i++)
    {
        if (badPanels[i]->isMillerInPanel(miller))
        {
            //     Logger::mainLogger->addString("Hit bad mask", LogLevelDetailed);
            return PanelPtr();
        }
    }
    
    for (int i = 0; i < panels.size(); i++)
    {
        if (panels[i]->isMillerInPanel(miller))
        {
            return panels[i];
        }
    }
    
    return PanelPtr();
}

Coord Panel::shiftSpot(Coord xy)
{
    Coord shift = bestShift;
    Coord swivelShift = getSwivelShift(xy, true);
    
    if (shift.first == FLT_MAX)
    {
        shift.first = 0;
        shift.second = 0;
    }
    
    Coord translated = std::make_pair(xy.first - shift.first, xy.second - shift.second);
    
 //   logged << "trans/swivel:\t" << bestShift.first << "\t" << bestShift.second << "\t" << swivelShift.first << "\t" << swivelShift.second << std::endl;
    sendLog();
    
    translated.first += swivelShift.first;
    translated.second += swivelShift.second;
    
    return translated;
}

PanelPtr Panel::panelForSpot(Spot *spot)
{
    Coord coord = spot->getRawXY();
    
    for (int i = 0; i < badPanels.size(); i++)
    {
        Coord shifted = badPanels[i]->shiftSpot(coord);
        
        if (badPanels[i]->isCoordInPanel(shifted))
        {
            return PanelPtr();
        }
    }
    
    for (int i = 0; i < panels.size(); i++)
    {
        Coord shifted = panels[i]->shiftSpot(coord);
      
  //      panels[i]->logged << "compare\t" << coord.first << "\t" << coord.second << "\t" << shifted.first << "\t" << shifted.second << "\t" << panels[i]->topLeft.first << "\t"  << panels[i]->topLeft.second << "\t"  << panels[i]->bottomRight.first << "\t"  << panels[i]->bottomRight.second << "\t"  << std::endl;
        panels[i]->sendLog();
        
        
        if (panels[i]->isCoordInPanel(shifted))
        {
            return panels[i];
        }
    }
    
    return PanelPtr();
    
    return panelForCoord(coord);
}

PanelPtr Panel::panelForCoord(Coord coord)
{
    for (int i = 0; i < badPanels.size(); i++)
    {
        if (badPanels[i]->isCoordInPanel(coord))
        {
       //     Logger::mainLogger->addString("Hit bad mask", LogLevelDetailed);
            return PanelPtr();
        }
    }
    
    for (int i = 0; i < panels.size(); i++)
    {
        if (panels[i]->isCoordInPanel(coord))
        {
            return panels[i];
        }
    }
    
    return PanelPtr();
}

double cartesian_to_distance(Coord dV)
{
    double distance = sqrt(pow(dV.first, 2) + pow(dV.second, 2));
    
    return distance;
}

double cartesian_to_angle(Coord dV)
{
    double angle = atan(dV.second / dV.first);
    
    if ((dV.first < 0 && dV.second > 0) ||
        (dV.first < 0 && dV.second < 0))
        angle += M_PI;
    
    if (dV.first > 0 && dV.second < 0)
        angle += M_PI * 2;
    
    return angle;
}

Coord Panel::getSwivelCoords(Coord coord)
{
    return coord;
   /*
    Coord fracCoords;
    fractionalCoordinates(miller, &fracCoords);
    
    Coord midCoords = std::make_pair(0.5, 0.5);
    // direction vector
    Coord dV = std::make_pair(fracCoords.first - midCoords.first,
                              fracCoords.second - midCoords.second);
    
    double distance = cartesian_to_distance(dV);
    double angle = cartesian_to_angle(dV);
    
    angle += swivel;
    double x = distance * cos(angle);
    double y = distance * sin(angle);
    
    double xReal = (x + 0.5) * width() + topLeft.first;
    double yReal = (y + 0.5) * height() + topLeft.second;
    
    return std::make_pair(xReal, yReal);*/
}

Coord Panel::getSwivelShift(Coord millerCoord, bool isSpot)
{
    Coord relative = relativeToMidPointForMiller(millerCoord, isSpot);
    
    if (swivel == 0)
        return std::make_pair(0, 0);
    
    double oldX = relative.first;
    double oldY = relative.second;
    
    Coord rotated = rotateCoordByAngle(relative, isSpot);
    
    rotated.first -= relative.first;
    rotated.second -= relative.second;
    
    return std::make_pair(rotated.first, rotated.second);
}

Coord Panel::rotateCoordByAngle(Coord newCoord, bool negative)
{
    double newSinSwivel = negative ? -sinSwivel : sinSwivel;
    
    double newX = cosSwivel * newCoord.first - newSinSwivel * newCoord.second;
    double newY = newSinSwivel * newCoord.first + cosSwivel * newCoord.second;

    return std::make_pair(newX, newY);
}

Coord Panel::getTiltShift(Coord millerCoord)
{
    Coord swivelCoords = millerCoord;
    Coord fracCoords;
    fractionalCoordinates(swivelCoords, &fracCoords);
    
    double xShift = tilt.first * (fracCoords.first - 0.5);
    double yShift = tilt.second * (fracCoords.second - 0.5);
    
    return std::make_pair(xShift, yShift);
}

Coord Panel::getTotalShift(Coord millerCoord, bool isMiller)
{
    Coord swivelShift = getSwivelShift(millerCoord);
    Coord tiltShift = getTiltShift(millerCoord);
    
    double x = bestShift.first + tiltShift.first + swivelShift.first;
    double y = bestShift.second + tiltShift.second + swivelShift.second;
    
    
    return std::make_pair(x, y);
}

Coord Panel::shiftForMiller(Miller *miller)
{
    PanelPtr panel = panelForMiller(miller);
    
    if (!panel)
        return std::make_pair(FLT_MAX, FLT_MAX);
    
    return panel->getTotalShift(miller->getLastXY());
}

Coord Panel::translationShiftForSpot(Spot *spot)
{
    PanelPtr panel = panelForSpot(spot);
    
    if (!panel)
        return std::make_pair(FLT_MAX, FLT_MAX);
    
    return panel->bestShift;
}

Coord Panel::swivelShiftForSpot(Spot *spot)
{
    PanelPtr panel = panelForSpot(spot);
    
    if (!panel)
        return std::make_pair(FLT_MAX, FLT_MAX);
    
    return panel->getSwivelShift(spot->getRawXY(), true);
}

double Panel::scaleForMiller(Miller *miller)
{
    PanelPtr panel = panelForMiller(miller);
    
    if (panel)
        return panel->gainScale;
    
    return 1;
}

void Panel::plotAll(PlotType plotType)
{
    for (int i = 0; i < panels.size(); i++)
    {
        panels[i]->plotVectors(i, plotType);
    }
}

Coord Panel::relativeToMidPointForMiller(Coord coord, bool isSpot)
{
    double pos_x = coord.first;
    double pos_y = coord.second;
    
    Coord panelMidPoint = midPoint();
    
    panelMidPoint.first += isSpot ? bestShift.first : 0;
    panelMidPoint.second += isSpot ? bestShift.second : 0;
    
    double rel_x = pos_x - panelMidPoint.first;
    double rel_y = pos_y - panelMidPoint.second;
    
    return std::make_pair(rel_x, rel_y);
}

double Panel::distanceFromMidPointForMiller(Miller *miller)
{
    Coord relative = relativeToMidPointForMiller(miller->getLastXY());
    
    return sqrt(relative.first * relative.first + relative.second * relative.second);
}

double Panel::angleForMiller(Miller *miller)
{
    Coord relative = relativeToMidPointForMiller(miller->getLastXY());
    
    double rel_x = relative.first;
    double rel_y = relative.second;
    
    double angle = atan(rel_y / rel_x) * 180 / M_PI;
    
    if ((rel_x < 0 && rel_y > 0) || (rel_x < 0 && rel_y < 0))
    {
        angle += 180;
    }
    else if (rel_x > 0 && rel_y < 0)
    {
        angle += 360;
    }
    
    return angle;
}

void Panel::plotVectors(int i, PlotType plotType)
{
#ifdef MAC
    bool splitResolution = Panel::panelCount() == 1;
    
    vector<double> bins;
    if (splitResolution)
    {
        StatisticsManager::generateResolutionBins(0, 1.0, 4, &bins);
    }
    else
    {
        StatisticsManager::generateResolutionBins(0, 1.0, 1, &bins);
    }
    
    for (int j = 0; j < bins.size(); j++)
    {
        if (bins.size() == 2 && j == 1)
            break;

   /*     vector<double> xRed, yRed;
        vector<double> xBlue, yBlue;
    */
        vector<MillerPtr> resMillers;
    
        
        double minD = 0;
        double maxD = 0;
        
        if (j == bins.size() - 1)
        {
            minD = 0;
            maxD = FLT_MAX;
        }
        else
        {
            minD = bins[j] == 0 ? 0 : 1 / bins[j];
            maxD = bins[j + 1] == 0 ? FLT_MAX : 1 / bins[j + 1];
        }
        
        std::cout << minD << ", " << maxD << std::endl;
        
        for (int k = 0; k < millers.size(); k++)
        {
            MillerPtr miller = millers[k];
            
            if (miller->getResolution() > minD && miller->getResolution() < maxD)
            {
                resMillers.push_back(miller);
            }
        }
        
        double aveIntensity = Miller::averageRawIntensity(resMillers);
        
        std::cout << "Average intensity: " << aveIntensity << std::endl;
        
        double minX = FLT_MAX; double minY = FLT_MAX; double maxX = -FLT_MAX; double maxY = -FLT_MAX;
        
        CSV csv = CSV(6, "shift_x", "shift_y", "intensity", "rel_x", "rel_y", "angle");
        
        for (int k = 0; k < resMillers.size(); k++)
        {
            MillerPtr miller = resMillers[k];
            
            double intensity = miller->getRawIntensity();
            
            Coord shift = miller->getShift();
            Coord expectedShift = this->getTotalShift((&*miller)->getLastXY());
         //   shift.first -= originalShift.first; // plot relative shift away from the last round of refinement
         //   shift.second -= originalShift.second;
            
            if (shift.first < minX)
                minX = shift.first;
            if (shift.first > maxX)
                maxX = shift.first;
            if (shift.second < minY)
                minY = shift.second;
            if (shift.second < maxY)
                maxY = shift.second;
            
            Coord difference;
            
            difference.first = shift.first;
            difference.second = shift.second;
            
            if (plotType == PlotTypeRelative)
            {
            difference.first = expectedShift.first - shift.first;
            difference.second = expectedShift.second - shift.second;
            }
            
            double pos_x = miller->getLastX();
            double pos_y = miller->getLastY();
            
            Coord panelMidPoint = midPoint();
            
            double rel_x = pos_x - panelMidPoint.first;
            double rel_y = pos_y - panelMidPoint.second;
            
            double angle = angleForMiller(&*miller);
            
            double shiftAngle = cartesian_to_angle(difference);
            shiftAngle *= 180 / M_PI;
            
            bool under = intensity < aveIntensity * 2.5;
            under = (miller->getRawIntensity() / miller->getCountingSigma() < 17);
            double strength = miller->getRawIntensity();// / miller->getCountingSigma();
            
            csv.addEntry(0, difference.first, difference.second, strength, rel_x, rel_y, angle);
            
      /*      xRed.push_back(difference.first);
            yRed.push_back(difference.second);
            xBlue.push_back(strength);
            yBlue.push_back(0); */
            
        }
        
        std::ostringstream resString;
        if (j < bins.size())
        {
            resString << bins[j];
        }
        else
        {
            resString << "all";
        }
        
        std::ostringstream filename_stream;
        filename_stream << std::string("panel_plot_") << i << "_" << j << "_" << resString.str() + ".csv";
        std::string filename = filename_stream.str();
        
        csv.writeToFile(filename);
        
        resMillers.clear();
        vector<MillerPtr>().swap(resMillers);
        /*
        vector<vector<double> > xs;
        vector<vector<double> > ys;
        xs.push_back(xRed); // actually green
        ys.push_back(yRed);
        xs.push_back(xBlue); // actually blue
        ys.push_back(yBlue);
        
        GraphDrawer drawer = GraphDrawer(MtzManager::getReferenceManager());
        GraphMap map = GraphMap();
        
        std::ostringstream title_stream;
        title_stream << "Panel plot ";
        title_stream << topLeft.first << " " << topLeft.second << " "
        << bottomRight.first << " " << bottomRight.second;
        std::string title = title_stream.str();
        
        map["title"] = title;
        map["xTitle"] = "X shift";
        map["yTitle"] = "Y shift";
        map["plotType"] = "point";
        
        if (plotType == PlotTypeRelative)
        {
            map["xMin"] = minX;
            map["yMin"] = minY;
            map["xMax"] = maxX;
            map["yMax"] = maxY;
        }
        map["colour_0"] = 3;
        map["colour_1"] = 9;
        
        std::ostringstream resString;
        if (j < bins.size())
            resString << bins[j];
        else
            resString << "all";
        
        std::ostringstream filename_stream;
        filename_stream << std::string("panel_plot_") << i << "_" << j << "_" << resString.str();
        std::string filename = filename_stream.str();
        
        drawer.plot(filename, map, xs, ys);*/
    }
#endif
}

double Panel::scoreBetweenResolutions(double minRes, double maxRes)
{
    double score = 0;
    
    for (int i = 0; i < panelCount(); i++)
    {
        score += panels[i]->stdevScore(minRes, maxRes);
    }
    
    return score;
}

double Panel::stdevScore(double minRes, double maxRes)
{
    vector<MillerPtr> resMillers;
    
    double minD, maxD;
    
    StatisticsManager::convertResolutions(minRes, maxRes, &minD, &maxD);
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        
        if (miller->getResolution() > minD && miller->getResolution() < maxD)
        {
            resMillers.push_back(miller);
        }
    }
    
    double aveIntensity = Miller::averageRawIntensity(resMillers);
    vector<double> xShifts, yShifts;
    
    for (int i = 0; i < resMillers.size(); i++)
    {
        MillerPtr miller = resMillers[i];
        
        if (miller->getRawIntensity() > aveIntensity)
        {
            Coord shift = miller->getShift();
            
            xShifts.push_back(shift.first);
            yShifts.push_back(shift.second);
        }
    }
    
    double xStdev = standard_deviation(&xShifts);
    double yStdev = standard_deviation(&yShifts);
    
    return xStdev + yStdev;
}

double Panel::scoreDetectorDistance(void *object)
{
    vector<double> xShifts, yShifts;
    
    double squares = 0;
    
    vector<double> oldBeamCentre = FileParser::getKey("BEAM_CENTRE", vector<double>());
    
    if (oldBeamCentre.size() == 0)
    {
        oldBeamCentre.push_back(BEAM_CENTRE_X);
        oldBeamCentre.push_back(BEAM_CENTRE_Y);
    }
    
//    Coord beamShift = std::make_pair(beamCentre.first - oldBeamCentre[0], beamCentre.second - oldBeamCentre[1]);
    
    for (int i = 0; i < panels.size(); i++)
    {
        Coord bestShift = panels[i]->getBestShift();
        bestShift = std::make_pair(-bestShift.first, -bestShift.second);
        Coord midPoint = panels[i]->midPoint();
        bestShift.first += midPoint.first;
        bestShift.second += midPoint.second;
        
        bestShift.first -= beamCentre.first;
        bestShift.second -= beamCentre.second;
        
        bestShift.first *= distanceMultiplier;
        bestShift.second *= distanceMultiplier;
        
        bestShift.first -= midPoint.first;
        bestShift.second -= midPoint.second;
        
        bestShift.first += beamCentre.first;
        bestShift.second += beamCentre.second;
        
        xShifts.push_back(bestShift.first);
        yShifts.push_back(bestShift.second);
        
        squares += pow(bestShift.first, 2);
        squares += pow(bestShift.second, 2);
    }
    
    double stdev = standard_deviation(&xShifts);
    stdev += standard_deviation(&yShifts);
    
    std::cout << "standard_deviation: " << stdev << std::endl;
    
    return stdev;
}

void Panel::refineDetectorDistance()
{
    double ddStep = 0.5;
    int count = 0;
    
    while (ddStep > 0.0001 && count < 20)
    {
        minimizeParameter(ddStep, &distanceMultiplier, scoreDetectorDistance, NULL);
        
        count++;
    }
    
    std::ostringstream logged;
    logged << "New detector distance multiplier: " << distanceMultiplier << std::endl;
    
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    
    double newDistance = detectorDistance * distanceMultiplier;
    
    logged << "Old detector distance: " << detectorDistance << std::endl;
    logged << "New detector distance: " << newDistance << std::endl;
    
    FileParser::setKey("DETECTOR_DISTANCE", newDistance);
    
    Logger::mainLogger->addStream(&logged);
}

void Panel::expectedBeamCentre()
{
    vector<double> newBeamCentre = FileParser::getKey("BEAM_CENTRE", vector<double>());
    
    if (newBeamCentre.size() == 0)
    {
        newBeamCentre.push_back(BEAM_CENTRE_X);
        newBeamCentre.push_back(BEAM_CENTRE_Y);
    }
    
    double totalShiftX = 0;
    double totalShiftY = 0;
    int num = 0;
    
    for (int i = 0; i < panels.size(); i++)
    {
        if (panels[i]->millers.size() == 0)
            continue;
        
        std::ostringstream progress;
        progress << "Calculating position for Panel " << i << std::endl;
        Logger::mainLogger->addStream(&progress);
        
        num++;
        
        panels[i]->findShift(2, 0.4);
        
        Coord bestShift = panels[i]->getBestShift();
        totalShiftX += bestShift.first;
        totalShiftY += bestShift.second;
    }
    
    totalShiftX /= num;
    totalShiftY /= num;
    
    std::ostringstream logged;
    logged << "Original beam centre: " << newBeamCentre[0] <<
    "\t" << newBeamCentre[1] << std::endl;
    
    newBeamCentre[0] += totalShiftX;
    newBeamCentre[1] += totalShiftY;
    
    logged << "New beam centre: " << newBeamCentre[0] << "\t"
    << newBeamCentre[1] << std::endl;
    
    beamCentre = std::make_pair(newBeamCentre[0], newBeamCentre[1]);
  
    vector<double> centreVector = vector<double>();
    centreVector.push_back(beamCentre.first);
    centreVector.push_back(beamCentre.second);
    
    FileParser::setKey("BEAM_CENTRE", centreVector);
    
    logged << "Please change your input script to new beam centre." << std::endl;
    
    
    
    Logger::mainLogger->addStream(&logged);
}

void Panel::fractionalCoordinates(Miller *miller, Coord *frac)
{
    Coord position = std::make_pair(miller->getLastX(),
                               miller->getLastY());
    
    fractionalCoordinates(position, frac);
}

void Panel::fractionalCoordinates(Coord coord, Coord *frac)
{
    double fractionX = (coord.first - topLeft.first)
    / width();
    
    double fractionY = (coord.second - topLeft.second)
    / height();
    
    *frac = std::make_pair(fractionX, fractionY);
}

void Panel::findAllParameters()
{
    if (millers.size() == 0)
    {
        return;
    }
    
    findShift(2, 0.4);
    this->findAxisDependence(3);
//    refineAllParameters(3);
//    findShift(2, 0.1);
}

void Panel::findShift(double windowSize, double step, double x, double y)
{
    std::ostringstream logged;
    logged << "**** NEW PANEL ****" << std::endl;
    
    vector<std::pair<Coord, double> > scores;
    logged << "Miller count: " << millers.size() << std::endl;
    
    double minX = -defaultShift + originalShift.first + x + 0.5;
    double maxX = defaultShift + originalShift.first + x - 0.5;
    double minY = -defaultShift + originalShift.second + y + 0.5;
    double maxY = defaultShift + originalShift.second + y - 0.5;
    
    double intensityThreshold = FileParser::getKey("INTENSITY_THRESHOLD", 12.);
    
    logged << "Finding best shift within window x(" << minX << ", " << maxX << "), y(" << minY << ", " << maxY << ")" << std::endl;
    
    for (double i = minX; i < maxX; i += step)
    {
        for (double j = minY; j < maxY; j += step)
        {
            Coord windowTopLeft = std::make_pair(i, j);
            Coord windowBottomRight = std::make_pair(i + windowSize, j + windowSize);
            
            Coord windowMidPoint = std::make_pair(i + windowSize / 2,
                                             j + windowSize / 2);
            
            double score = 0;
            
            for (int k = 0; k < millers.size(); k++)
            {
                if (!millers[k])
                {
                    continue;
                }
                
                Coord predictedPosition = millers[k]->getLastXY();
                Coord currentShift = bestShift;
                Coord newShift = rotateCoordByAngle(millers[k]->getShift(), true);
                newShift = millers[k]->getShift();
                
                
                Coord millerShift = std::make_pair(newShift.first + currentShift.first,
                                                   newShift.second + currentShift.second);
                
                bool inWindow = isCoordInPanel(millerShift, &windowTopLeft,
                                               &windowBottomRight);
                
                if (!inWindow)
                {
                    continue;
                }
                
                bool strong = millers[k]->getRawIntensity() / millers[k]->getCountingSigma() > intensityThreshold;
                
                if (strong)
                {
                    score++;
                }
            }
            
            scores.push_back(make_pair(windowMidPoint, score));
        }
    }
    
    std::sort(scores.begin(), scores.end(), scoreComparison);
    
    logged << "Changed best shift to " << scores[0].first.first << "\t" << scores[0].first.second << " with score of " << scores[0].second << std::endl;
    Logger::mainLogger->addStream(&logged);
    
    bestShift = scores[0].first;
}

double Panel::swivelShiftScoreWrapper(void *object)
{
    double threshold = 200;
    
    Panel *panel = static_cast<Panel *>(object);
    
    std::vector<double> shiftXs, shiftYs;
    
    for (int i = 0; i < panel->millers.size(); i++)
    {
        if (panel->millers[i]->getRawIntensity() > threshold)
        {
            MillerPtr miller = panel->millers[i];
            
            Coord shift = shiftForMiller(&*miller);
            
            shiftXs.push_back(shift.first);
            shiftYs.push_back(shift.second);
        }
    }
    
    double stdevX = standard_deviation(&shiftXs);
    double stdevY = standard_deviation(&shiftYs);
    
    double score = stdevX + stdevY;
    
 //   std::ostringstream logged;
 //   logged << "Swivel score: " << score << " for " << panel->swivel << std::endl;
 //   Logger::mainLogger->addStream(&logged);
    
    return score;
}

double Panel::tiltShiftScoreWrapper(void *object)
{
    return static_cast<Panel *>(object)->tiltShiftScore();
}

double Panel::tiltShiftScore(double stdev)
{
    double xMin = bestShift.first - tiltWindowSize / 2;
    double xMax = bestShift.first + tiltWindowSize / 2;
    
    double yMin = bestShift.second - tiltWindowSize / 2;
    double yMax = bestShift.second + tiltWindowSize / 2;
    
    vector<double> tiltShifts;
    vector<double> weights;
    
    for (int k = 0; k < millers.size(); k++)
    {
        Coord shift = millers[k]->getShift();
        
        if (shift.first < xMin || shift.first > xMax)
            continue;
        
        if (shift.second < yMin || shift.second > yMax)
            continue;
        
        if (millers[k]->getRawIntensity() < averageIntensity)
            continue;
        
        double rawIntensity = millers[k]->getRawIntensity();
        
        Coord totalShift = this->shiftForMiller(&*millers[k]);
        Coord shiftDifference = std::make_pair(totalShift.first - shift.first,
                                               totalShift.second - shift.second);
        
        double desiredAxis = tiltHorizontalAxis ? shiftDifference.first : shiftDifference.second;
        
        tiltShifts.push_back(desiredAxis);
        weights.push_back(log(rawIntensity));
    }
    
    double stdeviation = standard_deviation(&tiltShifts);
    double mean = weighted_mean(&tiltShifts, &weights);
    
    return stdev ? stdeviation : mean;
}

void Panel::refineAllParameters(double windowSize)
{
    double swivelStep = 0.2;
    bool minimised = false;
    int count = 0;
    std::ostringstream logged;
    
    double before = swivelShiftScoreWrapper(this);
    
    while (count < 20 && !minimised)
    {
        double swivelScore = minimizeParameter(swivelStep, &this->swivel, swivelShiftScoreWrapper, this);
        
        logged << swivel << ", " << swivelScore << std::endl;
        
        if (swivelStep < 0.001)
            break;
        
        count++;
    }

    double after = swivelShiftScoreWrapper(this);

    
    logged << "Swivel angle " << swivel << " - score " << before << ", after refinement " << after << "." << std::endl;
    Logger::mainLogger->addStream(&logged);
    
/*    double xStep = 0.2;
    double yStep = 0.2;
    bool minimized = false;
    int count = 0;
    
    tiltWindowSize = windowSize;
    
    logged << "Tilt refinement progress: X\tY\tsd(X)\tsd(Y)" << std::endl;
    
    averageIntensity = Miller::averageRawIntensity(millers);
    
    while (count < 25 && !minimized)
    {
        tiltHorizontalAxis = true;
        double hozScore = minimizeParameter(xStep, &this->tilt.first, tiltShiftScoreWrapper, this);
        
        tiltHorizontalAxis = false;
        double vertScore = minimizeParameter(yStep, &this->tilt.second, tiltShiftScoreWrapper, this);
        
        logged << tilt.first << "\t" << tilt.second << "\t" << hozScore << "\t" << vertScore << std::endl;
        
        if (xStep < 0.01 && yStep < 0.01)
            minimized = true;
        
        count++;
    }
    
    logged << std::endl;
    
    minimized = false;
    count = 0;*/
}

double Panel::detectorGain(double *error)
{
    vector<double> observedPartialities, observedWeights;
    MtzManager *reference = MtzManager::getReferenceManager();
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (!millers[i]->accepted())
            continue;
        
        count++;
        double observed = millers[i]->observedPartiality(reference);
        double intensity = millers[i]->intensity();
        
        if (observed == observed && std::isfinite(observed) && intensity == intensity)
        {
            observedPartialities.push_back(observed);
            observedWeights.push_back(fabs(intensity));
        }
    }
    
    double mean = weighted_mean(&observedPartialities, &observedWeights);
    *error = standard_deviation(&observedPartialities);
    
    gainScale = 1 / mean;
    
    return mean;
}

void Panel::checkDetectorGains()
{
    std::ostringstream logged;
    
    for (int i = 0; i < panels.size(); i++)
    {
        double error;
        double mean = panels[i]->detectorGain(&error);
        double millerCount = panels[i]->millers.size();
        
        logged << i << "\t" << mean << "\t" << error << "\t" << millerCount << std::endl;
    }
    
    std::cout << logged.str();
    
    Logger::mainLogger->addStream(&logged);
}

bool Panel::hasMillers()
{
    for (int i = 0; i < panels.size(); i++)
    {
        if (panels[i]->millers.size() > 0)
            return true;
    }
    
    return false;
}

void Panel::findAxisDependence(double windowSize)
{
    double fracStep = 0.1;
    
    for (int i = 0; i < 2; i++)
    {
        bool refineX = (i == 0);
        
        vector<double> fracs;
        vector<double> scores;
        
        double xMin = bestShift.first - windowSize / 2;
        double xMax = bestShift.first + windowSize / 2;
        
        double yMin = bestShift.second - windowSize / 2;
        double yMax = bestShift.second + windowSize / 2;
        
        double aveIntensity = Miller::averageRawIntensity(millers);
        
        for (double frac = 0; frac < 1; frac += fracStep)
        {
            double score = 0;
            double weights = 0;
            
            for (int k = 0; k < millers.size(); k++)
            {
                if (!millers[k])
                    continue;
                
                double rawIntensity = millers[k]->getRawIntensity();
                
                if (rawIntensity < aveIntensity)
                    continue;
                
                Coord fracCoords;
                Coord shift = millers[k]->getShift();
                Coord position = millers[k]->position();
                
                fractionalCoordinates(position, &fracCoords);
                
                double observedFrac = refineX ? fracCoords.first : fracCoords.second;
                
                if (observedFrac > frac && observedFrac < frac + fracStep)
                {
                    if (shift.first < xMin || shift.first > xMax)
                        continue;
                    
                    if (shift.second < yMin || shift.second > yMax)
                        continue;
                    
                    double weight = log(rawIntensity);
                    
                    double addition = refineX ? shift.first - bestShift.first : shift.second - bestShift.second;
                    addition *= weight;
                    score += addition;
                    weights += weight;
                }
            }
            
            if (weights > 0)
                score /= weights;
            
            fracs.push_back(frac - 0.5);
            scores.push_back(score);
        }
        
        double grad = gradient_between_vectors(&fracs, &scores);
        
        if (refineX)
            tilt.first = grad;
        else
            tilt.second = grad;
    }
}
