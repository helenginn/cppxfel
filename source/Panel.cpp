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
#include "GetterSetterMap.h"

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
    else if (tag != PanelTagBad)
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
    
    panel->setPanelNum((int)panels.size());
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
    
    translated.first += swivelShift.first;
    translated.second += swivelShift.second;
    
    return translated;
}

PanelPtr Panel::spotCoordFallsInMask(Coord shifted)
{
    for (int j = 0; j < badPanels.size(); j++)
    {
        if (badPanels[j]->isCoordInPanel(shifted))
        {
            return badPanels[j];
        }
    }

    return PanelPtr();
}

PanelPtr Panel::panelForSpotCoord(Coord coord, PanelPtr *anyBadPanel)
{
    for (int i = 0; i < panels.size(); i++)
    {
        Coord shifted = panels[i]->shiftSpot(coord);
        
        if (panels[i]->isCoordInPanel(shifted))
        {
            PanelPtr mask = spotCoordFallsInMask(shifted);
            
            if (mask)
            {
                if (anyBadPanel) *anyBadPanel = mask;
                return PanelPtr();
            }
            
            if (anyBadPanel) *anyBadPanel = PanelPtr();
            return panels[i];
        }
    }
    
    if (anyBadPanel) *anyBadPanel = PanelPtr();
    return PanelPtr();
}

PanelPtr Panel::panelForSpot(Spot *spot)
{
    Coord coord = spot->getRawXY();
    
    return panelForSpotCoord(coord);
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

Coord Panel::getSwivelShift(Coord millerCoord, bool isSpot)
{
    Coord relative = relativeToMidPointForMiller(millerCoord, isSpot);
    
    if (swivel == 0)
        return std::make_pair(0, 0);
    
    double oldX = relative.first;
    double oldY = relative.second;
    
    double newSinSwivel = isSpot ? -sinSwivel : sinSwivel;
    
    double newX = cosSwivel * oldX - newSinSwivel * oldY;
    double newY = newSinSwivel * oldX + cosSwivel * oldY;
    
    newX -= relative.first;
    newY -= relative.second;
    
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

Coord Panel::realCoordsForImageCoord(Coord xy)
{
    Coord shift = bestShift;
    Coord swivelShift = getSwivelShift(xy, true);
    
    Coord translated = std::make_pair(xy.first - shift.first, xy.second - shift.second);
    translated.first += swivelShift.first;
    translated.second += swivelShift.second;
    
    return translated;
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
            
            double rel_x = pos_x - panelMidPoint.first + expectedShift.first;
            double rel_y = pos_y - panelMidPoint.second + expectedShift.second;
            
            double angle = angleForMiller(&*miller);
            
            double shiftAngle = cartesian_to_angle(difference);
            shiftAngle *= 180 / M_PI;
            
            bool under = intensity < aveIntensity * 2.5;
            under = (miller->getRawIntensity() / miller->getCountingSigma() < 17);
            double strength = miller->getRawIntensity();// / miller->getCountingSigma();
            
            bool isStrong = IOMRefiner::millerReachesThreshold(miller);
            
            if (isStrong)
            {
                csv.addEntry(0, difference.first, difference.second, strength, rel_x, rel_y, angle);
            }
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
        
        logged << "Panel plot for (" << topLeft.first << ", " << topLeft.second << ") to ("
        << bottomRight.first << ", " << bottomRight.second << ")" << std::endl;
        sendLog();
        
        csv.plotColumns(0, 1);
        
        csv.writeToFile(filename);
        
        resMillers.clear();
        vector<MillerPtr>().swap(resMillers);
    
    }
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
    /*
    findShift(2, 0.4);
    this->findAxisDependence(3);*/
    
    this->stepSearch();
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
                Coord newShift = millers[k]->getShift();
                
                
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

bool Panel::hasMillers()
{
    for (int i = 0; i < panels.size(); i++)
    {
        if (panels[i]->millers.size() > 0)
            return true;
    }
    
    return false;
}

double Panel::scoreWrapper(void *object)
{
    return static_cast<Panel *>(object)->stepScore();
}

double Panel::stepScore()
{
    double intensityThreshold = FileParser::getKey("INTENSITY_THRESHOLD", 12.0);
    std::vector<double> shiftXs;
    std::vector<double> shiftYs;
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];

        if (miller->getRawIntensity() / miller->getRawCountingSigma() < intensityThreshold)
            continue;
        
        int x, y;
        
        count++;
        
        miller->positionOnDetector(miller->getMatrix(), &x, &y);
        
        double shiftX = miller->getShift().first;
        double shiftY = miller->getShift().second;
        
        shiftXs.push_back(shiftX);
        shiftYs.push_back(shiftY);
    }
    
    double stdevX = standard_deviation(&shiftXs, NULL, 0);
    double stdevY = standard_deviation(&shiftYs, NULL, 0);
    double aScore = stdevX + stdevY;
    
//    logged << "Panel " << panelNum << " on st dev " << aScore << std::endl;
//   sendLog();
    
    return aScore;
}

void Panel::stepSearch()
{
    usePanelInfo = true;
    
    GetterSetterMapPtr refinementMap = GetterSetterMapPtr(new GetterSetterMap());
    
    refinementMap->addParameter(this, getBestShiftX, setBestShiftX, 2.0, 0.2);
    refinementMap->addParameter(this, getBestShiftY, setBestShiftY, 2.0, 0.2);

    refinementMap->addParameter(this, getSwivel, setSwivel, 0.4, 0.01);
    
    refinementMap->setEvaluationFunction(scoreWrapper, this);
    refinementMap->setCycles(30);
    
    refinementMap->refine(GetterSetterStepSearch);
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
