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
#include "FileReader.h"

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
    
    Coord origTopLeft = millerToSpotCoord(topLeft);
    logged << "Original top left coord: " << origTopLeft.first << ", " << origTopLeft.second << std::endl;
    sendLog();
    
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

Coord Panel::midPoint()
{
    double x = topLeft.first + width() / 2;
    double y = topLeft.second + height() / 2;
    
    return std::make_pair(x, y);
}

bool Panel::isMillerCoordInPanel(Coord coord, Coord *topLeft, Coord *bottomRight)
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
    
    return true;
}

bool Panel::isMillerInPanel(Miller *miller)
{
    return isMillerCoordInPanel(miller->getLastXY());
}

bool Panel::addMiller(MillerPtr miller)
{
    Coord coord = std::make_pair(miller->getLastX(), miller->getLastY());
    
    if (!isMillerCoordInPanel(coord))
        return false;
    
    boost::thread::id thread_id = boost::this_thread::get_id();
    
    if (tempMillers.count(thread_id) == 0)
        tempMillers[thread_id] = vector<MillerPtr>();
    
    tempMillers[thread_id].push_back(miller);
    
    return true;
}

void Panel::addMillerToPanelArray(MillerPtr miller)
{
    if (!miller->reachesThreshold())
        return;
    
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
    std::string dirFilename = FileReader::addOutputDirectory(filename);
    
    std::ofstream file;
    file.open(dirFilename);
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

int Panel::panelCount()
{
    return (int)panels.size();
}

PanelPtr Panel::panelForMiller(Miller *miller)
{
    PanelPtr currentPanel = miller->getPanel();
    
    if (currentPanel) return currentPanel;
    
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
            miller->setPanel(panels[i]);
            
            return panels[i];
        }
    }
    
    return PanelPtr();
}

Coord Panel::spotToMillerCoord(Coord xy)
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
        if (badPanels[j]->isMillerCoordInPanel(shifted))
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
        Coord shifted = panels[i]->spotToMillerCoord(coord);
        
        if (panels[i]->isMillerCoordInPanel(shifted))
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
    
    PanelPtr panel = panelForSpotCoord(coord);
    
    return panel;
}

Coord Panel::getSwivelShift(Coord millerCoord, bool isSpot)
{
    Coord relative = relativeToMidPointForMiller(millerCoord, isSpot);
    
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

Coord Panel::millerToSpotCoordShift(Coord millerCoord)
{
    Coord swivelShift = getSwivelShift(millerCoord);
    Coord tiltShift = getTiltShift(millerCoord);
    
    double x = bestShift.first + tiltShift.first + swivelShift.first;
    double y = bestShift.second + tiltShift.second + swivelShift.second;
    
    return std::make_pair(x, y);
}

Coord Panel::millerToSpotCoord(Coord millerCoord)
{
    Coord theShift = millerToSpotCoordShift(millerCoord);
    
    return std::make_pair(millerCoord.first + theShift.first, millerCoord.second + theShift.second);
}

Coord Panel::shiftForMiller(Miller *miller)
{
    PanelPtr panel = panelForMiller(miller);
    
    if (!panel)
        return std::make_pair(FLT_MAX, FLT_MAX);
    
    return panel->millerToSpotCoordShift(miller->getLastXY());
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

void Panel::plotAll(PlotType plotType)
{
    for (int i = 0; i < panels.size(); i++)
    {
        panels[i]->plotVectors(i, plotType);
    }
}

Coord Panel::relativeToCoordForMiller(Coord coord, Coord relative, bool isSpot)
{
    double pos_x = coord.first;
    double pos_y = coord.second;
    
    relative.first += isSpot ? bestShift.first : 0;
    relative.second += isSpot ? bestShift.second : 0;
    
    double rel_x = pos_x - relative.first;
    double rel_y = pos_y - relative.second;
    
    return std::make_pair(rel_x, rel_y);
}

Coord Panel::relativeToMidPointForMiller(Coord coord, bool isSpot)
{
    Coord panelMidPoint = midPoint();
    return relativeToCoordForMiller(coord, panelMidPoint, isSpot);
}

Coord Panel::relativeToTopLeftForMiller(Coord coord, bool isSpot)
{
    return relativeToCoordForMiller(coord, topLeft, isSpot);
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
        
        for (int k = 0; k < millers.size(); k++)
        {
            MillerPtr miller = millers[k];
            
            if (miller->getResolution() > minD && miller->getResolution() < maxD)
            {
                resMillers.push_back(miller);
            }
        }
        
        double aveIntensity = Miller::averageRawIntensity(resMillers);
        
        logged << "Average intensity: " << aveIntensity << std::endl;
        sendLog();
        
        double minX = FLT_MAX; double minY = FLT_MAX; double maxX = -FLT_MAX; double maxY = -FLT_MAX;
        
        CSV csv = CSV(6, "shift_x", "shift_y", "intensity", "rel_x", "rel_y", "angle");
        
        for (int k = 0; k < resMillers.size(); k++)
        {
            MillerPtr miller = resMillers[k];
            
            double intensity = miller->getRawIntensity();
            
            Coord shift = miller->getShift();
            Coord expectedShift = this->millerToSpotCoordShift((&*miller)->getLastXY());
            
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
            
            double shiftAngle = cartesian_to_angle(difference.first, difference.second);
            shiftAngle *= 180 / M_PI;
            
            bool under = intensity < aveIntensity * 2.5;
            under = (miller->getRawIntensity() / miller->getCountingSigma() < 17);
            double strength = miller->getRawIntensity();// / miller->getCountingSigma();
            
            bool isStrong = miller->reachesThreshold();
            
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
    
    this->stepSearch();
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

double Panel::globalPseudoScoreWrapper(void *object)
{
    IndexManager *manager = static_cast<IndexManager *>(object);
    
    double score = IndexManager::pseudoScore(manager);
    
    return score;
}

double Panel::globalScoreWrapper(void *object)
{
    double scoreTotal = 0;
    
    for (int i = 0; i < panels.size(); i++)
    {
        panels[i]->refreshMillerPositions();
        scoreTotal += panels[i]->stepScore();
    }
    
    scoreTotal /= (double)panels.size();
    
    return scoreTotal;
}

double Panel::stepScore()
{
    std::vector<double> shiftXs;
    std::vector<double> shiftYs;
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        MillerPtr miller = millers[i];
        
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
    
    return aScore;
}

double Panel::scoreWrapper(void *object)
{
    return static_cast<Panel *>(object)->stepScore();
}

void Panel::stepSearch()
{
    usePanelInfo = true;
    
    Coord origBottomRight = millerToSpotCoord(bottomRight);
    
    Coord origTopLeft = millerToSpotCoord(topLeft);
    
    GetterSetterMapPtr refinementMap = GetterSetterMapPtr(new GetterSetterMap());
    
    refinementMap->addParameter(this, getBestShiftX, setBestShiftX, 2.0, 0.08);
    refinementMap->addParameter(this, getBestShiftY, setBestShiftY, 2.0, 0.08);

    refinementMap->addParameter(this, getSwivel, setSwivel, 0.4, 0.01);
    
    refinementMap->setEvaluationFunction(scoreWrapper, this);
    refinementMap->setCycles(30);
    
    refinementMap->refine(GetterSetterStepSearch);
    
    Coord newTopLeft = spotToMillerCoord(origTopLeft);
    Coord newBottomRight = spotToMillerCoord(origBottomRight);
    
    double topLeftCorrectionX = newTopLeft.first - topLeft.first;
    double topLeftCorrectionY = newTopLeft.second - topLeft.second;
    
    double bottomRightCorrectionX = newBottomRight.first - bottomRight.first;
    double bottomRightCorrectionY = newBottomRight.second - bottomRight.second;
    
    topLeft.first -= topLeftCorrectionX;
    topLeft.second -= topLeftCorrectionY;
    
    bottomRight.first -= bottomRightCorrectionX;
    bottomRight.second -= bottomRightCorrectionY;
}

void Panel::refreshMillerPositions()
{
    for (int i = 0; i < millers.size(); i++)
    {
        int xTmp, yTmp;
        millers[i]->positionOnDetector(MatrixPtr(), &xTmp, &yTmp, false);
    }
}

void Panel::detectorStepSearch(bool pseudo, std::vector<ImagePtr> images)
{
    std::ostringstream logged;
    IndexManager *powderManager = NULL;
    void *object = NULL;
    Getter getter = globalScoreWrapper;
    std::string absRel = pseudo ? "relative" : "absolute";
    
    
    logged << "*********************************" << std::endl;
    logged << "***   Detector step search    ***" << std::endl;
    logged << "***  Beam centre vs " << absRel << "  ***" << std::endl;
    logged << "*********************************" << std::endl;
    
    Logger::log(logged);
    
    GetterSetterMapPtr beamRefinementMap = GetterSetterMapPtr(new GetterSetterMap());
    
    beamRefinementMap->addParameter(NULL, Image::getGlobalBeamX, Image::setGlobalBeamX, 2.0, 0.05);
    beamRefinementMap->addParameter(NULL, Image::getGlobalBeamY, Image::setGlobalBeamY, 2.0, 0.05);
    beamRefinementMap->setVerbose(true);
    
    beamRefinementMap->setEvaluationFunction(globalScoreWrapper, NULL);
    
    beamRefinementMap->setCycles(30);
    
    beamRefinementMap->refine(GetterSetterStepSearch);
    
    if (pseudo)
    {
        powderManager = new IndexManager(images);
        object = powderManager;
        getter = globalPseudoScoreWrapper;
        
        logged << "*********************************" << std::endl;
        logged << "***   Detector step search    ***" << std::endl;
        logged << "***  D. distance vs " << absRel << "  ***" << std::endl;
        logged << "*********************************" << std::endl;
        
        Logger::log(logged);
        
        GetterSetterMapPtr distRefinementMap = GetterSetterMapPtr(new GetterSetterMap());
        
        distRefinementMap->addParameter(NULL, Image::getGlobalDetectorDistance, Image::setGlobalDetectorDistance, 5.0, 0.01);
        distRefinementMap->setVerbose(true);
        
        distRefinementMap->setEvaluationFunction(getter, object);
        
        distRefinementMap->setCycles(30);
        
        distRefinementMap->refine(GetterSetterStepSearch);
        
    }
    
    /*
    logged << "*********************************" << std::endl;
    logged << "****  Detector step search  *****" << std::endl;
    logged << "**** Distance + beam centre *****" << std::endl;
    logged << "*********************************" << std::endl;
    
    Logger::log(logged);
    
    
    GetterSetterMapPtr totalRefinementMap = GetterSetterMapPtr(new GetterSetterMap());
    
    totalRefinementMap->addParameter(NULL, Image::getGlobalBeamX, Image::setGlobalBeamX, 1.0, 0.01);
    totalRefinementMap->addParameter(NULL, Image::getGlobalBeamY, Image::setGlobalBeamY, 1.0, 0.01);
    totalRefinementMap->addParameter(NULL, Image::getGlobalDetectorDistance, Image::setGlobalDetectorDistance, 1.0, 0.01);
    totalRefinementMap->setVerbose(true);
    
    totalRefinementMap->setEvaluationFunction(getter, object);
    totalRefinementMap->setCycles(30);
    
    totalRefinementMap->refine(GetterSetterStepSearch);
    
    if (pseudo)
    {
        powderManager->powderPattern("powder_after.csv");
    }*/
    
    logged << "New parameters (please enter into your input file): " << std::endl << std::endl;
    logged << "DETECTOR_DISTANCE " << Image::getGlobalDetectorDistance(NULL) << std::endl;
    logged << "BEAM_CENTRE " << Image::getGlobalBeamX(NULL) << " " << Image::getGlobalBeamY(NULL) << std::endl << std::endl;
    Logger::log(logged);
    
}