//
//  IndexManager.cpp
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "IndexManager.h"
#include "FileParser.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "Logger.h"
#include <algorithm>
#include <iomanip>
#include "parameters.h"
#include <fstream>
#include "SpotVector.h"
#include "IndexingSolution.h"
#include "FreeLattice.h"
#include "CSV.h"
#include "Reflection.h"
#include "Miller.h"
#include "misc.h"
#include "Detector.h"

IndexManager::IndexManager(std::vector<ImagePtr> newImages)
{
    images = newImages;
    scoreType = PseudoScoreTypeIntraPanel;
    proportionDistance = FileParser::getKey("DISTANCE_VS_ANGLE_FRACTION", 1.0);
    _canLockVectors = false;
    _axisWeighting = PseudoScoreWeightingAxisNone;
    _maxFrequency = -1;
    interPanelDistance = 0.07;
    intraPanelDistance = 0.07;
    
    unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());
    
    if (unitCell.size() < 6)
    {
        std::cout << "Please supply target unit cell in keyword UNIT_CELL." << std::endl;
        exit(1);
    }
    
    unitCellOnly = Matrix::matrixFromUnitCell(unitCell);
    
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    maxDistance = 0;
    lattice = UnitCellLattice::getMainLattice();
 
    sendLog(LogLevelDebug);
}


void IndexManager::indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset)
{
    std::ostringstream logged;
    
    while (true)
    {
        ImagePtr image = indexer->getNextImage();
        
        if (!image)
        {
            logged << "Finishing thread " << offset << std::endl;
            Logger::mainLogger->addStream(&logged); logged.str("");
            return;
        }
        
        logged << "Starting image " << image->getFilename() << " on thread " << offset << std::endl;
        Logger::mainLogger->addStream(&logged); logged.str("");
        
        image->findIndexingSolutions();
        
        std::vector<MtzPtr> mtzs = image->currentMtzs();

		image->dropImage();

        mtzSubset->reserve(mtzSubset->size() + mtzs.size());
        mtzSubset->insert(mtzSubset->begin(), mtzs.begin(), mtzs.end());
    }
}

double IndexManager::debugPseudoScore(void *object)
{
    Detector::getMaster()->drawSpecialImage();
    double value = pseudoScore(object);
    std::ostringstream logged;
    logged << "Value " << value << std::endl;
    Logger::log(logged);
    
    return value;
}

bool IndexManager::checkVector(SpotVectorPtr vec, bool permissive)
{
    bool bc = vec->usesBeamCentre();
    bool isInterPanel = (scoreType == PseudoScoreTypeInterPanel || scoreType == PseudoScoreTypeAllInterPanel);
    DetectorPtr activeDetector = getActiveDetector();
    
    if (!activeDetector)
    {
        activeDetector = Detector::getMaster();
    }
    
    if (scoreType == PseudoScoreTypeAngleConsistency && vec->isIntraPanelVector())
    {
        return true;
    }
    
    if (isInterPanel && vec->isIntraPanelVector() && !permissive)
    {
        return false;
    }
    
    if (scoreType == PseudoScoreTypeAllInterPanel && bc)
    {
        return true;
    }
    
    if (scoreType == PseudoScoreTypeAllInterPanel)
    {
        return true;
    }
    
    if (scoreType == PseudoScoreTypeBeamCentre && bc)
    {
        return true;
    }
    
    if (bc)
    {
        return false;
    }
    
    if (isInterPanel && permissive)
    {
        return false;
    }

    if (scoreType == PseudoScoreTypeIntraPanel)
    {
        if (!vec->isIntraPanelVector())
        {
            return false;
        }
        
        if (vec->isOnlyFromDetector(activeDetector) && vec->originalDistanceLessThan(intraPanelDistance))
        {
            return true;
        }
    }
    
    if (scoreType == PseudoScoreTypeInterPanel && vec->spansChildrenOfDetector(activeDetector)
        && vec->originalDistanceLessThan(interPanelDistance))
    {
        return true;
    }

    return false;
}

double IndexManager::pseudoScore(void *object)
{
    IndexManager *me = static_cast<IndexManager *>(object);
 
    if (me->_axisWeighting != PseudoScoreWeightingAxisNone)
    {
        double distanceScore = me->pseudoDistanceConsistency();
        return distanceScore;
    }
    
    double angleScore = (me->proportionDistance < 0.999) ? pseudoAngleScore(object) : 0;
    double distanceScore = (me->proportionDistance > 0.0001) ? pseudoDistanceScore(object) : 0;
    
    angleScore *= (1 - me->proportionDistance);
    distanceScore *= me->proportionDistance;
    
    return angleScore + distanceScore;
}

void IndexManager::processConsistencyVector(SpotVectorPtr vect, CSVPtr distCSV, bool lock, bool skipCheck)
{
    if (!skipCheck && !checkVector(vect))
    {
        return;
    }

    vect->setUpdate();
    double realDistance = vect->distance();
    vec firstSpotVec = vect->getFirstSpot()->estimatedVector();
    vec secondSpotVec = vect->getFirstSpot()->estimatedVector();
    
    add_vector_to_vector(&firstSpotVec, secondSpotVec);
    double axisValue =  (_axisWeighting == PseudoScoreWeightingAxisH) ? firstSpotVec.h : firstSpotVec.k;
    
    double distanceWeight = lattice->weightForDistance(realDistance);
    double axisWeight = 1;//fabs(axisValue);
    distanceWeight *= sqrt(axisWeight);
//    distanceWeight *= vect->getFirstSpot()->getParentImage()->getSpotVectorWeight();
    
    int column = axisValue > 0 ? 1 : 2;
    
    distCSV->addOneToFrequency(realDistance, column, distanceWeight);
    
    if (lock)
    {
        goodVectors.push_back(vect);
    }
}

bool IndexManager::processVector(SpotVectorPtr vec, double *score, double *count, bool lock, bool skipCheck)
{
    if (!skipCheck && !checkVector(vec))
    {
        return false;
    }
    
    vec->setUpdate();
    double realDistance = vec->distance();
    
    double value = lattice->weightForDistance(realDistance);
    double weight = 1;
    
    *score += weight * value;
    (*count) += weight;
    
    if (lock)
    {
        goodVectors.push_back(vec);
    }
    
    return true;
}

double IndexManager::pseudoDistanceConsistency()
{
    double score = 0;
    double count = 0;
    
    bool locked = (goodVectors.size() > 0);
    CSVPtr distCSV = CSVPtr(new CSV(0));
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
    double maxDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    distCSV->setupHistogram(0, maxDistance, step, "Distance", 4, "half1", "half2", "convol1", "convol2");
    
    if (!locked)
    {
        for (int i = 0; i < images.size(); i++)
        {
            for (int j = 0; j < images[i]->spotVectorCount(); j++)
            {
                SpotVectorPtr vec = images[i]->spotVector(j);
                processConsistencyVector(vec, distCSV, _canLockVectors, false);
            }
        }
    }
    else
    {
        for (int i = 0; i < goodVectors.size(); i++)
        {
            SpotVectorPtr vec = goodVectors[i];
            processConsistencyVector(vec, distCSV, false, true);
        }
    }
    
    for (int i = 0; i < distCSV->entryCount(); i++)
    {
        distCSV->convolutedPeaks("Distance", "half1", "convol1", 0, ConvolutionTypeUniform);
        distCSV->convolutedPeaks("Distance", "half2", "convol2", 0, ConvolutionTypeUniform);
        double value1 = distCSV->valueForEntry("convol1", i);
        double value2 = distCSV->valueForEntry("convol2", i);
        
        score -= value1 * value2;
        count++;
    }
    
    return score / pow(count, 2);
}

std::string targetString(PseudoScoreType type)
{
    switch (type) {
        case PseudoScoreTypeAllInterPanel:
            return "all_interpanel";
            break;
        case PseudoScoreTypeInterPanel:
            return "interpanel";
            break;
        case PseudoScoreTypeIntraPanel:
            return "intrapanel";
            break;
 
        default:
            return "";
            break;
    }
    
}

void IndexManager::plotGoodVectors()
{
    if (!goodVectors.size())
    {
        return;
    }
    
    ImagePtr special = Detector::getSpecialImage();
    std::string filename = getActiveDetector()->getTag() + "_" + targetString(getPseudoScoreType()) + "_vec.png";
    
    special->plotVectorsOnPNG(goodVectors, filename);
}

double IndexManager::pseudoDistanceScore(void *object, bool writeToCSV, std::string tag)
{
    IndexManager *me = static_cast<IndexManager *>(object);
    double score = 0;
    double count = 0;
    
    CSVPtr csv = CSVPtr(new CSV());
    
    if (writeToCSV)
    {
        double maxDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15) + DISTANCE_BUFFER;
        double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005) * 20;

        csv->setupHistogram(0, maxDistance, step, "distance", 2, "observed", "model");
        
        for (int i = 0; i < csv->entryCount(); i++)
        {
            double distance = csv->valueForEntry("distance", i);
            double weight = me->lattice->weightForDistance(distance);
            csv->addOneToFrequency(distance, "model", weight);
        }
    }
    
    bool locked = me->goodVectors.size();
    
    if (!locked)
    {
        for (int i = 0; i < me->images.size(); i++)
        {
            for (int j = 0; j < me->images[i]->spotVectorCount(); j++)
            {
                SpotVectorPtr vec = me->images[i]->spotVector(j);
                bool good = me->processVector(vec, &score, &count, me->_canLockVectors, false);
                
                if (writeToCSV && good)
                {
                    csv->addOneToFrequency(vec->distance(), "observed", 1);
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < me->goodVectors.size(); i++)
        {
            SpotVectorPtr vec = me->goodVectors[i];
            me->processVector(vec, &score, &count, false, true);
            
            if (writeToCSV)
            {
                csv->addOneToFrequency(vec->distance(), "observed", 1);
            }

        }
    }
    
    if (writeToCSV)
    {
        if (me->_maxFrequency < 0)
        {
            double min;
            csv->minMaxCol(1, &min, &me->_maxFrequency, true);
        }
        
        std::map<std::string, std::string> plotMap;
        plotMap["xHeader0"] = "distance";
        plotMap["yMax0"] = f_to_str(me->_maxFrequency, -1);
        plotMap["yMin0"] = "0";
        plotMap["xMin0"] = "0";
        plotMap["xTitle0"] = "Reciprocal distance between vectors (Ang)";
        plotMap["yHeader0"] = "observed";
        plotMap["style0"] = "line";
        
        if (me->getPseudoScoreType() == PseudoScoreTypeInterPanel)
        {
            plotMap["xMax0"] = f_to_str(me->interPanelDistance);
        }
        else if (me->getPseudoScoreType() == PseudoScoreTypeIntraPanel)
        {
            plotMap["xMax0"] = f_to_str(me->intraPanelDistance);
        }
        
        plotMap["xHeader1"] = "distance";
        plotMap["yHeader1"] = "model";
        plotMap["style1"] = "line";
        plotMap["colour1"] = "blue";
        plotMap["xMin1"] = "0";
        plotMap["xMax1"] = plotMap["xMax0"];
        
        plotMap["filename"] = tag + me->getActiveDetector()->getTag() + "_" + targetString(me->getPseudoScoreType()) + "_" + i_to_str(me->getCycleNum());
        
        csv->plotPNG(plotMap);
        
        std::ostringstream logged;
       // logged << plotMap["filename"] << " used " << count << " vectors." << std::endl;
       // Logger::log(logged);
    }
    
    score /= count;
    
    return score;
}

void IndexManager::pseudoAngleCSV()
{
    double angleDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);
    
    if (angleDistance <= 0)
    {
        return;
    }
    
    angleCSV = CSVPtr(new CSV(0));
    angleCSV->setupHistogram(0, 90, 0.1, "Angle", 3, "all", "intra-angle", "inter-angle");
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr vec1 = images[i]->spotVector(j);
            
            if (!vec1->originalDistanceLessThan(angleDistance))
                continue;
            
            for (int k = 0; k < j; k++)
            {
                SpotVectorPtr vec2 = images[i]->spotVector(k);
                
                if (!vec2->originalDistanceLessThan(angleDistance))
                    continue;
                
                if (k == j)
                {
                    continue;
                }
                
                double angle = vec1->angleWithVector(vec2);
                
                angle *= 180 / M_PI;
                angle = (angle > 90) ? 180 - angle : angle;
                
                double weight = lattice->weightForDistance(vec1->distance());
                weight *= lattice->weightForDistance(vec2->distance());

                angleCSV->addOneToFrequency(angle, "all", weight);
                
                if (vec1->isIntraPanelVector() && vec2->isIntraPanelVector())
                {
                    angleCSV->addOneToFrequency(angle, "intra-angle", weight);
                }
                else
                {
                    angleCSV->addOneToFrequency(angle, "inter-angle", weight);
                }
            }
        }
    }
}

PseudoScoreType IndexManager::checkVectors(SpotVectorPtr vec1, SpotVectorPtr vec2)
{
    DetectorPtr activeDetector = getActiveDetector();
    DetectorPtr brother = activeDetector->getChild(0);
    DetectorPtr sister = activeDetector->getChild(1);
    bool isInterPanel = scoreType == PseudoScoreTypeAllInterPanel || scoreType == PseudoScoreTypeInterPanel;
    
    if (!checkVector(vec1, true))
    {
        return PseudoScoreTypeInvalid;
    }

    if (!checkVector(vec2, true))
    {
        return PseudoScoreTypeInvalid;
    }
    
    if (scoreType == PseudoScoreTypeAngleConsistency)
    {
        return scoreType;
    }
    
    if (!isInterPanel)
    {
        if (!(vec1->isIntraPanelVector() && vec2->isIntraPanelVector()))
        {
            return PseudoScoreTypeInvalid;
        }
        
        if (vec1->getFirstSpot()->getDetector() != vec2->getFirstSpot()->getDetector())
        {
            return PseudoScoreTypeInvalid;
        }
        
        return PseudoScoreTypeIntraPanel;
    }
    else
    {
        bool vec1Okay = false;
        bool vec2Okay = false;
        
        if (!vec1->isIntraPanelVector())
        {
            vec1Okay = true;
        }

        if (!vec2->isIntraPanelVector())
        {
            vec2Okay = true;
        }
        
        if (vec1Okay || vec2Okay)
        {
            return PseudoScoreTypeInterPanel;
        }
        
        DetectorPtr dec1 = vec1->getFirstSpot()->getDetector();
        DetectorPtr dec2 = vec2->getFirstSpot()->getDetector();
        
        vec1Okay = true;
        vec2Okay = true;
        
        if (!brother || (brother->isAncestorOf(dec1) && brother->isAncestorOf(dec2)))
        {
            vec1Okay = false;
        }
        
        if (!sister || (sister->isAncestorOf(dec1) && sister->isAncestorOf(dec2)))
        {
            vec2Okay = false;
        }
        
        if (dec1 != dec2 && scoreType == PseudoScoreTypeAllInterPanel)
        {
            return PseudoScoreTypeAllInterPanel;
        }
        
        if (vec1Okay || vec2Okay)
        {
            return PseudoScoreTypeInterPanel;
        }
    }
    
    return PseudoScoreTypeInvalid;
}

double IndexManager::pseudoAngleScore(void *object)
{
    IndexManager *me = static_cast<IndexManager *>(object);
    DetectorPtr activeDetector = me->getActiveDetector();
    
    double angleDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);
    double score = 0;
    double count = 0;
    int panelCount = 0;
    int staticCount = 0;
    
    bool firstCSV = (me->angleConsistencyCSV == CSVPtr());
    
    if (me->scoreType == PseudoScoreTypeAngleConsistency)
    {
        if (firstCSV)
        {
            me->angleConsistencyCSV = CSVPtr(new CSV(0));
            
            me->angleConsistencyCSV->setupHistogram(0, 90, 0.2, "Angle", 4, "Cosine", "Static", "Convol", "Panel");
            
            for (int i = 0; i < me->angleConsistencyCSV->entryCount(); i++)
            {
                double angle = me->angleConsistencyCSV->valueForEntry("Angle", i);
                double cosine = cos(angle * M_PI / 180);
                me->angleConsistencyCSV->setValueForEntry(i, "Cosine", cosine);
            }
        }
        else
        {
            if (me->scoreType == PseudoScoreTypeAngleConsistency)
            {
                me->angleConsistencyCSV->resetColumn("Panel", 0);
            }
        }
    }
    
    for (int i = 0; i < me->images.size(); i++)
    {
        for (int j = 0; j < me->images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr vec1 = me->images[i]->spotVector(j);
            
            if (!vec1->originalDistanceLessThan(angleDistance))
                continue;
            
            vec1->calculateDistance();
            double distanceWeight1 = me->lattice->weightForDistance(vec1->distance());
            
            if (distanceWeight1 <= 0)
                continue;
            
            if (vec1->usesBeamCentre())
                continue;
            
            if (me->scoreType == PseudoScoreTypeAngleConsistency)
            {
                bool isPanel = vec1->isOnlyFromDetector(activeDetector);
                int *counter = isPanel ? &panelCount : &staticCount;
                (*counter)++;
                
                std::string chosenVecType = isPanel ? "Panel" : "Static";
                double cosine = vec1->cosineWithVertical();
                
                if (firstCSV || (!firstCSV && isPanel))
                {
                    me->angleConsistencyCSV->addOneToFrequency(cosine, chosenVecType, distanceWeight1, "Cosine");
                }
                
                continue;
            }
            
            for (int k = 0; k < j; k++)
            {
                SpotVectorPtr vec2 = me->images[i]->spotVector(k);
                vec2->calculateDistance();
                
                if (vec2->usesBeamCentre())
                    continue;

                if (!vec2->originalDistanceLessThan(angleDistance))
                    continue;

                if (k == j)
                {
                    continue;
                }
                
                if (me->checkVectors(vec1, vec2) != me->scoreType)
                {
                    continue;
                }
                
                double vec2Distance = vec2->distance();
                
                double distanceWeight2 = me->lattice->weightForDistance(vec2Distance);
                
                double cosine = vec1->cosineWithVector(vec2);
                cosine = fabs(cosine);
                
                double angleWeight = me->lattice->weightForCosine(cosine);
                
                if (distanceWeight2 <= 0)
                    continue;
                
                score += distanceWeight1 * distanceWeight2 * angleWeight / 10000000;
                count += distanceWeight1 * distanceWeight2 / 10000000;
            }
        }
    }
    
    if (me->scoreType == PseudoScoreTypeAngleConsistency)
    {
        if (firstCSV)
        {
            me->angleConsistencyCSV->convolutedPeaks("Angle", "Static", "Convol", 0.5);
        }
        
        for (int i = 0; i < me->angleConsistencyCSV->entryCount(); i++)
        {
            double panelValue = me->angleConsistencyCSV->valueForEntry("Panel", i);
            double staticValue = me->angleConsistencyCSV->valueForEntry("Convol", i);
            
            score += panelValue * staticValue;
            count++;
        }
    }
    
    score /= count;
    return -sqrt(score);
}

void IndexManager::powderPattern(std::string csvName, bool force)
{
    bool alwaysFilterSpots = FileParser::getKey("ALWAYS_FILTER_SPOTS", false);
    
    std::ostringstream pdbLog;


    for (int i = 0; i < images.size(); i++)
    {
        if (force)
        {
            images[i]->compileDistancesFromSpots(maxDistance, smallestDistance, alwaysFilterSpots);
        }
        
        for (int j = 0; j < images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr spotVec = images[i]->spotVector(j);
            spotVec->calculateDistance();
        }
    }

    
    PowderHistogram allFrequencies = generatePowderHistogram();
    PowderHistogram intraFrequencies = generatePowderHistogram(1);
    PowderHistogram interFrequencies = generatePowderHistogram(0);
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr spotVec = images[i]->spotVector(j);
            
            vec spotDiff = copy_vector(spotVec->getVector());
            spotDiff.h *= 106 * 20;
            spotDiff.k *= 106 * 20;
            spotDiff.l *= 106 * 20;

            if (i == 0)
            {
                pdbLog << "HETATM";
                pdbLog << std::fixed;
                pdbLog << std::setw(5) << j << "                   ";
                pdbLog << std::setprecision(2) << std::setw(8)  << spotDiff.h;
                pdbLog << std::setprecision(2) << std::setw(8) << spotDiff.k;
                pdbLog << std::setprecision(2) << std::setw(8) << spotDiff.l;
                pdbLog << "                       O" << std::endl;
            }
        }
    }
    
    if (images.size() == 0)
    {
        logged << "No images specified." << std::endl;
        sendLog();
        return;
    }
    
    pdbLog << "HETATM";
    pdbLog << std::fixed;
    pdbLog << std::setw(5) << images[0]->spotVectorCount() + 1 << "                   ";
    pdbLog << std::setw(8) << std::setprecision(2) << 0;
    pdbLog << std::setw(8) << std::setprecision(2) << 0;
    pdbLog << std::setw(8) << std::setprecision(2) << 0;
    pdbLog << "                       N" << std::endl;
    
    std::ofstream pdbLog2;
    pdbLog2.open("projection.pdb");
    
    for (int j = 0; j < images[0]->spotCount(); j++)
    {
        vec estimatedVector = images[0]->spot(j)->estimatedVector();
        
        pdbLog2 << "HETATM";
        pdbLog2 << std::fixed;
        pdbLog2 << std::setw(5) << j << "                   ";
        pdbLog2 << std::setw(8) << std::setprecision(2) << estimatedVector.h * 106;
        pdbLog2 << std::setw(8) << std::setprecision(2) << estimatedVector.k * 106;
        pdbLog2 << std::setw(8) << std::setprecision(2) << estimatedVector.l * 106;
        pdbLog2 << "                       O" << std::endl;
    }
    
    pdbLog2 << "HETATM";
    pdbLog2 << std::fixed;
    pdbLog2 << std::setw(5) << images[0]->spotCount() + 1 << "                   ";
    pdbLog2 << std::setw(8) << std::setprecision(2) << 0;
    pdbLog2 << std::setw(8) << std::setprecision(2) << 0;
    pdbLog2 << std::setw(8) << std::setprecision(2) << 0;
    pdbLog2 << "                       N" << std::endl;
    
    pdbLog2.close();
    
    logged << "******* DISTANCE FREQUENCY *******" << std::endl;
    
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);

    CSV powder(5, "Distance", "Frequency", "Perfect frequency", "Intra-panel", "Inter-panel");
    
    for (PowderHistogram::iterator it = allFrequencies.begin(); it != allFrequencies.end(); it++)
    {
        double distance = it->first * step;
        double freq = it->second.first;
        double perfect = it->second.second;
        double intrapanel = intraFrequencies[it->first].first;
        double interpanel = interFrequencies[it->first].first;
        
        powder.addEntry(0, distance, freq, perfect, intrapanel, interpanel);
    }
    
    if (angleCSV)
    {
        angleCSV->writeToFile("angle_" + csvName);
    }
    powder.writeToFile(csvName);
    
    std::map<std::string, std::string> plotMap;
    plotMap["xHeader0"] = "Distance";
    plotMap["xTitle0"] = "Reciprocal distance between vectors (Ang)";
    plotMap["yHeader0"] = "Intra-panel";
    plotMap["style0"] = "line";
    plotMap["round0"] = "true";
    plotMap["xHeader1"] = "Distance";
    plotMap["xTitle1"] = "Reciprocal distance between vectors (Ang)";
    plotMap["yHeader1"] = "Inter-panel";
    plotMap["style1"] = "line";
    plotMap["colour1"] = "blue";
    plotMap["filename"] = csvName;


    powder.plotPNG(plotMap);
    logged << "Written to " << csvName << std::endl;
    sendLog();
    
    if (force)
    {
        powder.plotColumns(0, 1);
    }
    
    lattice->powderPattern(false, "unitCellLatticePowder.csv");
    
  //  powderLog.close();
    
    /// angles
    
    std::vector<double> probeDistances = FileParser::getKey("PROBE_DISTANCES", std::vector<double>());
    
    if (probeDistances.size() < 2)
        return;
    
    double distance1 = probeDistances[0]; double distance2 = probeDistances[1];
    
    double angleStep = FileParser::getKey("POWDER_PATTERN_STEP_ANGLE", 2.0) * M_PI / 180;
    double distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 500.);
    std::map<int, int> angleHistogram;
    std::map<int, int> perfectAngleHistogram;
    
    for (double i = 0; i < M_PI; i += angleStep)
    {
        int angleCategory = i / angleStep;
        angleHistogram[angleCategory] = 0;
        perfectAngleHistogram[angleCategory] = 0;
    }
    
    for (int i = 0; i < images.size(); i++)
    {
        std::vector<double> newAngles = images[i]->anglesBetweenVectorDistances(distance1, distance2, distanceTolerance);
        
        for (int j = 0; j < newAngles.size(); j++)
        {
            int angleCategory = newAngles[j] / angleStep;
            angleHistogram[angleCategory]++;
        }
    }
    
    std::vector<VectorDistance> perfectProbe1, perfectProbe2;
    
    for (int i = 0; i < vectorDistances.size(); i++)
    {
        double distance = vectorDistances[i].second;
        bool chosen = false;
        
        if (1 / fabs(probeDistances[0] - distance) > distanceTolerance)
        {
            perfectProbe1.push_back(vectorDistances[i]);
            chosen = true;
        }

        if (1 / fabs(probeDistances[1] - distance) > distanceTolerance)
        {
            perfectProbe2.push_back(vectorDistances[i]);
            chosen = true;
        }
        
        if (chosen)
        {
            logged << "Matching vector: " << vectorDistances[i].first.h << "\t" << vectorDistances[i].first.k << "\t" << vectorDistances[i].first.l << std::endl;
            sendLog();
        }

    }
    
    logged << "Comparing " << perfectProbe1.size() << " against " << perfectProbe2.size() << std::endl;
    sendLog();
    
    for (int i = 0; i < perfectProbe1.size(); i++)
    {
        for (int j = 0; j < perfectProbe2.size(); j++)
        {
            vec probeTransformed1 = copy_vector(perfectProbe1[i].first);
            vec probeTransformed2 = copy_vector(perfectProbe2[j].first);
            
            if (perfectProbe1[i].first.h == perfectProbe2[i].first.h &&
                perfectProbe1[i].first.k == perfectProbe2[i].first.k &&
                perfectProbe1[i].first.l == perfectProbe2[i].first.l)
            {
         //       continue;
            }
            
            unitCellMatrix->multiplyVector(&probeTransformed1);
            unitCellMatrix->multiplyVector(&probeTransformed2);
            
            double angle = angleBetweenVectors(probeTransformed1, probeTransformed2);
            double mirror = M_PI - angle;
            
            int angleCategory = angle / angleStep;
            perfectAngleHistogram[angleCategory]++;
            
            angleCategory = mirror / angleStep;
            perfectAngleHistogram[angleCategory]++;
        }
    }
    
    
    CSV csv(3, "Angle", "Frequency", "Perfect frequency");
    
    for (std::map<int, int>::iterator it = angleHistogram.begin(); it != angleHistogram.end(); it++)
    {
        double angle = it->first * angleStep * 180 / M_PI;
        double freq = it->second;
        double perfectFreq = 0;
        
        if (perfectAngleHistogram.count(it->first))
            perfectFreq = perfectAngleHistogram[it->first];
        
        csv.addEntry(0, angle, freq, perfectFreq);
        
 //       angleLog << angle << "," << freq << "," << perfectFreq << "," << std::endl;
    }
    
    csv.writeToFile("angle.csv");
    csv.plotColumns(0, 1);
    
   // angleLog.close();
}

void IndexManager::refineUnitCell()
{
    bool alwaysFilterSpots = FileParser::getKey("ALWAYS_FILTER_SPOTS", false);
    
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->compileDistancesFromSpots(maxDistance, smallestDistance, alwaysFilterSpots);
    }
    
    PowderHistogram frequencies = generatePowderHistogram();
    
    lattice->refineUnitCell(frequencies);
}

ImagePtr IndexManager::getNextImage()
{
    std::lock_guard<std::mutex> lock(indexMutex);
    
    nextImage++;
    
    if (nextImage >= images.size())
    {
        return ImagePtr();
    }
    else
    {
        return images[nextImage];
    }
}

void IndexManager::index()
{
    int maxThreads = FileParser::getMaxThreads();
    IndexingSolution::setupStandardVectors();
    
    boost::thread_group threads;
    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    nextImage = -1;
    
    for (int num = 0; num < 1; num++)
    {
        time_t startcputime;
        time(&startcputime);

        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(indexThread, this, &managerSubsets[i], i);
            threads.add_thread(thr);
        }

        threads.join_all();
    
        time_t endcputime;
        time(&endcputime);
        
        nextImage = -1;
    }
    
    
    int total = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }
    
    mtzs.reserve(total);
    int lastPos = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        mtzs.insert(mtzs.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }
}

PowderHistogram IndexManager::generatePowderHistogram(int intraPanel, int perfectPadding)
{
    PowderHistogram frequencies;
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
    
    double maxCatDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    int maxCategory = maxCatDistance / step;
    
    for (int i = 0; i < maxCategory; i++)
    {
        frequencies[i] = std::make_pair(0, 0);
    }
    
    for (int i = 0; i < lattice->standardVectorCount(); i++)
    {
        double distance = lattice->standardVector(i)->distance();
        int categoryNum = distance / step;
        
        if (categoryNum > maxCategory)
            continue;
        
        for (int i = categoryNum - perfectPadding; i <= categoryNum + perfectPadding; i++)
        {
            if (i < 0 || i > maxCategory)
                continue;
            
            frequencies[i].second++;
        }
    }
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr spotVec = images[i]->spotVector(j);
            
            double distance = spotVec->distance();
            int categoryNum = distance / step;
            
            if (categoryNum > maxCategory || categoryNum < 0)
                continue;

            bool isIntraPanel = spotVec->isIntraPanelVector();
            
            if (intraPanel == -1)
            {
                frequencies[categoryNum].first++;
            }
            else if (intraPanel == 0 && !isIntraPanel)
            {
                frequencies[categoryNum].first++;
            }
            else if (intraPanel == 1 && isIntraPanel)
            {
                frequencies[categoryNum].first++;
            }
        }
    }
    
    return frequencies;
}

void IndexManager::updateAllSpots()
{
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->updateAllSpots();
    }
}

void IndexManager::indexingParameterAnalysis()
{
    std::vector<double> goodNetworkCount, badNetworkCount;
    std::vector<double> goodDistanceTrusts, badDistanceTrusts;
    std::vector<double> goodAngleTrusts, badAngleTrusts;
    std::vector<double> goodDistances, badDistances;
    
    logged << "Indexing parameter analysis." << std::endl;
    sendLog();
    
    for (int i = 0; i < images.size(); i++)
    {
        logged << "Processing image " << images[i]->getFilename() << std::endl;
        sendLog();
        
        for (int good = 0; good < 2; good++)
        {
            std::vector<double> *networkCount = good ? &goodNetworkCount : &badNetworkCount;
            std::vector<double> *distanceTrusts = good ? &goodDistanceTrusts : &badDistanceTrusts;
            std::vector<double> *angles = good ? &goodAngleTrusts : &badAngleTrusts;
            std::vector<double> *distances = good ? &goodDistances : &badDistances;
            
            for (int j = 0; j < images[i]->goodOrBadSolutionCount(good); j++)
            {
                IndexingSolutionPtr solution = images[i]->getGoodOrBadSolution(j, good);
                
                networkCount->push_back(solution->spotVectorCount());
                
                std::vector<double> additionalTrusts = solution->totalDistanceTrusts();
                distanceTrusts->reserve(distanceTrusts->size() + additionalTrusts.size());
                distanceTrusts->insert(distanceTrusts->begin(), additionalTrusts.begin(), additionalTrusts.end());
                
                std::vector<double> additionalDistances = solution->totalDistances();
                distances->reserve(distances->size() + additionalDistances.size());
                distances->insert(distances->begin(), additionalDistances.begin(), additionalDistances.end());
                
                std::vector<double> additionalAngles = solution->totalAngles();
                angles->reserve(angles->size() + additionalAngles.size());
                angles->insert(angles->begin(), additionalAngles.begin(), additionalAngles.end());
            }
        }
    }
    /*
    std::vector<double> goodNetworkCount, badNetworkCount;
    std::vector<double> goodDistanceTrusts, badDistanceTrusts;
    std::vector<double> goodAngleTrusts, badAngleTrusts;
    std::vector<double> goodDistances, badDistances;
    */
    
    std::map<double, int> goodNetworkCountHistogram = histogram(goodNetworkCount, 2);
    std::map<double, int> badNetworkCountHistogram = histogram(badNetworkCount, 2);
    histogramCSV("networkCountAnalysis.csv", goodNetworkCountHistogram, badNetworkCountHistogram);

    std::map<double, int> goodDistanceTrustsHistogram = histogram(goodDistanceTrusts, 100);
    std::map<double, int> badDistanceTrustsHistogram = histogram(badDistanceTrusts, 100);
    histogramCSV("distanceTrustAnalysis.csv", goodDistanceTrustsHistogram, badDistanceTrustsHistogram);

    double powderPatternStep = FileParser::getKey("POWDER_PATTERN_STEP", 0.001);
    std::map<double, int> goodDistancesHistogram = histogram(goodDistances, powderPatternStep);
    std::map<double, int> badDistancesHistogram = histogram(badDistances, powderPatternStep);
    histogramCSV("distanceAnalysis.csv", goodDistancesHistogram, badDistancesHistogram);

    std::map<double, int> goodAnglesHistogram = histogram(goodAngleTrusts, 0.05);
    std::map<double, int> badAnglesHistogram = histogram(badAngleTrusts, 0.05);
    histogramCSV("angleTrustAnalysis.csv", goodAnglesHistogram, badAnglesHistogram);
}

std::vector<MtzPtr> IndexManager::consolidateOrientations(ImagePtr image1,
ImagePtr image2, int *oneHand, int *otherHand, int *both)
{
    std::vector<MtzPtr> refiners;
    
    for (int i = 0; i < image1->mtzCount(); i++)
    {
        bool found = false;
        MtzPtr refiner1 = image1->mtz(i);
        
        for (int j = 0; j < image2->mtzCount(); j++)
        {
            MtzPtr refiner2 = image2->mtz(j);
            
            MatrixPtr mat1 = refiner1->getMatrix();
            MatrixPtr mat2 = refiner2->getMatrix();
            
            if (IndexingSolution::matrixSimilarToMatrix(mat1, mat2))
            {
                found = true;
                refiners.push_back(refiner1);
                (*both)++;

                // ignore refiner2 as it is duplicate
            }
        }
        
        if (found == false)
        {
            refiners.push_back(refiner1);
            (*oneHand)++;
        }
    }
    
    for (int i = 0; i < image2->mtzCount(); i++)
    {
        MtzPtr refiner1 = image2->mtz(i);
        bool found = false;
        
        for (int j = 0; j < image1->mtzCount(); j++)
        {
            MtzPtr refiner2 = image1->mtz(j);
            
            MatrixPtr mat1 = refiner1->getMatrix();
            MatrixPtr mat2 = refiner2->getMatrix();
            
            if (IndexingSolution::matrixSimilarToMatrix(mat1, mat2))
            {
                found = true;
                // ignore both as it would have been done in the previous loop
                break;
            }
        }
        
        if (!found)
        {
            refiners.push_back(refiner1);
            (*otherHand)++;
        }
    }
    
    return refiners;
}

void IndexManager::combineLists()
{
    IndexingSolution::setupStandardVectors();
    
    std::vector<ImagePtr> mergedImages;
    logged << "Combining lists from " << images.size() << " images and " << mergeImages.size() << " images." << std::endl;
    
    int oneHand = 0;
    int otherHand = 0;
    int both = 0;
    
    // first we iterate through array 1
    for (int i = 0; i < images.size(); i++)
    {
        ImagePtr image1 = images[i];
        bool found = false;
        
        for (int j = 0; j < mergeImages.size(); j++)
        {
            ImagePtr image2 = mergeImages[j];
            
            // Are image1 and image2 the same
            
            if (image1->getFilename() == image2->getFilename())
            {
                found = true;
                // we have a duplicate! check individual orientations
                
                std::vector<MtzPtr> newRefiners = consolidateOrientations(image1, image2, &oneHand, &otherHand, &both);
                image1->clearMtzs();
                image1->addMtzs(newRefiners);
                mergedImages.push_back(image1);
                
                break;
            }
        }
        
        if (!found)
        {
            mergedImages.push_back(image1);
            oneHand += image1->mtzCount();
        }
    }
    
    for (int i = 0; i < mergeImages.size(); i++)
    {
        ImagePtr image1 = mergeImages[i];
        bool found = false;
        
        for (int j = 0; j < images.size(); j++)
        {
            ImagePtr image2 = images[j];
            // Are image1 and image2 the same
            
            if (image1->getFilename() == image2->getFilename())
            {
                // we do not consider this image because it was covered in the last loop
                found = true;
                break;
            }
        }
        
        if (!found)
        {
            mergedImages.push_back(image1);
            otherHand += image1->mtzCount();
        }
    }
    
    logged << "Unique to list 1: " << oneHand << std::endl;
    logged << "Unique to list 2: " << otherHand << std::endl;
    logged << "Present in both lists: " << both << std::endl;

    sendLog();
    
    images = mergedImages;
    mergeImages.clear();
    std::vector<ImagePtr>().swap(mergeImages);
}