//
//  GeometryRefiner.cpp
//  cppxfel
//
//  Created by Helen Ginn on 28/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "GeometryRefiner.h"
#include "GeometryParser.h"
#include "Detector.h"
#include "RefinementStepSearch.h"
#include "RefinementGridSearch.h"
#include "NelderMead.h"
#include "IndexManager.h"
#include "misc.h"

RefinementStrategyPtr GeometryRefiner::makeRefiner(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr strategy = RefinementStrategy::userChosenStrategy();
    
    strategy->setVerbose(true);
    strategy->setCycles(30);
    strategy->setJobName("Detector " + detector->getTag());

    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    indexManagers.push_back(aManager);
    aManager->setActiveDetector(detector);
    
    switch (type)
    {
        case GeometryScoreTypeMiller:
            strategy->setEvaluationFunction(Detector::millerScoreWrapper, &*detector);
            break;
        case GeometryScoreTypeMillerStdev:
            strategy->setEvaluationFunction(Detector::millerStdevScoreWrapper, &*detector);
            break;
        case GeometryScoreTypeInterpanel:
            aManager->setPseudoScoreType(PseudoScoreTypeInterPanel);
            strategy->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
            break;
        case GeometryScoreTypeIntrapanel:
            aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
            strategy->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
            break;
    }
    
    return strategy;
}

GeometryRefiner::GeometryRefiner()
{
    refinementEvent = 0;
    cycleNum = 0;
    lastIntraScore = 0;
    lastInterScore = 0;
}

void GeometryRefiner::refineGeometry()
{
    Detector::setNoisy(true);
    
    reportProgress();
    std::vector<DetectorPtr> detectors;
    detectors.push_back(Detector::getMaster());
    
    for (int i = 0; i < 2; i++)
    {
        cycleNum++;
        geometryCycleForDetector(detectors);
    }
}

void GeometryRefiner::reportProgress()
{
    std::string filename = "special_image_" + i_to_str(refinementEvent) + ".png";
    Detector::drawSpecialImage(filename);
    
    double intraScore = IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeInterPanel);
    double interScore = IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    
    double intraIncrease = 100;
    double interIncrease = 100;
    double mScore = Detector::getMaster()->millerScore(true, false);
    double sScore = Detector::getMaster()->millerScore(false, true);
    
    interIncrease = 100 * (interScore - lastInterScore) / interScore;
    intraIncrease = 100 * (intraScore - lastIntraScore) / intraScore;
    
    lastInterScore = interScore;
    lastIntraScore = intraScore;
    
    if (refinementEvent > 0 && fabs(interIncrease) < 0.00001 && fabs(intraIncrease)  < 0.00001)
    {
        return;
    }
    
    logged << "N: Progress score (event " << refinementEvent << ", intra-panel): " << intraScore
    << " (" << (intraIncrease > 0 ? "+" : "") << intraIncrease << "% from last round) " << std::endl;
    logged << "N: Progress score (event " << refinementEvent << ", inter-panel): " << interScore
    << " (" << (interIncrease > 0 ? "+" : "") << interIncrease << "% from last round)" << std::endl;
    logged << "N: Progress score (event " << refinementEvent << ", miller mean): " << mScore << std::endl;
    logged << "N: Progress score (event " << refinementEvent << ", miller stdev): " << sScore << std::endl;
    sendLog();

    manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
    GeometryParser geomParser = GeometryParser("whatever", GeometryFormatCppxfel);
    geomParser.writeToFile("new_" + i_to_str(refinementEvent) + ".cppxfel_geom");
    refinementEvent++;
    
    sendLog();
}

void GeometryRefiner::geometryCycleForDetector(std::vector<DetectorPtr> detectors)
{
    std::ostringstream detectorList;
    
    for (int j = 0; j < detectors.size(); j++)
    {
        detectorList << detectors[j]->getTag() << " ";
    }

    
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining detector" << (detectors.size() == 1 ? "" : "s") << ": " << detectorList.str() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();

    if (detectors[0]->isLUCA())
    {
        refineMasterDetector();
    }
    else
    {
        for (int j = 0; j < 4; j++)
        {
            int maxThreads = FileParser::getMaxThreads();
            boost::thread_group threads;
            
            for (int i = 0; i < maxThreads; i++)
            {
                boost::thread *thr = new boost::thread(refineDetectorWrapper, this, detectors, i, j);
                threads.add_thread(thr);
            }
            
            threads.join_all();
            reportProgress();
        }
    }
    
    std::vector<DetectorPtr> nextDetectors;
    
    for (int i = 0; i < detectors.size(); i++)
    {
        for (int j = 0; j < detectors[i]->childrenCount(); j++)
        {
            nextDetectors.push_back(detectors[i]->getChild(j));
        }
    }
    
    if (nextDetectors.size())
    {
        geometryCycleForDetector(nextDetectors);
    }
}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset, int strategy)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < detectors.size(); i += maxThreads)
    {
        me->refineDetectorStrategy(detectors[i], strategy);
    }
}

void GeometryRefiner::refineMidPointXY(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr refinementMap = makeRefiner(detector, type);
    
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointX,
                                Detector::setArrangedMidPointX, 1.0, 0.01, "midPointX");
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointY,
                                Detector::setArrangedMidPointY, 1.0, 0.01, "midPointY");

    refinementMap->refine();
}

void GeometryRefiner::refineMidPointZ(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr refinementMap = makeRefiner(detector, type);
    
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                                Detector::setArrangedMidPointZ, 1.0, 0.01, "midPointZ");
    
    refinementMap->refine();
}

void GeometryRefiner::refineTiltXY(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr refinementMap = makeRefiner(detector, type);
    
    refinementMap->addParameter(&*detector, Detector::getAlpha, Detector::setAlpha, 0.001, 0.00001, "alpha");
    refinementMap->addParameter(&*detector, Detector::getBeta, Detector::setBeta, 0.001, 0.00001, "beta");
    
    refinementMap->refine();
}

void GeometryRefiner::refineTiltZ(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr refinementMap = makeRefiner(detector, type);
    
    refinementMap->addParameter(&*detector, Detector::getGamma, Detector::setGamma, 0.001, 0.00001, "gamma");
    
    refinementMap->refine();
}

void GeometryRefiner::refineVarious(DetectorPtr detector, GeometryScoreType type)
{
    RefinementGridSearchPtr gridMap = RefinementGridSearchPtr(new RefinementGridSearch());
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    indexManagers.push_back(aManager);
    aManager->setActiveDetector(detector);

    double origAlpha = Detector::getAlpha(&*detector);
    double origMidPointX = Detector::getArrangedMidPointX(&*detector);
    gridMap->setJobName("intrapanel_tiltX_posX");
    
    gridMap->addParameter(&*detector, Detector::getAlpha, Detector::setAlpha, 0.002, 0.00001);
    gridMap->addParameter(&*detector, Detector::getArrangedMidPointX, Detector::setArrangedMidPointX, 1.0, 0.1);
    gridMap->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    gridMap->refine();
    
    Detector::setAlpha(&*detector, origAlpha);
    Detector::setArrangedMidPointX(&*detector, origMidPointX);
    gridMap->setEvaluationFunction(Detector::millerScoreWrapper, &*Detector::getMaster());
    gridMap->refine();
    
}

void GeometryRefiner::refineMasterDetector()
{
    DetectorPtr detector = Detector::getMaster();

    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining master only: " << detector->getTag() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();
    
    refineVarious(detector, GeometryScoreTypeIntrapanel);
    reportProgress();
    refineMidPointXY(detector, GeometryScoreTypeMiller);
    reportProgress();
    refineMidPointZ(detector, GeometryScoreTypeIntrapanel);
    reportProgress();
    refineTiltXY(detector, GeometryScoreTypeIntrapanel);
    reportProgress();
    refineMidPointXY(detector, GeometryScoreTypeMiller);
    reportProgress();
}

void GeometryRefiner::refineDetectorStrategy(DetectorPtr detector, int strategy)
{
    if (strategy == 0)
    {
        if (!detector->hasChildren() && cycleNum == 1 && Detector::isCSPAD())
        {
       //     return;
        }
        
        refineMidPointXY(detector, GeometryScoreTypeMiller);
        detector->millerScore(true, false);
    }
    else if (strategy == 1)
    {
        if (!detector->hasChildren() && Detector::isCSPAD())
        {
            return;
        }
        
        refineVarious(detector, GeometryScoreTypeIntrapanel);
        
        refineMidPointZ(detector, GeometryScoreTypeIntrapanel);
        detector->millerScore(true, false);
    }
    else if (strategy == 2)
    {
        if (!detector->hasChildren() && Detector::isCSPAD())
        {
            return;
        }
        
        refineTiltXY(detector, GeometryScoreTypeIntrapanel);
        detector->millerScore(true, false);
    }
    else if (strategy == 3 && cycleNum == 1)
    {
        if (!detector->hasChildren() && Detector::isCSPAD())
        {
            return;
        }
        
    //    refineTiltZ(detector, GeometryScoreTypeMillerStdev);
    //    detector->millerScore(true, false);
    }
    
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));
    manager->powderPattern("geom_refinement_event_0_start.csv", false);
    
}