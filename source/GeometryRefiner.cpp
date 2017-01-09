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
#include "NelderMead.h"
#include "IndexManager.h"
#include "misc.h"

RefinementStrategyPtr GeometryRefiner::makeRefiner(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr strategy = RefinementStrategy::userChosenStrategy();
    
    strategy->setVerbose(true);
    strategy->setCycles(30);
    strategy->setJobName("Detector" + detector->getTag());

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
    reportProgress();
    
    for (int i = 0; i < 2; i++)
    {
        cycleNum++;
        refineMasterDetector();
    }

    for (int i = 0; i < 1; i++)
    {
        refineGeometryCycle();
    }
    
    GeometryParser geomParser = GeometryParser("test.cppxfel_geom", GeometryFormatCppxfel);
    geomParser.writeToFile("new.cppxfel_geom");
}

void GeometryRefiner::reportProgress()
{
    double intraScore = IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeInterPanel);
    double interScore = IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    
    double intraIncrease = 100;
    double interIncrease = 100;
    double mScore = Detector::getMaster()->millerScore(true);
    
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
    logged << "N: Progress score (event " << refinementEvent << ", miller): " << mScore << std::endl;
    sendLog();

    sendLog();
}

void GeometryRefiner::refineGeometryCycle()
{
    std::vector<DetectorPtr> sameLevelDetectors;
    std::vector<DetectorPtr> nextLevelDetectors;
    
    DetectorPtr master = Detector::getMaster();
    nextLevelDetectors.push_back(master);
    
    reportProgress();
    cycleNum++;
    refinementEvent++;
    std::ostringstream detectorList;
    
    for (int j = 0; j < master->childrenCount(); j++)
    {
        detectorList << master->getChild(j)->getTag() << " ";
        sameLevelDetectors.push_back(master->getChild(j));
    }
    
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining detector: " << (sameLevelDetectors.size() == 1 ? "" : "s") << ": " << detectorList.str() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    
    sendLog();
    
    int maxThreads = FileParser::getMaxThreads();
    boost::thread_group threads;

    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(refineDetectorWrapper, this, sameLevelDetectors, i, 1);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    reportProgress();
    manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
    
    boost::thread_group threads2;
    refinementEvent++;

    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(refineDetectorWrapper, this, sameLevelDetectors, i, 2);
        threads2.add_thread(thr);
    }
    
    threads2.join_all();

    reportProgress();
    manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);

    boost::thread_group threads3;
    refinementEvent++;
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(refineDetectorWrapper, this, sameLevelDetectors, i, 3);
        threads3.add_thread(thr);
    }
    
    threads3.join_all();

    
    reportProgress();
    manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset, int strategy)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < detectors.size(); i += maxThreads)
    {
        if (strategy == 1)
        {
            me->refineDetectorStrategyOne(detectors[i]);
        }
        if (strategy == 2)
        {
            me->refineDetectorStrategyTwo(detectors[i]);
        }
        if (strategy == 3)
        {
            me->refineDetectorStrategyThree(detectors[i]);
        }
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
    RefinementStrategyPtr refinementMap = makeRefiner(detector, type);
    
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointX,
                                Detector::setArrangedMidPointX, 2.0, 0.01, "midPointX");
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointY,
                                Detector::setArrangedMidPointY, 2.0, 0.01, "midPointY");
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                                Detector::setArrangedMidPointZ, 2.0, 0.01, "midPointZ");
    refinementMap->addParameter(&*detector, Detector::getAlpha, Detector::setAlpha, 0.01, 0.00001, "alpha");
    refinementMap->addParameter(&*detector, Detector::getBeta, Detector::setBeta, 0.01, 0.00001, "beta");
    refinementMap->addParameter(&*detector, Detector::getGamma, Detector::setGamma, 0.01, 0.00001, "gamma");
    
    refinementMap->refine();
}

void GeometryRefiner::refineMasterDetector()
{
    DetectorPtr detector = Detector::getMaster();
    refinementEvent++;

    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining detector: " << detector->getTag() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();
    
    refineMidPointXY(detector, GeometryScoreTypeMiller);
    reportProgress();
    refineMidPointZ(detector, GeometryScoreTypeIntrapanel);
    reportProgress();
    refineTiltXY(detector, GeometryScoreTypeMillerStdev);
    reportProgress();
    refineMidPointXY(detector, GeometryScoreTypeMiller);
    reportProgress();
    
    manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
}

void GeometryRefiner::refineDetectorStrategyOne(DetectorPtr detector)
{
    refineMidPointXY(detector, GeometryScoreTypeMiller);
    detector->millerScore(true, false);
}

void GeometryRefiner::refineDetectorStrategyTwo(DetectorPtr detector)
{
    refineTiltXY(detector, GeometryScoreTypeMillerStdev);
}

void GeometryRefiner::refineDetectorStrategyThree(DetectorPtr detector)
{
    refineMidPointXY(detector, GeometryScoreTypeMiller);
    //    refineVarious(detector, GeometryScoreTypeMiller);
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));
    manager->powderPattern("geom_refinement_event_0_start.csv", false);
    
}