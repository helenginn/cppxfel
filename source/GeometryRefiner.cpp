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
    RefinementStrategyPtr strategy;
    
    int methodInt = FileParser::getKey("MINIMIZATION_METHOD", 0);
    
    if (methodInt == 0)
    {
        RefinementStepSearchPtr stepSearch = RefinementStepSearchPtr(new RefinementStepSearch());
        strategy = boost::static_pointer_cast<RefinementStrategy>(stepSearch);
    }
    else
    {
        NelderMeadPtr nelderMead = NelderMeadPtr(new NelderMead());
        strategy = boost::static_pointer_cast<RefinementStrategy>(nelderMead);
    }
    
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
    for (int i = 0; i < 2; i++)
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

    while (nextLevelDetectors.size())
    {
        sameLevelDetectors = nextLevelDetectors;
        nextLevelDetectors.clear();
        
        std::ostringstream detectorList;
        
        for (int i = 0; i < sameLevelDetectors.size(); i++)
        {
            detectorList << sameLevelDetectors[i]->getTag() << " ";
            
            for (int j = 0; j < sameLevelDetectors[i]->childrenCount(); j++)
            {
                nextLevelDetectors.push_back(sameLevelDetectors[i]->getChild(j));
            }
        }
        
        refinementEvent++;

        logged << "***************************************************" << std::endl;
        logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
        logged << "  Refining detector" << (sameLevelDetectors.size() == 1 ? "" : "s") << ": " << detectorList.str() << std::endl;
        logged << "***************************************************" << std::endl << std::endl;
        
        sendLog();
        
        double maxThreads = FileParser::getMaxThreads();
        boost::thread_group threads;
        
        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(refineDetectorWrapper, this, sameLevelDetectors, i);
            threads.add_thread(thr);
        }
        
        threads.join_all();
        reportProgress();
        manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
    }
    
}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < detectors.size(); i += maxThreads)
    {
        me->refineDetector(detectors[i]);
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

    reportProgress();
}

void GeometryRefiner::refineMidPointZ(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr refinementMap = makeRefiner(detector, type);
    
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                                Detector::setArrangedMidPointZ, 1.0, 0.01, "midPointZ");
    
    refinementMap->refine();
    
    reportProgress();
}

void GeometryRefiner::refineTiltXY(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr refinementMap = makeRefiner(detector, type);
    
    refinementMap->addParameter(&*detector, Detector::getAlpha, Detector::setAlpha, 0.001, 0.00001, "alpha");
    refinementMap->addParameter(&*detector, Detector::getBeta, Detector::setBeta, 0.001, 0.00001, "beta");
    
    refinementMap->refine();
    
    reportProgress();
}

void GeometryRefiner::refineMasterDetector()
{
    DetectorPtr detector = Detector::getMaster();
    
    refineMidPointXY(detector, GeometryScoreTypeMiller);
    refineMidPointZ(detector, GeometryScoreTypeIntrapanel);
 //   refineTiltXY(detector, GeometryScoreTypeIntrapanel);
 //   refineMidPointXY(detector, GeometryScoreTypeMiller);
    
    reportProgress();
}

void GeometryRefiner::refineDetector(DetectorPtr detector)
{
    if (detector->isLUCA())
    {
        refineMasterDetector();
        return;
    }
    
    RefinementStrategyPtr refinementMap = makeRefiner(detector, GeometryScoreTypeMiller);
    
/*
    refinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                               Detector::setArrangedMidPointZ, 0.05, 0.001, "midPointZ");
    refinementMap->addParameter(&*detector, Detector::getAlpha, Detector::setAlpha, 0.00005, 0.000001, "alpha");
    refinementMap->addParameter(&*detector, Detector::getBeta, Detector::setBeta, 0.00005, 0.000001, "beta");
*/
    
    refinementMap->refine();
    
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));
    manager->powderPattern("geom_refinement_event_0_start.csv", false);
    
}