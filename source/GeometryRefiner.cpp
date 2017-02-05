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
#include "UnitCellLattice.h"

RefinementStrategyPtr GeometryRefiner::makeRefiner(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr strategy = RefinementStrategy::userChosenStrategy();
    
    strategy->setVerbose(true);
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
        case GeometryScoreTypeBeamCentre:
            aManager->setPseudoScoreType(PseudoScoreTypeBeamCentre);
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
    Detector::getMaster()->enableNudge();

    reportProgress();
    std::vector<DetectorPtr> detectors;
    detectors.push_back(Detector::getMaster());
    
    for (int i = 0; i < 5; i++)
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
    
    if (intraScore > 0)
    {
        intraIncrease *= -1;
    }

    if (interScore > 0)
    {
        interIncrease *= -1;
    }
    
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
    
    refineUnitCell();

    if (detectors[0]->isLUCA())
    {
        refineMasterDetector();
        reportProgress();
    }
    else
    {
        int maxThreads = FileParser::getMaxThreads();
        boost::thread_group threads;
        
        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(refineDetectorWrapper, this, detectors, i, 0);
            threads.add_thread(thr);
        }
        /*
        for (int i = 0; i < detectors.size(); i++)
        {
            refineDetectorStrategy(detectors[i], 0);
            reportProgress();
        }*/
        
        threads.join_all();
        reportProgress();
        
        
        RefinementStrategyPtr strategy = makeRefiner(Detector::getMaster(), GeometryScoreTypeInterpanel);
        strategy->setJobName("nudge");
        
        for (int i = 0; i < detectors.size(); i++)
        {
            strategy->addParameter(&*detectors[i], Detector::getNudgeX, Detector::setNudgeX, 1.0, 0.01, "nudge_x");
            strategy->addParameter(&*detectors[i], Detector::getNudgeY, Detector::setNudgeY, 1.0, 0.01, "nudge_y");
            strategy->addParameter(&*detectors[i], Detector::getNudgeTiltZ, Detector::setNudgeTiltZ, 0.005, 0.000001, "nudge_tz");
        }
        
        if (detectors.size())
        {
            strategy->refine();
            Detector::getMaster()->lockNudges();
            
            for (int i = 0; i < detectors.size(); i++)
            {
                detectors[i]->lockNudges();
            }
            
            reportProgress();
        }
        
    }

    refineBeamCentre();
    reportProgress();
    
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

void GeometryRefiner::refineDetector(DetectorPtr detector, GeometryScoreType type)
{
    std::string typeString = (type == GeometryScoreTypeIntrapanel ? "intrapanel" : "beam_centre");

    RefinementStrategyPtr strategy = makeRefiner(detector, type);
    strategy->setJobName("Refining " + detector->getTag() + " (" + typeString + ")");

    if (type == GeometryScoreTypeIntrapanel)
    {
        strategy->addParameter(&*detector, Detector::getNudgeTiltX, Detector::setNudgeTiltX, 0.1, 0.00001, "nudge_tx");
        strategy->addParameter(&*detector, Detector::getNudgeTiltY, Detector::setNudgeTiltY, 0.1, 0.00001, "nudge_ty");
        strategy->addParameter(&*detector, Detector::getNudgeZ, Detector::setNudgeZ, 20.0, 0.001, "nudge_z");
    }
    else if (type == GeometryScoreTypeBeamCentre)
    {
        strategy->addParameter(&*detector, Detector::getNudgeX, Detector::setNudgeX, 5.0, 0.00001, "nudge_x");
        strategy->addParameter(&*detector, Detector::getNudgeY, Detector::setNudgeY, 5.0, 0.00001, "nudge_y");
    }
    
    strategy->refine();
    
    detector->lockNudges();
}

void GeometryRefiner::refineMasterDetector()
{
    DetectorPtr detector = Detector::getMaster();

    refineDetector(detector, GeometryScoreTypeIntrapanel);
}

void GeometryRefiner::refineBeamCentre()
{
    DetectorPtr detector = Detector::getMaster();
    
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining beam centre: " << detector->getTag() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();
    
    refineDetector(detector, GeometryScoreTypeBeamCentre);
    reportProgress();


}

void GeometryRefiner::refineUnitCell()
{
    bool optimisingUnitCellA = FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
    bool optimisingUnitCellB = FileParser::getKey("OPTIMISING_UNIT_CELL_B", false);
    bool optimisingUnitCellC = FileParser::getKey("OPTIMISING_UNIT_CELL_C", false);
    bool optimisingUnitCellAlpha = FileParser::getKey("OPTIMISING_UNIT_CELL_ALPHA", false);
    bool optimisingUnitCellBeta = FileParser::getKey("OPTIMISING_UNIT_CELL_BETA", false);
    bool optimisingUnitCellGamma = FileParser::getKey("OPTIMISING_UNIT_CELL_GAMMA", false);
    
    if (!optimisingUnitCellA && !optimisingUnitCellAlpha &&
        !optimisingUnitCellB && !optimisingUnitCellBeta &&
        !optimisingUnitCellC && !optimisingUnitCellGamma)
    {
        return;
    }
    
    DetectorPtr detector = Detector::getMaster();
    
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining unit cell: " << detector->getTag() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();
    
    RefinementStrategyPtr strategy = makeRefiner(detector, GeometryScoreTypeIntrapanel);
    strategy->setJobName("Refining unit cell");
    IndexManager *manager = static_cast<IndexManager *>(strategy->getEvaluationObject());

    if (optimisingUnitCellA)
    {
        strategy->addParameter(&*(manager->getLattice()), UnitCellLattice::getUnitCellA, UnitCellLattice::setUnitCellA, 0.1, 0.00001, "unit_cell_a");
    }
    if (optimisingUnitCellB)
    {
        strategy->addParameter(&*(manager->getLattice()), UnitCellLattice::getUnitCellB, UnitCellLattice::setUnitCellB, 0.1, 0.00001, "unit_cell_b");
    }
    if (optimisingUnitCellC)
    {
        strategy->addParameter(&*(manager->getLattice()), UnitCellLattice::getUnitCellC, UnitCellLattice::setUnitCellC, 0.1, 0.00001, "unit_cell_c");
    }

    if (optimisingUnitCellAlpha)
    {
        strategy->addParameter(&*(manager->getLattice()), UnitCellLattice::getUnitCellAlpha, UnitCellLattice::setUnitCellAlpha, 0.1, 0.00001, "unit_cell_alpha");
    }
    if (optimisingUnitCellBeta)
    {
        strategy->addParameter(&*(manager->getLattice()), UnitCellLattice::getUnitCellBeta, UnitCellLattice::setUnitCellBeta, 0.1, 0.00001, "unit_cell_beta");
    }
    if (optimisingUnitCellGamma)
    {
        strategy->addParameter(&*(manager->getLattice()), UnitCellLattice::getUnitCellGamma, UnitCellLattice::setUnitCellGamma, 0.1, 0.00001, "unit_cell_gamma");
    }
    
    strategy->setCycles(15);
    
    strategy->refine();
    detector->lockNudges();
    reportProgress();
}

void GeometryRefiner::gridSearch(DetectorPtr detector)
{
    RefinementGridSearchPtr strategy = RefinementGridSearchPtr(new RefinementGridSearch());
    Detector::setAlpha(&*Detector::getMaster(), 0.001);
    Detector::setBeta(&*Detector::getMaster(), 0.001);
    
    strategy->addParameter(&*detector, Detector::getGamma, Detector::setGamma, 0.05, 0.00001, "nudge_tx");
//    strategy->addParameter(&*detector, Detector::getNudgeY, Detector::setNudgeY, 0.001, 0.00001, "nudge_ty");
//    strategy->addParameter(&*detector, Detector::getNudgeTiltZ, Detector::setNudgeTiltZ, 0.05, 0.001, "nudge_tz");
    
    strategy->setVerbose(true);
    strategy->setJobName("Detector " + detector->getTag());
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    indexManagers.push_back(aManager);
    aManager->setActiveDetector(detector);
    aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    strategy->setEvaluationFunction(IndexManager::debugPseudoScore, &*aManager);
    
    strategy->refine();
    
    exit(0);
}

void GeometryRefiner::refineDetectorStrategy(DetectorPtr detector, int strategy)
{
    if (strategy == 0)
    {
        refineDetector(detector, GeometryScoreTypeIntrapanel);
        detector->millerScore(true, false);
    }
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));
    manager->powderPattern("geom_refinement_event_0_start.csv", false);
    
}