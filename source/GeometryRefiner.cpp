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

PseudoScoreType pseudoScoreTypeForGeometryType(GeometryScoreType type)
{
    switch (type)
    {
        case GeometryScoreTypeInterpanel:
            return PseudoScoreTypeInterPanel;
        case GeometryScoreTypeBeamCentre:
            return PseudoScoreTypeBeamCentre;
        case GeometryScoreTypeIntrapanel:
            return PseudoScoreTypeIntraPanel;
        default:
            return PseudoScoreTypeInvalid;
    }
}

RefinementGridSearchPtr GeometryRefiner::makeGridRefiner(DetectorPtr detector, GeometryScoreType type)
{
    RefinementGridSearchPtr strategy = RefinementGridSearchPtr(new RefinementGridSearch());
    strategy->setGridLength(31);
    strategy->setVerbose(false);
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    aManager->setActiveDetector(detector);
    aManager->lockVectors();
    aManager->setPseudoScoreType(pseudoScoreTypeForGeometryType(type));
    strategy->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    
    return strategy;
}

RefinementStrategyPtr GeometryRefiner::makeRefiner(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr strategy = RefinementStrategy::userChosenStrategy();
    
    std::string typeString = "";
    
    switch (type) {
        case GeometryScoreTypeBeamCentre:
            typeString = "beam_centre";
            break;
        case GeometryScoreTypeInterpanel:
            typeString = "interpanel";
            break;
        case GeometryScoreTypeIntrapanel:
            typeString = "intrapanel";
            break;
        case GeometryScoreTypeIntrapanelParent:
            typeString = "intrapanel_parent";
            break;
        default:
            break;
    }
    
    strategy->setJobName("Refining " + detector->getTag() + " (" + typeString + ")");

    strategy->setVerbose(true);
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    aManager->setActiveDetector(detector);
    aManager->setCycleNum(refinementEvent + 1);
    aManager->lockVectors();
    
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
        default:
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
    lastInterAngleScore = 0;
    lastIntraAngleScore = 0;

}

void GeometryRefiner::refineGeometry()
{
    Detector::getMaster()->enableNudge();
    
  //  startingGraphs();
    reportProgress();
    
    std::vector<DetectorPtr> detectors;
    detectors.push_back(Detector::getMaster());
    
    int maxCycles = FileParser::getKey("MAXIMUM_CYCLES", 6);
    
    std::vector<double> sweepDetectorDistance = FileParser::getKey("SWEEP_DETECTOR_DISTANCE", std::vector<double>());
    
    if (sweepDetectorDistance.size() >= 2)
    {
        
        double start = sweepDetectorDistance[0];
        double end = sweepDetectorDistance[1];
        
        if (end < start)
        {
            end = sweepDetectorDistance[0];
            start = sweepDetectorDistance[1];
        }
        
        logged << "***********************************************" << std::endl;
        logged << "Doing a detector distance sweep from " << start << " to " << end << " mm." << std::endl;
        logged << "***********************************************" << std::endl;
        sendLog();
        
        double mmPerPixel = FileParser::getKey("MM_PER_PIXEL", 0.11);

        end /= mmPerPixel;
        start /= mmPerPixel;
        
        gridSearchDetectorDistance(Detector::getMaster(), start, end);
        double newDistance = Detector::getArrangedMidPointZ(&*Detector::getMaster());
        
        logged << "**** Grid search done ****" << std::endl;
        logged << "**** New detector distance: " << newDistance * mmPerPixel << " mm. ****" << std::endl;
        sendLog();
        
        reportProgress();
    }
    
    for (int i = 0; i < maxCycles; i++)
    {
        cycleNum++;
        geometryCycleForDetector(detectors, false);
        refineBeamCentre();
    }
}

void GeometryRefiner::startingGraphs()
{
    logged << "Generating starting graphs..." << std::endl;
    sendLog();
    
    std::vector<DetectorPtr> allDetectors;
    Detector::getMaster()->getAllSubDetectors(allDetectors, false);
    
    for (int i = 0; i < allDetectors.size(); i++)
    {
        DetectorPtr detector = allDetectors[i];
        IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
        aManager->setActiveDetector(detector);
        aManager->setCycleNum(0);
        aManager->lockVectors();
        
        if (detector->isRefinable(GeometryScoreTypeIntrapanel))
        {
            aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
            IndexManager::pseudoDistanceScore(&*aManager, true, "start_");
            aManager->plotGoodVectors();
        }

        if (detector->isRefinable(GeometryScoreTypeInterpanel))
        {
            aManager->clearGoodVectors();
            aManager->setPseudoScoreType(PseudoScoreTypeInterPanel);
            aManager->resetMaxFrequency();
            IndexManager::pseudoDistanceScore(&*aManager, true, "start_");
            aManager->plotGoodVectors();
        }
    }
}

void GeometryRefiner::reportProgress()
{
    double maxAngleDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.);
    std::string filename = "special_image_" + i_to_str(refinementEvent) + ".png";
    Detector::drawSpecialImage(filename);

    manager->setProportionDistance(1.0);
    manager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    double intraScore = IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeAllInterPanel);
    double interScore = IndexManager::pseudoScore(&*manager);
    
    manager->setProportionDistance(0.0);
    manager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    double intraAngle = -IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeAllInterPanel);
    double interAngle = -IndexManager::pseudoScore(&*manager);
    
    manager->pseudoAngleCSV();

    double intraIncrease = 100;
    double interIncrease = 100;
    double intraAngleIncrease = 100;
    double interAngleIncrease = 100;
    double mScore = Detector::getMaster()->millerScore(true, false);
    double sScore = Detector::getMaster()->millerScore(false, true);
    
    interIncrease = -100 * (interScore - lastInterScore) / interScore;
    intraIncrease = -100 * (intraScore - lastIntraScore) / intraScore;
    interAngleIncrease = 100 * (interAngle - lastInterAngleScore) / interAngle;
    intraAngleIncrease = 100 * (intraAngle - lastIntraAngleScore) / intraAngle;
    
    lastInterScore = interScore;
    lastIntraScore = intraScore;
    lastInterAngleScore = interAngle;
    lastIntraAngleScore = intraAngle;
    

    logged << "N: Progress score (event " << refinementEvent << ", intra-panel-dist): " << intraScore
    << " (" << (intraIncrease > 0 ? "+" : "") << intraIncrease << "% from last round) " << std::endl;
    logged << "N: Progress score (event " << refinementEvent << ", inter-panel-dist): " << interScore
    << " (" << (interIncrease > 0 ? "+" : "") << interIncrease << "% from last round)" << std::endl;
    
    if (maxAngleDistance > 0)
    {
        logged << "N: Progress score (event " << refinementEvent << ", intra-panel-angle): " << intraAngle
        << " (" << (intraAngleIncrease > 0 ? "+" : "") << intraAngleIncrease << "% from last round) " << std::endl;
        logged << "N: Progress score (event " << refinementEvent << ", inter-panel-angle): " << interAngle
        << " (" << (interAngleIncrease > 0 ? "+" : "") << interAngleIncrease << "% from last round)" << std::endl;
    }
    
    if (mScore > 0 || sScore > 0)
    {
        logged << "N: Progress score (event " << refinementEvent << ", miller mean): " << mScore << std::endl;
        logged << "N: Progress score (event " << refinementEvent << ", miller stdev): " << sScore << std::endl;
    }
    sendLog();

    manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
    GeometryParser geomParser = GeometryParser("whatever", GeometryFormatCppxfel);
    geomParser.writeToFile("new_" + i_to_str(refinementEvent) + ".cppxfel_geom", refinementEvent);
    refinementEvent++;
    
    sendLog();
}

void GeometryRefiner::printHeader(std::vector<DetectorPtr> detectors, std::string detectorList)
{
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining detector" << (detectors.size() == 1 ? "" : "s") << ": " << detectorList << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();

}

bool GeometryRefiner::geometryCycleForDetector(std::vector<DetectorPtr> detectors, bool interPanelOnly)
{
    std::ostringstream detectorList;
    
    bool hasRefineableDetectors = false;
    GeometryScoreType type = GeometryScoreTypeIntrapanel;

    for (int j = 0; j < detectors.size(); j++)
    {
        detectorList << detectors[j]->getTag() << " ";
        
        if (detectors[j]->isRefinable(GeometryScoreTypeIntrapanel) && !hasRefineableDetectors)
        {
            hasRefineableDetectors = true;
            type = GeometryScoreTypeIntrapanel;
        }
        
        if (interPanelOnly)
        {
            hasRefineableDetectors = false;
        }
    }
    
    if (hasRefineableDetectors && type == GeometryScoreTypeIntrapanel)
    {
        printHeader(detectors, detectorList.str());
    }
    
    refineUnitCell();

    if (hasRefineableDetectors && type == GeometryScoreTypeIntrapanel)
    {
        for (int i = 0; i < 3; i++)
        {
            refineDetectorStrategyWrapper(this, detectors, type, 0);
            refineDetectorStrategyWrapper(this, detectors, type, 1);
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
        bool somethingHappened = geometryCycleForDetector(nextDetectors, interPanelOnly);
        
        if (!somethingHappened)
        {
            return true;
        }
    }
    
    for (int i = 0; i < 1; i++)
    {
        refineDetectorStrategyWrapper(this, detectors, GeometryScoreTypeInterpanel, 0);
        refineDetectorStrategyWrapper(this, detectors, GeometryScoreTypeInterpanel, 1);
    }
    
    refineUnitCell();
    
    return true;
}

void GeometryRefiner::refineDetectorStrategyWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, GeometryScoreType type, int strategyType)
{
    int maxThreads = FileParser::getMaxThreads();
   
    {
        boost::thread_group threads;
        
        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(refineDetectorWrapper, me, detectors, i, type, strategyType);
            threads.add_thread(thr);
        }
        
        threads.join_all();
    }
    
    me->reportProgress();
    
    me->firstCycle = false;
}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset, GeometryScoreType type, int strategyType)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < detectors.size(); i += maxThreads)
    {
        me->refineDetectorStrategy(detectors[i], type, strategyType);
    }
}

void GeometryRefiner::interPanelGridSearch(DetectorPtr detector)
{
    RefinementGridSearchPtr strategy = makeGridRefiner(detector, GeometryScoreTypeInterpanel);
    
    detector->resetPoke();
    double distFromSample = detector->distanceFromSample();
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);

    strategy->addParameter(&*detector, Detector::getPokeX, Detector::setPokeX, interNudge / 2, 0.001, "poke_x");
    strategy->addParameter(&*detector, Detector::getPokeY, Detector::setPokeY, interNudge / 2, 0.001, "poke_y");
    strategy->setGridLength(99);
    
    double pixelPerTiltNudge = sin(nudgeTiltY) * distFromSample;
    double pixelsToCheck = 1.9;
    int gridPointsToCheck = pixelsToCheck / pixelPerTiltNudge + 0.5;
    gridPointsToCheck *= 2;
    
    strategy->setCheckGridNum(gridPointsToCheck);
    
    strategy->setJobName("interpanel_grid_" + detector->getTag() + "_" + i_to_str(refinementEvent));
    
    strategy->refine();
    strategy->assignInterpanelMinimum();
}

void GeometryRefiner::interPanelNormalSearch(DetectorPtr detector)
{
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);
    
    RefinementStrategyPtr strategy = makeRefiner(detector, GeometryScoreTypeInterpanel);
    detector->resetPoke();
  
    strategy->addParameter(&*detector, Detector::getPokeX, Detector::setPokeX, interNudge / 20, 0.000001, "poke_x");
    strategy->addParameter(&*detector, Detector::getPokeY, Detector::setPokeY, interNudge / 20, 0.000001, "poke_y");
    strategy->addParameter(&*detector, Detector::getPokeZ, Detector::setPokeZ, interNudge / 20, 0.000001, "poke_z");
    strategy->refine();
    
    detector->resetPoke();
}

void GeometryRefiner::intraPanelNormalSearch(DetectorPtr detector)
{
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);

    RefinementStrategyPtr strategy = makeRefiner(detector, GeometryScoreTypeIntrapanel);
    detector->prepareInterNudges();
    strategy->addParameter(&*detector, Detector::getNudgeTiltX, Detector::setSmartTiltX, nudgeTiltX, 0.1, "nudgetilt_x");
    strategy->addParameter(&*detector, Detector::getNudgeTiltY, Detector::setSmartTiltY, nudgeTiltY, 0.1, "nudgetilt_y");
    strategy->refine();
    
    detector->prepareInterNudges();

}

void GeometryRefiner::refineBeamCentre()
{
    DetectorPtr detector = Detector::getMaster();
    
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining beam centre: " << detector->getTag() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();
    
    RefinementGridSearchPtr gridSearch = makeGridRefiner(detector, GeometryScoreTypeBeamCentre);
    gridSearch->setJobName("grid_search_beam_centre_" + i_to_str(refinementEvent));
    gridSearch->setGridLength(99);
    
    double nudgeStep, nudgeTiltX, nudgeTiltY;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep);
    
    gridSearch->addParameter(&*detector, Detector::getInterNudgeX, Detector::setInterNudgeX, nudgeTiltX, 0.001, "poke_x");
    gridSearch->addParameter(&*detector, Detector::getInterNudgeY, Detector::setInterNudgeY, nudgeTiltY, 0.001, "poke_y");
    
    detector->getIndexManager()->setPseudoScoreType(PseudoScoreTypeBeamCentre);
    gridSearch->refine();
    
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
    reportProgress();
}

void GeometryRefiner::nudgeZGridSearch(DetectorPtr detector)
{
    int num = 50;
    
    double step = FileParser::getKey("EXPECTED_GEOMETRY_MOVEMENT", 0.5);
    
    double confidence = 2 * num + 1;
    
    {
        RefinementGridSearchPtr strategy = makeGridRefiner(detector, GeometryScoreTypeIntrapanel);
        strategy->addParameter(&*detector, Detector::getNudgeZ, Detector::setNudgeZ, step, 0.1, "nudge_z");
        strategy->setGridLength(confidence);
        strategy->setJobName("nudge_Z_wide_sweep_" + detector->getTag());
        strategy->refine();
    }
}

void GeometryRefiner::gridSearchDetectorDistance(DetectorPtr detector, double start, double end)
{
    double step = 0.2;
    double middle = (end + start) / 2;
    Detector::setArrangedMidPointZ(&*detector, middle);
    double confidence = (end - start) / step / 2;
    
    RefinementGridSearchPtr strategy = RefinementGridSearchPtr(new RefinementGridSearch());
    strategy->addParameter(&*detector, Detector::getArrangedMidPointZ, Detector::setArrangedMidPointZ, step, 0.1, "nudge_z");
    
    strategy->setGridLength(confidence * 2 + 1);
    strategy->setVerbose(true);
    strategy->setJobName("Wide sweep detector " + detector->getTag());
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    aManager->setActiveDetector(detector);
    aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    strategy->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    
    strategy->refine();
}

void GeometryRefiner::refineDetectorStrategy(DetectorPtr detector, GeometryScoreType type, int strategyType)
{
    {
        RefinementStrategyPtr strategy = makeRefiner(detector, type);
        IndexManager::pseudoDistanceScore(&*detector->getIndexManager(), true, "before_");
    }

    if (type == GeometryScoreTypeInterpanel)
    {
        if (!detector->isRefinable(type))
            return;
       
        if (strategyType == 0)
        {
            interPanelGridSearch(detector);
        }
        else if (strategyType == 1)
        {
            interPanelNormalSearch(detector);
        }
    }
    else if (type == GeometryScoreTypeIntrapanel)
    {
        if (!detector->isRefinable(type))
            return;
        
        if (strategyType == 0)
        {
            nudgeZGridSearch(detector);
        }
        else if (strategyType == 1)
        {
            intraPanelNormalSearch(detector);
        }
    }

    {
        RefinementStrategyPtr strategy = makeRefiner(detector, type);
        IndexManager::pseudoDistanceScore(&*detector->getIndexManager(), true, "after_");
    }
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));    
}
