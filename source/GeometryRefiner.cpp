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
    aManager->setActiveDetector(detector);
    detector->setIndexManager(aManager);
    
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
        case GeometryScoreTypeAngleConsistency:
            aManager->setPseudoScoreType(PseudoScoreTypeAngleConsistency);
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
    lastInterAngleScore = 0;
    lastIntraAngleScore = 0;
}

void GeometryRefiner::refineGeometry()
{
    Detector::setNoisy(true);
    Detector::getMaster()->enableNudge();
    
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
        
        gridSearch(Detector::getMaster(), start, end);
        double newDistance = Detector::getArrangedMidPointZ(&*Detector::getMaster());
        
        logged << "**** Grid search done ****" << std::endl;
        logged << "**** New detector distance: " << newDistance * mmPerPixel << " mm. ****" << std::endl;
        sendLog();
        
        reportProgress();
    }
    
    for (int i = 0; i < maxCycles; i++)
    {
        cycleNum++;
        geometryCycleForDetector(detectors);
    }
}

void GeometryRefiner::reportProgress()
{
    std::string filename = "special_image_" + i_to_str(refinementEvent) + ".png";
    Detector::drawSpecialImage(filename);

    manager->setProportionDistance(1.0);
    manager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    double intraScore = -IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeAllInterPanel);
    double interScore = -IndexManager::pseudoScore(&*manager);
    
    manager->setProportionDistance(0.0);
    manager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    double intraAngle = -IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeAllInterPanel);
    double interAngle = -IndexManager::pseudoScore(&*manager);
    
    double intraIncrease = 100;
    double interIncrease = 100;
    double intraAngleIncrease = 100;
    double interAngleIncrease = 100;
    double mScore = Detector::getMaster()->millerScore(true, false);
    double sScore = Detector::getMaster()->millerScore(false, true);
    
    interIncrease = 100 * (interScore - lastInterScore) / interScore;
    intraIncrease = 100 * (intraScore - lastIntraScore) / intraScore;
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
    
    logged << "N: Progress score (event " << refinementEvent << ", intra-panel-angle): " << intraAngle
    << " (" << (intraAngleIncrease > 0 ? "+" : "") << intraAngleIncrease << "% from last round) " << std::endl;
    logged << "N: Progress score (event " << refinementEvent << ", inter-panel-angle): " << interAngle
    << " (" << (interAngleIncrease > 0 ? "+" : "") << interAngleIncrease << "% from last round)" << std::endl;
    
    if (mScore > 0 || sScore > 0)
    {
        logged << "N: Progress score (event " << refinementEvent << ", miller mean): " << mScore << std::endl;
        logged << "N: Progress score (event " << refinementEvent << ", miller stdev): " << sScore << std::endl;
    }
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
        refineBeamCentre();
    }
    else
    {
        refineDetectorStrategyWrapper(this, detectors, 0);
        refineDetectorStrategyWrapper(this, detectors, 1);
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
    
    refineDetectorStrategyWrapper(this, detectors, 2);
    
    if (detectors[0]->isLUCA())
    {
        refineMasterDetector();
        reportProgress();
        refineBeamCentre();
    }
    
    refineUnitCell();
}

void GeometryRefiner::refineDetectorStrategyWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int strategy)
{
    int maxThreads = FileParser::getMaxThreads();
    
    if (strategy == 1)
    {
        maxThreads = 1;
    }
    
    boost::thread_group threads;
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(refineDetectorWrapper, me, detectors, i, strategy);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    me->reportProgress();

}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset, int strategy)
{
    int maxThreads = FileParser::getMaxThreads();
    
    if (strategy == 1)
    {
        maxThreads = 1;
    }
    
    for (int i = offset; i < detectors.size(); i += maxThreads)
    {
        me->refineDetectorStrategy(detectors[i], strategy);
    }
}

void GeometryRefiner::refineDetector(DetectorPtr detector, GeometryScoreType type)
{
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
        default:
            break;
    }
    
    RefinementStrategyPtr strategy = makeRefiner(detector, type);
    strategy->setJobName("Refining " + detector->getTag() + " (" + typeString + ")");
    
    double transNudge = detector->isLUCA() ? 2.0 : 0.2;
    double tiltNudge = 0.001;
    double zTiltNudge = 0.0002;
    
    if (type == GeometryScoreTypeAngleConsistency)
    {
        detector->getIndexManager()->setProportionDistance(0.0);
        strategy->addParameter(&*detector, Detector::getGamma, Detector::setGamma, zTiltNudge, 0.000001, "tilt_z");
    }
    else if (type == GeometryScoreTypeIntrapanel)
    {
        strategy->addParameter(&*detector, Detector::getNudgeTiltX, Detector::setNudgeTiltX, tiltNudge, 0.00001, "nudge_tx");
        strategy->addParameter(&*detector, Detector::getNudgeTiltY, Detector::setNudgeTiltY, tiltNudge, 0.00001, "nudge_ty");
        strategy->addParameter(&*detector, Detector::getNudgeZ, Detector::setNudgeZ, transNudge, 0.001, "nudge_z");
        
    }
    else if (type == GeometryScoreTypeInterpanel)
    {
        detector->resetPoke();
        strategy->addParameter(&*detector, Detector::getPokeX, Detector::setPokeX, 0.1, 0.01, "poke_x");
        strategy->addParameter(&*detector, Detector::getPokeY, Detector::setPokeY, 0.1, 0.01, "poke_y");
        strategy->addParameter(&*detector, Detector::getPokeZ, Detector::setPokeZ, zTiltNudge, 0.000001, "poke_z");
    }
    else if (type == GeometryScoreTypeBeamCentre)
    {
        strategy->addParameter(&*detector, Detector::getNudgeX, Detector::setNudgeX, transNudge, 0.00001, "nudge_x");
        strategy->addParameter(&*detector, Detector::getNudgeY, Detector::setNudgeY, transNudge, 0.00001, "nudge_y");
    }
    
    strategy->refine();
    detector->resetPoke();
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
    reportProgress();
}

void GeometryRefiner::gridSearch(DetectorPtr detector, double start, double end)
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

void GeometryRefiner::refineDetectorStrategy(DetectorPtr detector, int strategy)
{
    if (strategy == 0)
    {
        if (!detector->hasChildren())
        {
            return;
        }
        
        refineDetector(detector, GeometryScoreTypeIntrapanel);
    }
    else if (strategy == 1)
    {
        if (!detector->hasChildren())
        {
            return;
        }
        
        refineDetector(detector, GeometryScoreTypeAngleConsistency);
    }
    else if (strategy == 2)
    {
        if (!detector->hasChildren() || !detector->getChild(0)->hasChildren())
        {
            return;
        }
        
        refineDetector(detector, GeometryScoreTypeInterpanel);
    }
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));    
}