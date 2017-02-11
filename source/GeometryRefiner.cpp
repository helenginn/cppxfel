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
    
    
     logged << "********* TEST **********" << std::endl;
    sendLog();
    
//    Detector::setAlpha(&*(Detector::getMaster()), 0.1);
//    Detector::setArrangedMidPointX(&*(Detector::getMaster()), 5.0);
/*    Detector::getMaster()->description();
    Detector::getMaster()->getChild(0)->getChild(0)->getChild(1)->getChild(0)->getChild(0)->getChild(0)->description();
    Detector::setNudgeTiltX(&*(Detector::getMaster()), 0.005);
    Detector::setNudgeTiltY(&*(Detector::getMaster()), 0.5);
    Detector::setNudgeTiltZ(&*(Detector::getMaster()), 0.005);
    logged << "Setting nudge stuff" << std::endl;
    sendLog();
    Detector::getMaster()->description();
    //Detector::getMaster()->getChild(0)->getChild(0)->getChild(1)->getChild(0)->getChild(0)->getChild(0)->description();
    logged << "Locking nudges now" << std::endl;
    sendLog();
    Detector::getMaster()->lockNudges();
    Detector::getMaster()->description();
    //Detector::getMaster()->getChild(0)->getChild(0)->getChild(1)->getChild(0)->getChild(0)->getChild(0)->description();
    */
    
  //  gridSearch(Detector::getMaster());
    
    
 //   exit(0);
     
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
    manager->setPseudoScoreType(PseudoScoreTypeAllInterPanel);
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
    int maxThreads = FileParser::getMaxThreads();

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
        Detector::getMaster()->lockNudges();
        reportProgress();
    }
    else
    {
        boost::thread_group threads;
        
        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(refineDetectorWrapper, this, detectors, i, 0);
            threads.add_thread(thr);
        }
        
        threads.join_all();
        reportProgress();
        
        Detector::getMaster()->lockNudges();
        
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
    
    boost::thread_group threads;
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(refineDetectorWrapper, this, detectors, i, 1);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    reportProgress();
    
    if (detectors[0]->isLUCA())
    {
        refineMasterDetector();
        Detector::getMaster()->lockNudges();
        reportProgress();
    }
    
    refineUnitCell();
    
    reportProgress();
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
    
    double zNudge = (detector->isLUCA() ? 2.0 : 0.2);

    if (type == GeometryScoreTypeIntrapanel)
    {
        strategy->addParameter(&*detector, Detector::getNudgeTiltX, Detector::setNudgeTiltX, 0.001, 0.00001, "nudge_tx");
        strategy->addParameter(&*detector, Detector::getNudgeTiltY, Detector::setNudgeTiltY, 0.001, 0.00001, "nudge_ty");
        strategy->addParameter(&*detector, Detector::getNudgeZ, Detector::setNudgeZ, zNudge, 0.001, "nudge_z");
        strategy->addParameter(&*detector, Detector::getNudgeTiltZ, Detector::setNudgeTiltZ, 0.001, 0.000001, "nudge_tz");
    }
    else if (type == GeometryScoreTypeInterpanel)
    {
        for (int i = 0; i < detector->childrenCount(); i++)
        {
            DetectorPtr child = detector->getChild(i);
            strategy->addParameter(&*child, Detector::getNudgeX, Detector::setNudgeX, 0.2, 0.01, "nudge_x");
            strategy->addParameter(&*child, Detector::getNudgeY, Detector::setNudgeY, 0.2, 0.01, "nudge_y");
     //       strategy->addParameter(&*child, Detector::getNudgeTiltZ, Detector::setNudgeTiltZ, 0.001, 0.000001, "nudge_tz");
        }
    }
    else if (type == GeometryScoreTypeBeamCentre)
    {
        strategy->addParameter(&*detector, Detector::getNudgeX, Detector::setNudgeX, 2.0, 0.00001, "nudge_x");
        strategy->addParameter(&*detector, Detector::getNudgeY, Detector::setNudgeY, 2.0, 0.00001, "nudge_y");
    }
    
    strategy->refine();
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
 //   Detector::setAlpha(&*Detector::getMaster(), 0.001);
 //   Detector::setBeta(&*Detector::getMaster(), 0.001);
    
//    strategy->addParameter(&*detector, Detector::getGamma, Detector::setGamma, 0.05, 0.00001, "nudge_tx");
    strategy->addParameter(&*detector, Detector::getNudgeY, Detector::setNudgeY, 2.0, 0.00001, "nudge_x");
    strategy->addParameter(&*detector, Detector::getNudgeX, Detector::setNudgeX, 2.0, 0.00001, "nudge_y");
//    strategy->addParameter(&*detector, Detector::getNudgeTiltZ, Detector::setNudgeTiltZ, 0.05, 0.001, "nudge_tz");
    
    strategy->setGridLength(5);
    strategy->setVerbose(true);
    strategy->setJobName("Detector " + detector->getTag());
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    indexManagers.push_back(aManager);
    aManager->setActiveDetector(detector);
    aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    strategy->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    
    strategy->refine();
    
    logged << "****** PRELOCK ******" << std::endl;
    sendLog();
    Detector::getMaster()->fullDescription();
    logged << "Intrapanel: " << IndexManager::pseudoScore(&*aManager);
    sendLog();
    
    Detector::getMaster()->lockNudges();

    logged << "****** POSTLOCK ******" << std::endl;
    sendLog();
    Detector::setNudgeTiltX(&*(Detector::getMaster()), 0);
    Detector::getMaster()->fullDescription();
    logged << "Intrapanel: " << IndexManager::pseudoScore(&*aManager);
    sendLog();
    
    strategy->setJobName("job2");
    strategy->refine();
    
    exit(0);
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
        detector->millerScore(true, false);
    }
    else if (strategy == 1)
    {
        if (!detector->hasChildren() || !detector->getChild(0)->hasChildren())
        {
            return;
        }
        
        refineDetector(detector, GeometryScoreTypeInterpanel);
        detector->millerScore(true, false);
    }
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));
    manager->powderPattern("geom_refinement_event_0_start.csv", false);
    
}