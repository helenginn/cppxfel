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

std::string stringForScoreType(GeometryScoreType type)
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
        case GeometryScoreTypeIntrapanelParent:
            typeString = "intrapanel_parent";
            break;
        case GeometryScoreTypeInterMiller:
            typeString = "intra_miller";
            break;
        case GeometryScoreTypeIntraMiller:
            typeString = "inter_miller";
            break;
        default:
            break;
    }

    return typeString;
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
    
    std::string typeString = stringForScoreType(type);
    
    strategy->setJobName("Refining " + detector->getTag() + " (" + typeString + ")");

    strategy->setVerbose(false);
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    aManager->setActiveDetector(detector);
    aManager->setCycleNum(refinementEvent + 1);
    aManager->lockVectors();
    
    switch (type)
    {
        case GeometryScoreTypeInterMiller:
            strategy->setEvaluationFunction(Detector::millerScoreWrapper, &*detector);
            break;
        case GeometryScoreTypeIntraMiller:
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
    _changed = true;
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
    if (!_changed)
    {
        return;
    }
    
    _changed = false;
    
    std::string filename = "special_image_" + i_to_str(refinementEvent) + ".png";
    Detector::drawSpecialImage(filename);
    
    manager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    double intraScore = IndexManager::pseudoScore(&*manager);
    manager->setPseudoScoreType(PseudoScoreTypeAllInterPanel);
    double interScore = IndexManager::pseudoScore(&*manager);

    double intraIncrease = 100;
    double interIncrease = 100;
    double mScore = Detector::getMaster()->millerScore(false, false);
	Detector::getMaster()->millerScore(false, true);
    Detector::getMaster()->reportMillerScores(refinementEvent);

	intraIncrease = (intraScore / lastIntraScore - 1) * 100;
	interIncrease = (interScore / lastInterScore - 1) * 100;

    lastInterScore = interScore;
    lastIntraScore = intraScore;

    
    logged << "N: Progress score (event " << refinementEvent << ", intra-panel-dist): " << intraScore
    << " (" << (intraIncrease > 0 ? "+" : "") << intraIncrease << "% from last round) " << std::endl;
    logged << "N: Progress score (event " << refinementEvent << ", inter-panel-dist): " << interScore
    << " (" << (interIncrease > 0 ? "+" : "") << interIncrease << "% from last round)" << std::endl;

    if (mScore != 0)
    {
        logged << "N: Progress score (event " << refinementEvent << ", miller mean): " << mScore << std::endl;
    }
    sendLog();
    
    manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
    GeometryParser geomParser = GeometryParser("whatever", GeometryFormatCppxfel);
    geomParser.writeToFile("new_" + i_to_str(refinementEvent) + ".cppxfel_geom", refinementEvent);
    refinementEvent++;
    
    sendLog();
}

void GeometryRefiner::refineGeometry()
{
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
		intraPanelCycle();
		interPanelCycle();
	    refineBeamCentre();
    }

	logged << "**** Geometry refinement complete ****" << std::endl;
	logged << "**** Now setting METROLOGY_SEARCH_SIZE to 0 for convenience. ****" << std::endl;
	sendLog();

	FileParser::setKey("METROLOGY_SEARCH_SIZE", 0);
}


void GeometryRefiner::printHeader(std::vector<DetectorPtr> detectors)
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
}

void GeometryRefiner::intraPanelCycle()
{
	bool hasMillers = Detector::getMaster()->millerCount() > 0;
	std::vector<DetectorPtr> detectors;
	Detector::getMaster()->getAllSubDetectors(detectors, true);
	printHeader(detectors);

	GeometryScoreType type = GeometryScoreTypeIntrapanel;
	int cycles = 4;

	if (hasMillers)
	{
		type = GeometryScoreTypeIntraMiller;
		cycles = 4;
	}

	for (int i = 0; i < cycles; i++)
	{
		refineDetectorStrategyWrapper(this, detectors, type, 0);
	}
}

void GeometryRefiner::interPanelCycle()
{
	bool hasMillers = Detector::getMaster()->millerCount() > 0;

	GeometryScoreType type = GeometryScoreTypeInterpanel;

	if (hasMillers)
	{
		type = GeometryScoreTypeInterMiller;
		std::vector<DetectorPtr> detectors;
		Detector::getMaster()->getAllSubDetectors(detectors, true);
		printHeader(detectors);

		for (int i = 0; i < 4; i++)
		{
			refineDetectorStrategyWrapper(this, detectors, type, 0);
		}
	}
	else
	{
		int i = 64;
		while (i >= 0)
		{
			std::vector<DetectorPtr> detectors;
			detectors = Detector::getMaster()->getSubDetectorsOnLevel(i);
			i--;

			if (detectors.size() == 0)
			{
				continue;
			}

			refineDetectorStrategyWrapper(this, detectors, type, 0);
		}
	}

}

bool GeometryRefiner::geometryCycleForDetector(std::vector<DetectorPtr> detectors, bool interPanelOnly)
{
    bool hasRefineableDetectors = false;
    bool hasMillers = Detector::getMaster()->millerCount() > 0;
    GeometryScoreType type = GeometryScoreTypeIntrapanel;
    
    if (interPanelOnly)
    {
        type = GeometryScoreTypeInterpanel;
    }
    
    if (hasMillers && type == GeometryScoreTypeIntrapanel)
    {
        type = GeometryScoreTypeIntraMiller;
    }
    else if (hasMillers && type == GeometryScoreTypeInterpanel)
    {
        type = GeometryScoreTypeInterMiller;
    }

    refineUnitCell();

	for (int i = 0; i < 4; i++)
	{
		refineDetectorStrategyWrapper(this, detectors, type, 0);
		refineDetectorStrategyWrapper(this, detectors, type, 1);
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

    refineUnitCell();
    
    return hasRefineableDetectors;
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

void GeometryRefiner::refineDetectorStrategy(DetectorPtr detector, GeometryScoreType type, int strategyType)
{
    {
        RefinementStrategyPtr strategy = makeRefiner(detector, type);
        IndexManager::pseudoDistanceScore(&*detector->getIndexManager(), true, "before_");
    }
    
    std::string typeString = stringForScoreType(type);
    
    if (!detector->isRefinable(type))
    {
        return;
    }

	_changed = true;

	if (type == GeometryScoreTypeInterMiller)
	{
		if (detector->millerCount())
		{
			interPanelMillerSearch(detector, type);
		}
	}
	else if (type == GeometryScoreTypeIntraMiller)
	{
		intraPanelMillerSearch(detector, type);
	}
	else if (type == GeometryScoreTypeInterpanel)
	{
		interPanelGridSearch(detector);
		interPanelMillerSearch(detector, type);
	}
	else if (type == GeometryScoreTypeIntrapanel)
	{
		intraPanelMillerSearch(detector, type);
	}

	{
		RefinementStrategyPtr strategy = makeRefiner(detector, type);
        IndexManager::pseudoDistanceScore(&*detector->getIndexManager(), true, "after_");
    }
}

void GeometryRefiner::interPanelGridSearch(DetectorPtr detector)
{
	double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
	detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);

	RefinementGridSearchPtr strategy = makeGridRefiner(detector, GeometryScoreTypeInterpanel);
	detector->resetPoke();
	strategy->setJobName(detector->getTag() + "_powder");
	strategy->addParameter(&*detector->getChild(0), Detector::getInterNudgeX, Detector::setInterNudgeX, interNudge * 5, 0, "internudge_x");
	strategy->addParameter(&*detector->getChild(0), Detector::getInterNudgeY, Detector::setInterNudgeY, interNudge * 5, 0, "internudge_y");
	strategy->setGridLength(19);

	strategy->refine();

	detector->resetPoke();
}

void GeometryRefiner::interPanelNormalSearch(DetectorPtr detector)
{
 /*   double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);
    RefinementStrategyPtr strategy = makeRefiner(detector, type);
    detector->prepareInterNudges();
    strategy->addParameter(&*detector, Detector::getPokeX, Detector::setPokeX, interNudge / 10, 0.001, "poke_x");
    strategy->addParameter(&*detector, Detector::getPokeY, Detector::setPokeY, interNudge / 10, 0.001, "poke_y");
    strategy->addParameter(&*detector, Detector::getPokeZ, Detector::setPokeZ, interNudge / 10, 0.001, "poke_y");
    strategy->refine();
    
    detector->prepareInterNudges();
    
    
    return;*/
}

void GeometryRefiner::intraPanelMillerSearch(DetectorPtr detector, GeometryScoreType type)
{
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);
	nudgeStep /= 2;
	nudgeTiltX /= 2;
	nudgeTiltY /= 2;

    RefinementStrategyPtr strategy = makeRefiner(detector, type);
    detector->prepareInterNudges();
    strategy->setJobName(detector->getTag() + "_miller_stdev_z");
    strategy->addParameter(&*detector, Detector::getNudgeZ, Detector::setNudgeZ, nudgeStep, 0, "nudge_z");
    strategy->refine();
    strategy->clearParameters();
    detector->prepareInterNudges();
    strategy->addParameter(&*detector, Detector::getNudgeTiltX, Detector::setNudgeTiltX, interNudge, 0, "nudgetilt_x");
    strategy->addParameter(&*detector, Detector::getNudgeTiltY, Detector::setNudgeTiltY, interNudge, 0, "nudgetilt_y");
    strategy->setJobName(detector->getTag() + "_miller_stdev_tilt");
    strategy->refine();
    
    detector->prepareInterNudges();
}

void GeometryRefiner::interPanelMillerSearch(DetectorPtr detector, GeometryScoreType type)
{
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);
    
    RefinementStrategyPtr strategy = makeRefiner(detector, type);
    strategy->setJobName(detector->getTag() + "_miller");
    detector->resetPoke();

	if (type == GeometryScoreTypeInterMiller)
	{
		strategy->addParameter(&*detector, Detector::getInterNudgeX, Detector::setInterNudgeX, interNudge, 0, "internudge_x");
		strategy->addParameter(&*detector, Detector::getInterNudgeY, Detector::setInterNudgeY, interNudge, 0, "internudge_y");
		strategy->addParameter(&*detector, Detector::getInterNudgeZ, Detector::setInterNudgeZ, interNudge, 0, "internudge_z");
	}
	else
	{
		strategy->setJobName(detector->getTag() + "_powder");
		strategy->addParameter(&*detector->getChild(0), Detector::getInterNudgeX, Detector::setInterNudgeX, interNudge, 0, "internudge_x");
		strategy->addParameter(&*detector->getChild(0), Detector::getInterNudgeY, Detector::setInterNudgeY, interNudge, 0, "internudge_y");
		strategy->addParameter(&*detector->getChild(0), Detector::getInterNudgeZ, Detector::setInterNudgeZ, interNudge, 0, "internudge_z");
	}

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

	return;

    if (detector->millerCount() > 0)
    {
        return;
    }
    
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining beam centre: " << detector->getTag() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();
    
    RefinementGridSearchPtr gridSearch = makeGridRefiner(detector, GeometryScoreTypeBeamCentre);
    gridSearch->setJobName("grid_search_beam_centre_" + i_to_str(refinementEvent));
    gridSearch->setGridLength(99);
    
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);
    
    gridSearch->addParameter(&*detector, Detector::getInterNudgeX, Detector::setInterNudgeX, interNudge / 2, 0.001, "poke_x");
    gridSearch->addParameter(&*detector, Detector::getInterNudgeY, Detector::setInterNudgeY, interNudge / 2, 0.001, "poke_y");
    
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
    int num = 150;
    
    double step = FileParser::getKey("EXPECTED_GEOMETRY_MOVEMENT", 0.02);
    
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
    strategy->setVerbose(false);
    strategy->setJobName("Wide sweep detector " + detector->getTag());
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    aManager->setActiveDetector(detector);
    aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    strategy->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    
    strategy->refine();
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));    
}
