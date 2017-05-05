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

	if (type == GeometryScoreTypePeakSearch)
	{
		strategy->setEvaluationFunction(Detector::peakScoreWrapper, &*detector);
		return strategy;
	}

	if (type == GeometryScoreTypeInterpanel || type == GeometryScoreTypeIntrapanel
		|| type == GeometryScoreTypeBeamCentre)
	{
		IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));

		if (detector->getIndexManager(type))
		{
			aManager = detector->getIndexManager(type);
		}
		else
		{
			aManager->setActiveDetector(detector, type);
			aManager->setCycleNum(refinementEvent);
			aManager->lockVectors();
			aManager->setPseudoScoreType(pseudoScoreTypeForGeometryType(type));
		}

		strategy->setEvaluationFunction(IndexManager::pseudoAngleScore, &*aManager);

		return strategy;
	}
	else
	{
		switch (type)
		{
			case GeometryScoreTypeInterMiller:
				strategy->setEvaluationFunction(Detector::millerScoreWrapper, &*detector);
				break;
			case GeometryScoreTypeIntraMiller:
				strategy->setEvaluationFunction(Detector::millerStdevScoreWrapper, &*detector);
				break;
			default:
				break;
		}

		return strategy;
	}
}

RefinementStrategyPtr GeometryRefiner::makeRefiner(DetectorPtr detector, GeometryScoreType type)
{
    RefinementStrategyPtr strategy = RefinementStrategy::userChosenStrategy();
    
    std::string typeString = stringForScoreType(type);
    
    strategy->setJobName("Refining " + detector->getTag() + " (" + typeString + ")");

    strategy->setVerbose(false);
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));

	// Reuse an old manager if the score type is correct
	if (detector->getIndexManager(type))
	{
		aManager = detector->getIndexManager(type);
	}
	else
	{
		aManager->setActiveDetector(detector, type);
		aManager->setCycleNum(refinementEvent);
		aManager->lockVectors();
	}

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
        aManager->setActiveDetector(detector, GeometryScoreTypeIntrapanel);
        aManager->setCycleNum(0);
        aManager->lockVectors();
        
        if (detector->isRefinable(GeometryScoreTypeIntrapanel))
        {
            aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
            IndexManager::pseudoAngleScore(&*aManager);
            aManager->plotGoodVectors();
        }
        
        if (detector->isRefinable(GeometryScoreTypeInterpanel))
        {
            aManager->clearGoodVectors();
            aManager->setPseudoScoreType(PseudoScoreTypeInterPanel);
            IndexManager::pseudoAngleScore(&*aManager);
            aManager->plotGoodVectors();
        }

		detector->setIndexManager(IndexManagerPtr(), GeometryScoreTypeInterpanel);
		detector->setIndexManager(IndexManagerPtr(), GeometryScoreTypeIntrapanel);
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

	IndexManagerPtr intraManager = Detector::getMaster()->getIndexManager(GeometryScoreTypeIntrapanel);
	if (!intraManager)
	{
		intraManager = IndexManagerPtr(new IndexManager(images));
		intraManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
		Detector::getMaster()->setIndexManager(intraManager, GeometryScoreTypeIntrapanel);
	}

	IndexManagerPtr interManager = Detector::getMaster()->getIndexManager(GeometryScoreTypeInterpanel);
	if (!interManager)
	{
		interManager = IndexManagerPtr(new IndexManager(images));
		interManager->setPseudoScoreType(PseudoScoreTypeAllInterPanel);
		Detector::getMaster()->setIndexManager(interManager, GeometryScoreTypeInterpanel);
	}

    double intraScore = IndexManager::pseudoScore(&*intraManager);
    double interScore = IndexManager::pseudoScore(&*interManager);

    double intraIncrease = 100;
    double interIncrease = 100;
    double mScore = Detector::getMaster()->millerScore(true, false);
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

	bool hasMillers = Detector::getMaster()->millerCount() > 0;
	std::vector<DetectorPtr> allDetectors;
	Detector::getMaster()->getAllSubDetectors(allDetectors, true);
	bool approximate = FileParser::getKey("GEOMETRY_IS_APPROXIMATE", false);

	if (hasMillers && approximate)
	{
		// we need to search around for the peaks
		addToQueue(allDetectors);
		refineDetectorStrategyWrapper(this, GeometryScoreTypePeakSearch, 0);
	}

	intraPanelCycle();

	for (int i = 0; i < maxCycles; i++)
    {
        cycleNum++;
		interPanelCycle();
    }

	refineBeamCentre();
/*
	logged << "**** Geometry refinement complete ****" << std::endl;
	logged << "**** Now setting METROLOGY_SEARCH_SIZE to 0 for convenience. ****" << std::endl;
	sendLog();

	FileParser::setKey("METROLOGY_SEARCH_SIZE", 0);
 */
}

void GeometryRefiner::printHeader(std::vector<DetectorPtr> detectors, GeometryScoreType type)
{
	std::ostringstream detectorList;

	for (int j = 0; j < detectors.size(); j++)
	{
		if (!detectors[j]->isRefinable(type))
		{
			continue;
		}

		detectorList << detectors[j]->getTag() << " ";
	}

	if (detectorList.str().length())
	{
		logged << "***************************************************" << std::endl;
		logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
		logged << "  Refining detector" << (detectors.size() == 1 ? "" : "s") << ": " << detectorList.str() << std::endl;
		logged << "***************************************************" << std::endl << std::endl;
		sendLog();
	}
}

void GeometryRefiner::intraPanelCycle()
{
	bool hasMillers = Detector::getMaster()->millerCount() > 0;
	std::vector<DetectorPtr> detectors;
	Detector::getMaster()->getAllSubDetectors(detectors, true);

	GeometryScoreType type = GeometryScoreTypeIntrapanel;

	if (hasMillers)
	{
		type = GeometryScoreTypeIntraMiller;
	}

	printHeader(detectors, type);
	addToQueue(detectors);
	refineDetectorStrategyWrapper(this, type, 0);
}

void GeometryRefiner::interPanelCycle()
{
	bool hasMillers = Detector::getMaster()->millerCount() > 0;

	GeometryScoreType type = GeometryScoreTypeInterpanel;

	if (hasMillers)
	{
		type = GeometryScoreTypeInterMiller;
	}

	int i = 64;
	while (i >= 0)
	{
		std::vector<DetectorPtr> detectors;
		detectors = Detector::getMaster()->getSubDetectorsOnLevel(i);
		printHeader(detectors, type);
		i--;

		if (detectors.size() == 0)
		{
			continue;
		}

		addToQueue(detectors);
		refineDetectorStrategyWrapper(this, type, 0);
	}

}

void GeometryRefiner::addToQueue(std::vector<DetectorPtr> dets)
{
	std::lock_guard<std::mutex> lg(queueMutex);

	refineQueue.reserve(refineQueue.size() + dets.size());
	refineQueue.insert(refineQueue.begin(), dets.begin(), dets.end());

	for (int i = 0; i < dets.size(); i++)
	{
		dets[i]->setCycleNum(0);
	}
}

void GeometryRefiner::addToQueue(DetectorPtr det)
{
	std::lock_guard<std::mutex> lg(queueMutex);

	refineQueue.push_back(det);
	det->setCycleNum(det->getCycleNum() + 1);

	logged << "Queue has " << refineQueue.size() << " detectors left." << std::endl;
	sendLog();
}

DetectorPtr GeometryRefiner::getNextDetector()
{
	if (refineQueue.size() == 0)
	{
		return DetectorPtr();
	}

	std::lock_guard<std::mutex> lg(queueMutex);

	if (refineQueue.size() == 0)
	{
		return DetectorPtr();
	}

	DetectorPtr det = refineQueue[0];
	refineQueue.erase(refineQueue.begin());

	return det;
}

void GeometryRefiner::refineDetectorStrategyWrapper(GeometryRefiner *me, GeometryScoreType type, int strategyType)
{
	int maxThreads = FileParser::getMaxThreads();

	boost::thread_group threads;

	for (int i = 0; i < maxThreads; i++)
	{
		boost::thread *thr = new boost::thread(refineDetectorWrapper, me, i, type, strategyType);
		threads.add_thread(thr);
	}

	threads.join_all();

	std::ostringstream logged;
	logged << "Finished a round." << std::endl;
	Logger::log(logged);
	me->reportProgress();

	me->firstCycle = false;
}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, int offset, GeometryScoreType type, int strategyType)
{
	while (true)
	{
		DetectorPtr det = me->getNextDetector();

		if (!det)
		{
			break;
		}

		bool finished = me->refineDetectorStrategy(det, type, strategyType);

		if (!finished)
		{
			me->addToQueue(det);
		}
		else if (det->isRefinable(type))
		{
			me->logged << "Finished detector " << det->getTag() << "!" << std::endl;
			me->sendLog();
		}
	}
}

bool GeometryRefiner::refineDetectorStrategy(DetectorPtr detector, GeometryScoreType type, int strategyType)
{
    std::string typeString = stringForScoreType(type);
    
    if (!detector->isRefinable(type))
    {
        return true;
    }

	_changed = true;

	bool modifiedDetector = false;

	if (type == GeometryScoreTypeInterMiller)
	{
		if (detector->millerCount())
		{
			if (strategyType == 1)
			{
				interPanelGridSearch(detector, type);
			}

			modifiedDetector = interPanelMillerSearch(detector, type);
		}
	}
	else if (type == GeometryScoreTypeIntraMiller)
	{
		modifiedDetector = intraPanelMillerSearch(detector, type);
	}
	else if (type == GeometryScoreTypeInterpanel)
	{
		interPanelGridSearch(detector, type);
		modifiedDetector = interPanelMillerSearch(detector, type);
	}
	else if (type == GeometryScoreTypeIntrapanel)
	{
		modifiedDetector = intraPanelMillerSearch(detector, type);
	}
	else if (type == GeometryScoreTypePeakSearch)
	{
		peakSearchDetector(detector);
	}

	return !modifiedDetector;
}

void GeometryRefiner::peakSearchDetector(DetectorPtr detector)
{
	int searchSize = FileParser::getKey("METROLOGY_SEARCH_SIZE", 3);
	int maxMovement = 50;
	int gridLength = (double)maxMovement / (double)searchSize;

	RefinementGridSearchPtr strategy = makeGridRefiner(detector, GeometryScoreTypePeakSearch);
	strategy->setGridLength(gridLength * 2);

	strategy->addParameter(&*detector, Detector::getAddPixelOffsetX, Detector::setAddPixelOffsetX, searchSize, 0.1);
	strategy->addParameter(&*detector, Detector::getAddPixelOffsetY, Detector::setAddPixelOffsetY, searchSize, 0.1);
	strategy->setJobName(detector->getTag() + "_peaksearch");

	strategy->refine();
}

bool GeometryRefiner::interPanelGridSearch(DetectorPtr detector, GeometryScoreType type)
{
	double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
	detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);
	interNudge /= 2;

	if (detector->getCycleNum() > 0)
	{
		return false;
	}

	RefinementGridSearchPtr strategy = makeGridRefiner(detector, type);
	detector->resetPoke();

	if (type == GeometryScoreTypeInterMiller)
	{
		return false;
		strategy->setJobName(detector->getTag() + "_miller_" + i_to_str(rand() % 10000));
		strategy->setEvaluationFunction(Detector::millerScoreWrapper, &*detector);
		strategy->setFinishFunction(NULL);
		strategy->addParameter(&*detector, Detector::getInterNudgeX, Detector::setInterNudgeX, interNudge / nudgeStep * 3.0, 0, "internudge_x");
		strategy->addParameter(&*detector, Detector::getInterNudgeY, Detector::setInterNudgeY, interNudge / nudgeStep * 3.0, 0, "internudge_y");
		strategy->setGridLength(101);
	}
	else
	{
		double proportion = 1.0 / nudgeStep;
		double onePixRot = proportion * interNudge;
		double totalMovements = 2 * nudgeStep + 1;
		onePixRot /= 4;
		totalMovements *= 4;

		strategy->setJobName(detector->getTag() + "_powder");
		strategy->addParameter(&*detector, Detector::getPokeX, Detector::setPokeX, onePixRot, 0, "poke_x");
		strategy->addParameter(&*detector, Detector::getPokeY, Detector::setPokeY, onePixRot, 0, "poke_y");
		strategy->setGridLength(totalMovements);
		strategy->refine();
		strategy->assignInterpanelMinimum();

	}

	detector->resetPoke();

	return strategy->didChange();
}

bool GeometryRefiner::intraPanelMillerSearch(DetectorPtr detector, GeometryScoreType type)
{
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);

	double totalMovement = nudgeStep;
	int intervals = 2 * totalMovement / 0.25;
	interNudge /= nudgeStep;

	bool changeHappened = false;


    RefinementGridSearchPtr strategy = makeGridRefiner(detector, type);

    detector->prepareInterNudges();
    strategy->setJobName(detector->getTag() + "_miller_stdev_z");
    strategy->addParameter(&*detector, Detector::getNudgeZ, Detector::setNudgeZ, 0.25, 0, "nudge_z");
	strategy->setGridLength(intervals);
    strategy->refine();

	changeHappened = (strategy->didChange() || changeHappened);

	{
		RefinementStrategyPtr strategy = makeRefiner(detector, type);
		strategy->clearParameters();
		detector->prepareInterNudges();
		strategy->addParameter(&*detector, Detector::getNudgeTiltX, Detector::setNudgeTiltX, interNudge, 0, "nudgetilt_x");
		strategy->addParameter(&*detector, Detector::getNudgeTiltY, Detector::setNudgeTiltY, interNudge, 0, "nudgetilt_y");
		strategy->setJobName(detector->getTag() + "_miller_stdev_tilt");
		strategy->refine();
		detector->prepareInterNudges();

		return (strategy->didChange() || changeHappened);
	}
}

bool GeometryRefiner::interPanelMillerSearch(DetectorPtr detector, GeometryScoreType type)
{
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);
    
    RefinementStrategyPtr strategy = makeRefiner(detector, type);
    strategy->setJobName(detector->getTag() + "_miller");
    detector->resetPoke();

	if (type == GeometryScoreTypeInterMiller)
	{
	//	detector->quickJumpToOrigin();

		strategy->addParameter(&*detector, Detector::getInterNudgeX, Detector::setInterNudgeX, interNudge, 0, "internudge_x");
		strategy->addParameter(&*detector, Detector::getInterNudgeY, Detector::setInterNudgeY, interNudge, 0, "internudge_y");
		strategy->addParameter(&*detector, Detector::getInterNudgeZ, Detector::setInterNudgeZ, interNudge, 0, "internudge_z");
	}
	else
	{
		strategy->setJobName(detector->getTag() + "_powder");
		strategy->addParameter(&*detector, Detector::getPokeX, Detector::setPokeX, interNudge, 0, "internudge_x");
		strategy->addParameter(&*detector, Detector::getPokeY, Detector::setPokeY, interNudge, 0, "internudge_y");
		strategy->addParameter(&*detector, Detector::getPokeZ, Detector::setPokeZ, interNudge, 0, "internudge_z");
	}

	_changed = true;
    strategy->refine();
    
    detector->resetPoke();

	return strategy->didChange();
}

void GeometryRefiner::refineBeamCentre()
{
    DetectorPtr detector = Detector::getMaster();

    if (detector->millerCount() > 0)
    {
        return;
    }
    
    logged << "***************************************************" << std::endl;
    logged << "  Cycle " << cycleNum << ", event " << refinementEvent << std::endl;
    logged << "  Refining beam centre: " << detector->getTag() << std::endl;
    logged << "***************************************************" << std::endl << std::endl;
    sendLog();
    
    RefinementStrategyPtr gridSearch = makeRefiner(detector, GeometryScoreTypeBeamCentre);
    gridSearch->setJobName("grid_search_beam_centre_" + i_to_str(refinementEvent));
    gridSearch->setCycles(99);
    
    double nudgeStep, nudgeTiltX, nudgeTiltY, interNudge;
    detector->nudgeTiltAndStep(&nudgeTiltX, &nudgeTiltY, &nudgeStep, &interNudge);

    gridSearch->addParameter(&*detector, Detector::getInterNudgeX, Detector::setInterNudgeX, interNudge, 0.001, "poke_x");
    gridSearch->addParameter(&*detector, Detector::getInterNudgeY, Detector::setInterNudgeY, interNudge, 0.001, "poke_y");

    gridSearch->refine();
	_changed = true;
    
    reportProgress();
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
    aManager->setActiveDetector(detector, GeometryScoreTypeIntrapanel);
    aManager->setPseudoScoreType(PseudoScoreTypeIntraPanel);
    strategy->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    
    strategy->refine();
}

bool vecToSpotRatioLessThan(ImagePtr one, ImagePtr two)
{
	double ratio1 = one->spotVectorCount() / (one->spotCount() * one->spotCount());
	double ratio2 = two->spotVectorCount() / (one->spotCount() * one->spotCount());

	return (ratio1 < ratio2);
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
	int cherryPick = FileParser::getKey("CHERRY_PICK", 0);

    if (cherryPick > 0)
	{
		for (int i = 0; i < newImages.size(); i++)
		{
			if (!newImages[i]->acceptableSpotCount())
			{
				newImages.erase(newImages.begin() + i);
				i--;
			}
		}

		if (cherryPick > newImages.size())
		{
			cherryPick = (int)newImages.size();
		}

		logged << "Cherry picking " << cherryPick << " images from " << newImages.size() << std::endl;
		std::sort(newImages.begin(), newImages.end(), vecToSpotRatioLessThan);
		images.clear();
		images.reserve(cherryPick);
		images.insert(images.begin(), newImages.begin(), newImages.begin() + cherryPick);
	}
	else
	{
		images = newImages;
	}

    manager = IndexManagerPtr(new IndexManager(images));

	if (cherryPick > 0)
	{
		logged << "Cherry picked powder pattern..." << std::endl;
		sendLog();
		manager->powderPattern();
	}
}
