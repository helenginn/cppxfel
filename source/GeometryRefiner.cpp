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
#include "GetterSetterMap.h"
#include "IndexManager.h"
#include "misc.h"

GeometryRefiner::GeometryRefiner()
{
    refinementEvent = 0;
    cycleNum = 0;
}

void GeometryRefiner::refineGeometry()
{
    for (int i = 0; i < 3; i++)
    {
        refineGeometryCycle();
    }
}

void GeometryRefiner::refineGeometryCycle()
{
    std::vector<DetectorPtr> sameLevelDetectors;
    std::vector<DetectorPtr> nextLevelDetectors;
    
    DetectorPtr master = Detector::getMaster();
    nextLevelDetectors.push_back(master);
    
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
        
        refinementEvent++;
        
        manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv", false);
    }
    
    GeometryParser geomParser = GeometryParser("test.cppxfel_geom", GeometryFormatCppxfel);
    geomParser.writeToFile("new.cppxfel_geom");
}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < detectors.size(); i += maxThreads)
    {
        me->refineDetector(detectors[i]);
    }
}

void GeometryRefiner::refineDetector(DetectorPtr detector)
{
    Logger::log(logged);

    GetterSetterMapPtr beamRefinementMap = GetterSetterMapPtr(new GetterSetterMap());
    
    double zMovement = detector->isLUCA() ? 2.0 : 0.4;
    
    beamRefinementMap->setJobName("Detector " + detector->getTag());
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointX,
                                    Detector::setArrangedMidPointX, 2.0, 0.02, "midPointX");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointY,
                                    Detector::setArrangedMidPointY, 2.0, 0.02, "midPointY");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                                    Detector::setArrangedMidPointZ, 0.4, 0.02, "midPointZ");
    
    beamRefinementMap->setVerbose(true);
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    aManager->setActiveDetector(detector);
    
    beamRefinementMap->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    
    beamRefinementMap->setCycles(15);
    
    beamRefinementMap->refine(GetterSetterStepSearch);
    
    beamRefinementMap->clearParameters();
    
    beamRefinementMap->addParameter(&*detector, Detector::getAlpha, Detector::setAlpha, 0.005, 0.00001, "alpha");
    beamRefinementMap->addParameter(&*detector, Detector::getBeta, Detector::setBeta, 0.005, 0.00001, "beta");
    
    if (!detector->isLUCA())
    {
        beamRefinementMap->addParameter(&*detector, Detector::getGamma, Detector::setGamma, 0.005, 0.00001, "gamma");
    }
    
    zMovement = detector->isLUCA() ? 0.5 : 0.1;
    
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointX,
                                    Detector::setArrangedMidPointX, 0.2, 0.01, "midPointX");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointY,
                                    Detector::setArrangedMidPointY, 0.2, 0.01, "midPointY");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                                    Detector::setArrangedMidPointZ, zMovement, 0.01, "midPointZ");
    
    beamRefinementMap->refine(GetterSetterStepSearch);
    
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));
    manager->powderPattern("geom_refinement_event_0_start.csv", false);
    
}