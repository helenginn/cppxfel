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
}

void GeometryRefiner::refineGeometry()
{
    std::vector<DetectorPtr> sameLevelDetectors;
    std::vector<DetectorPtr> nextLevelDetectors;
    
    DetectorPtr master = Detector::getMaster();
    nextLevelDetectors.push_back(master);
    
    while (nextLevelDetectors.size())
    {
        sameLevelDetectors = nextLevelDetectors;
        nextLevelDetectors.clear();
        /*
        double maxThreads = FileParser::getMaxThreads();
        boost::thread_group threads;
        
        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(refineDetectorWrapper, this, sameLevelDetectors, i);
            threads.add_thread(thr);
        }
        
        threads.join_all();
        */
        
        
        for (int i = 0; i < sameLevelDetectors.size(); i++)
        {
            refineDetector(sameLevelDetectors[i]);
            
            for (int j = 0; j < sameLevelDetectors[i]->childrenCount(); j++)
            {
                nextLevelDetectors.push_back(sameLevelDetectors[i]->getChild(j));
            }
        }

        refinementEvent++;
        
        manager->powderPattern("geom_refinement_event_" + i_to_str(refinementEvent) + ".csv");

    }
    
    GeometryParser geomParser = GeometryParser("test.cppxfel_geom", GeometryFormatCppxfel);
    geomParser.writeToFile("new.cppxfel_geom");
}

void GeometryRefiner::refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset)
{
    int maxThreads = 1;//FileParser::getMaxThreads();
    
    for (int i = offset; i < detectors.size(); i += maxThreads)
    {
        me->refineDetector(detectors[i]);
    }
}

void GeometryRefiner::refineDetector(DetectorPtr detector)
{
    logged << "*********************************" << std::endl;
    logged << "  Refining detector: " << detector->getTag() << std::endl;
    logged << "*********************************" << std::endl << std::endl;
    
    Logger::log(logged);

    GetterSetterMapPtr beamRefinementMap = GetterSetterMapPtr(new GetterSetterMap());
    
    beamRefinementMap->setJobName("Detector " + detector->getTag());
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointX,
                                    Detector::setArrangedMidPointX, 2.0, 0.02, "midPointX");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointY,
                                    Detector::setArrangedMidPointY, 2.0, 0.02, "midPointY");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                                    Detector::setArrangedMidPointZ, 2.0, 0.02, "midPointZ");
    
    beamRefinementMap->setVerbose(true);
    
    IndexManagerPtr aManager = IndexManagerPtr(new IndexManager(images));
    
    beamRefinementMap->setEvaluationFunction(IndexManager::pseudoScore, &*aManager);
    
    beamRefinementMap->setCycles(15);
    
    beamRefinementMap->refine(GetterSetterStepSearch);
    
    beamRefinementMap->clearParameters();
    
    beamRefinementMap->addParameter(&*detector, Detector::getAlpha, Detector::setAlpha, 0.005, 0.0002, "alpha");
    beamRefinementMap->addParameter(&*detector, Detector::getBeta, Detector::setBeta, 0.005, 0.0002, "beta");
    
    if (!detector->isLUCA())
    {
        beamRefinementMap->addParameter(&*detector, Detector::getGamma, Detector::setGamma, 0.005, 0.0002, "gamma");
    }
    
    double zMovement = detector->isLUCA() ? 0.5 : 0.1;
    
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointX,
                                    Detector::setArrangedMidPointX, 0.2, 0.1, "midPointX");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointY,
                                    Detector::setArrangedMidPointY, 0.2, 0.1, "midPointY");
    beamRefinementMap->addParameter(&*detector, Detector::getArrangedMidPointZ,
                                    Detector::setArrangedMidPointZ, zMovement, 0.1, "midPointZ");
    
    beamRefinementMap->refine(GetterSetterStepSearch);
    
}

void GeometryRefiner::setImages(std::vector<ImagePtr> newImages)
{
    images = newImages;
    manager = IndexManagerPtr(new IndexManager(images));
    manager->powderPattern("geom_refinement_event_0_start.csv");
    
}