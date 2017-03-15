//
//  GeometryRefiner.h
//  cppxfel
//
//  Created by Helen Ginn on 28/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__GeometryRefiner__
#define __cppxfel__GeometryRefiner__

#include "LoggableObject.h"
#include "parameters.h"
#include <stdio.h>

class GeometryRefiner : public LoggableObject
{
private:
    std::vector<ImagePtr> images;
    IndexManagerPtr manager;
    int refinementEvent;
    int cycleNum;
    bool firstCycle;
    void refineDetectorStrategy(DetectorPtr detector, GeometryScoreType type, int strategyType);
    static void refineDetectorStrategyWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, GeometryScoreType type, int strategyType);
    static void refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors,
                                      int offset, GeometryScoreType type, int strategyType);
    void refineGeometryCycle();
    bool geometryCycleForDetector(std::vector<DetectorPtr> detectors, bool interPanelOnly);
    double lastInterScore;
    double lastIntraScore;
    double lastIntraAngleScore;
    double lastInterAngleScore;
    
    void refineBeamCentre();

    void printHeader(std::vector<DetectorPtr> detectors, std::string detectorList);
    
    void interPanelGridSearch(DetectorPtr detector);
    void interPanelNormalSearch(DetectorPtr detector);
    void interPanelMillerSearch(DetectorPtr detector);
    void intraPanelMillerSearch(DetectorPtr detector);
    void intraPanelNormalSearch(DetectorPtr detector);
    void nudgeZGridSearch(DetectorPtr detector);    
    
    RefinementGridSearchPtr makeGridRefiner(DetectorPtr detector, GeometryScoreType type);
    RefinementStrategyPtr makeRefiner(DetectorPtr detector, GeometryScoreType type);
    void gridSearchDetectorDistance(DetectorPtr detector, double start, double end);
    void refineDetector(DetectorPtr detector, GeometryScoreType type);
    
public:
    GeometryRefiner();
    
    void refineGeometry();
    void reportProgress();
    void refineUnitCell();
    void startingGraphs();
    
    void setImages(std::vector<ImagePtr> newImages);
    
};

#endif /* defined(__cppxfel__GeometryRefiner__) */
