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

typedef enum
{
    GeometryScoreTypeMiller,
    GeometryScoreTypeMillerStdev,
    GeometryScoreTypeIntrapanel,
    GeometryScoreTypeInterpanel,
    GeometryScoreTypeInterAngle,
    GeometryScoreTypeIntraAngle,
    GeometryScoreTypeAngleConsistency,
    GeometryScoreTypeBeamCentre,
    
} GeometryScoreType;

class GeometryRefiner : public LoggableObject
{
private:
    std::vector<ImagePtr> images;
    IndexManagerPtr manager;
    int refinementEvent;
    int cycleNum;
    void refineDetectorStrategy(DetectorPtr detector, int strategy);
    static void refineDetectorStrategyWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int strategy);
    static void refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset, int strategy);
    void refineGeometryCycle();
    void geometryCycleForDetector(std::vector<DetectorPtr> detectors);
    double lastInterScore;
    double lastIntraScore;
    double lastIntraAngleScore;
    double lastInterAngleScore;
    void refineMasterDetector();
    void refineBeamCentre();
    
    RefinementStrategyPtr makeRefiner(DetectorPtr detector, GeometryScoreType type);
    void gridSearch(DetectorPtr detector, double start, double end);
    void refineDetector(DetectorPtr detector, GeometryScoreType type);
    
public:
    GeometryRefiner();
    
    void refineGeometry();
    void reportProgress();
    void refineUnitCell();
    
    void setImages(std::vector<ImagePtr> newImages);
    
};

#endif /* defined(__cppxfel__GeometryRefiner__) */
