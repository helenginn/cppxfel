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
    void refineDetector(DetectorPtr detector);
    static void refineDetectorWrapper(GeometryRefiner *me, std::vector<DetectorPtr> detectors, int offset);
    
public:
    GeometryRefiner();
    
    void refineGeometry();
    
    void setImages(std::vector<ImagePtr> newImages);
    
};

#endif /* defined(__cppxfel__GeometryRefiner__) */
