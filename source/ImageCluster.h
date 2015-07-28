//
//  ImageCluster.h
//  cppxfel
//
//  Created by Helen Ginn on 28/07/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__ImageCluster__
#define __cppxfel__ImageCluster__

#include <stdio.h>
#include "parameters.h"
#include "Logger.h"

class Image;

class ImageCluster
{
private:
    std::vector<Image *>images;
    std::map<Image *, std::map<Image *, MatrixPtr> > imageToImageMatrices;
    std::vector<ImageClusterPtr> clusters;
    void findSeedingClusters();
    
    std::ostringstream logged;
    void sendLog(LogLevel priority = LogLevelNormal);
    void setSeed();
public:
    ImageCluster(std::vector<Image *>theImages);
    void process();
    void imageStatistics();
    
};

#endif /* defined(__cppxfel__ImageCluster__) */
