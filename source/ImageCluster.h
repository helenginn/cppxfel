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
#include <boost/tuple/tuple.hpp>
#include <string>
#include <iostream>

class Image;

class ImageCluster
{
private:
    std::vector<Image *>images;
    std::map<Image *, std::map<Image *, MatrixPtr> > imageToImageMatrices;
    std::vector<ImageClusterPtr> clusters;
    std::map<Image *, std::vector<std::vector<Image *> > >goodPermutations;
    ImageCluster *parentCluster;
    std::vector<std::pair<SpotPtr, SpotPtr> > commonSpots;
    void findSeedingClusters();
    int unexpectedMatches;
    double lastScore;
    int lastImageUse;
    int trust;
    std::vector<std::vector<double> > wobbles;
    
    bool flaggedForMerging;
    double distanceBetweenSpotPair(std::pair<SpotPtr, SpotPtr> pair, bool addWobble);
    static bool clusterMoreUnexpected(ImageClusterPtr cluster1, ImageClusterPtr cluster2);
    static bool clusterMoreTrustworthy(ImageClusterPtr cluster1, ImageClusterPtr cluster2);
    
    bool containsImage(Image *image);
    std::ostringstream logged;
    void sendLog(LogLevel priority = LogLevelNormal);
    void setSeed();
    void growSeeds();
    void clusterStatistics();
    void expandMatrices(Image *newImage);
    void merge(ImageClusterPtr subA, ImageClusterPtr subB, std::pair<Image *, Image *> sharedImages);
    void calculateMatrices(Image *one, Image *two);
    void changeMatrix(Image *firstImage, Image *secondImage, MatrixPtr newMatrix, bool setInverse = true);
    void fillMatrixGaps();
    void findMergableClusters(int minUnexpected = 0);
    Image *sharesCommonImagesWithCluster(ImageClusterPtr otherCluster);
    void copyMatricesToParent();
    void writeRotations(std::string filename);
    void setSharedSpots(Image *one, Image *two);
    void findSpotPairs();
    
    int uniqueUnexpectedMatches();
    bool wobbled;
    void applyWobbles();
    void resetWobbles();
    void randomiseImages();
    void shakeUp();
    void minimizeWobble();
    void addTrust();
    static double scoreWrapper(void *object);
    bool nextPermutation(std::vector<Image *> **testImagePointer, Image *image, int count);
    double selfConsistencyScore(bool addWobble = false, int maxSuccesses = 1, int maxTrials = 5000000);
    void breakUp(ImageClusterPtr toBreak, bool recursive = false);
    void breakUp(bool recursive = false);
    void setParentCluster(ImageCluster *newParent);
    MatrixPtr getImageToImageMatrix(Image *one, Image *two, bool addWobble = false);
public:
    ImageCluster(std::vector<Image *>theImages, ImageCluster *parent = NULL);
    void process();
    void imageStatistics();
    void writeSpots(std::string filename);
    
};

#endif /* defined(__cppxfel__ImageCluster__) */
