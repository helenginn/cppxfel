//
//  ImageCluster.cpp
//  cppxfel
//
//  Created by Helen Ginn on 28/07/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "ImageCluster.h"
#include "CommonCircle.h"
#include "Image.h"
#include "Logger.h"
#include "Vector.h"
#include "FileParser.h"
#include <fstream>
#include "misc.h"
#include <iomanip>

void ImageCluster::resetWobbles()
{
    for (int i = 0; i < images.size(); i++)
    {
        wobbles[i] = std::vector<double>(3, 0);
    }
}

void ImageCluster::applyWobbles()
{
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images.size(); j++)
        {
            if (i == j) continue;
            
            double alpha = wobbles[j][0];
            double beta = wobbles[j][1];
            double gamma = wobbles[j][2];
            
            double alphaOne = wobbles[i][0];
            double betaOne = wobbles[i][1];
            double gammaOne = wobbles[i][2];
            
            imageToImageMatrices[images[i]][images[j]]->rotate(-alphaOne, -betaOne, -gammaOne);
            imageToImageMatrices[images[i]][images[j]]->rotate(alpha, beta, gamma);
        }
    }

    wobbled = true;
    resetWobbles();
}

void ImageCluster::randomiseImages()
{
    std::random_shuffle(images.begin(), images.end());
    resetWobbles();
}

ImageCluster::ImageCluster(std::vector<Image *> theImages, ImageCluster *parent)
{
    images = theImages;
    wobbles.resize(images.size());
    
    unexpectedMatches = 0;
    lastScore = 0;
    lastImageUse = 0;
    flaggedForMerging = false;
    trust = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images.size(); j++)
        {
            imageToImageMatrices[images[i]][images[j]] = MatrixPtr(new Matrix());
        }
    }
    
    assert(imageToImageMatrices.size() == images.size());
    
    wobbled = false;
    parentCluster = parent;
}

void ImageCluster::expandMatrices(Image *newImage)
{
    for (int i = 0; i < images.size(); i++)
    {
        imageToImageMatrices[images[i]][newImage] = MatrixPtr(new Matrix());
        imageToImageMatrices[newImage][images[i]] = MatrixPtr(new Matrix());
    }
    
    wobbles.resize(images.size());
    
    assert(imageToImageMatrices.size() == images.size());
}

void ImageCluster::copyMatricesToParent()
{
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = i + 1; j < images.size(); j++)
        {
            if (i == j)
                continue;
            
            Image *imgI = images[i];
            Image *imgJ = images[j];
            
            if (!imageToImageMatrices[imgI][imgJ]->isIdentity())
            {
                parentCluster->changeMatrix(imgI, imgJ, imageToImageMatrices[imgI][imgJ]);
            }
        }
    }
}

void ImageCluster::changeMatrix(Image *firstImage, Image *secondImage, MatrixPtr newMatrix, bool setInverse)
{
    imageToImageMatrices[firstImage][secondImage] = newMatrix;
    
    if (setInverse)
        imageToImageMatrices[secondImage][firstImage] = newMatrix->inverse3DMatrix();
    
    if (parentCluster)
    {
        parentCluster->changeMatrix(firstImage, secondImage, newMatrix, setInverse);
    }
}

void ImageCluster::fillMatrixGaps()
{
    int totalSuccesses = 0;
    int successes = -1;
    int rounds = 0;
    
    logged << "Filling gaps in matrix map" << std::endl;
    sendLog(LogLevelDetailed);
    
    while (successes != 0)
    {
        successes = 0;
        
        for (int i = 0; i < images.size(); i++)
        {
            for (int j = 0; j < images.size(); j++)
            {
                if (i == j)
                    continue;
                
                MatrixPtr newMat = MatrixPtr(new Matrix());
                
                Image *imgI = images[i];
                Image *imgJ = images[j];
                
                if (imageToImageMatrices[imgI][imgJ]->isIdentity())
                {
                    // find another matrix which has data with both I and J
                    bool found = false;
                    
                    for (int k = 0; k < images.size() && !found; k++)
                    {
                        Image *imgK = images[k];
                        
                        if (imgK == imgI || imgK == imgJ)
                            continue;
                        
                        if (images[k]->isUnexpectedMatch(this))
                            continue;
                        
                        MatrixPtr mapIOntoK = imageToImageMatrices[imgI][imgK]->copy();
                        if (mapIOntoK->isIdentity())
                            continue;
                        MatrixPtr mapKOntoJ = imageToImageMatrices[imgK][imgJ]->copy();
                        if (mapKOntoJ->isIdentity())
                            continue;
                        
                        MatrixPtr mapKOntoI = imageToImageMatrices[imgK][imgI]->copy();
                        if (mapIOntoK->isIdentity())
                            continue;
                        MatrixPtr mapJOntoK = imageToImageMatrices[imgJ][imgK]->copy();
                        if (mapKOntoJ->isIdentity())
                            continue;
                        
                        MatrixPtr mapIOntoJ = mapIOntoK->copy();
                        mapIOntoJ->multiply(*mapKOntoJ);
                        
                        if (!mapIOntoJ->isIdentity())
                        {
                            logged << "Filled matrix " << i << ", " << j << " using " << k << std::endl;
                            sendLog(LogLevelDebug);
                            successes++;
                        }
                        
                        found = true;
                        changeMatrix(imgI, imgJ, mapIOntoJ, false);
                    }
                    
                    if (found == false)
                    {
                 //       std::cout << "!" << std::endl;
                    }
                }
            }
        }
        rounds++;
        
        totalSuccesses += successes;
        logged << "Finished round " << rounds << ", changed " << successes << " image pair matrices." << std::endl;
        sendLog(LogLevelDebug);
        
        if (rounds > 10 && false)
            break;
    }
    
    int identities = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images.size(); j++)
        {
            Image *imgI = images[i];
            Image *imgJ = images[j];
            
            if (imageToImageMatrices[imgI][imgJ]->isIdentity())
            {
                identities++;
            }
        }
    }
    
    logged << "Filled in " << totalSuccesses << " matrices for " << images.size() << " images." << std::endl;
    sendLog(LogLevelDebug);
}

void ImageCluster::calculateMatrices(Image *one, Image *two)
{
    CommonCirclePair circlePair = one->getCommonCircleWithImage(two);
    
    double angleSeparation1 = circlePair.first->getAngleSeparation();
    double angleSeparation2 = circlePair.second->getAngleSeparation();
    double angleSeparation = (angleSeparation1 + angleSeparation2) / 2;
    
    double swivelA = circlePair.first->swivelToHorizontalAxis(false);
    double swivelB = circlePair.second->swivelToHorizontalAxis(true);
    

    logged << "Angle separation: " << angleSeparation << ", swivelA: " << swivelA << ", swivelB: " << swivelB << std::endl;
    
    vec axis = new_vector(0, 0, -1);
    MatrixPtr rotB = MatrixPtr(new Matrix());
    rotB->rotate(0, -angleSeparation, 0);
    rotB->multiplyVector(&axis);
    rotB->rotateRoundUnitVector(axis, swivelB);
    
    MatrixPtr rotA = MatrixPtr(new Matrix());
    rotA->rotate(0, 0, swivelA);
    
    logged << "axis: " << axis.h << " " << axis.k << " " << axis.l << std::endl;
    
    MatrixPtr invRotA = rotA->inverse3DMatrix();
    MatrixPtr rotInvAB = rotB->copy();
    rotInvAB->preMultiply(*invRotA);
    rotB->rotate(0, 0, swivelA);
    
    changeMatrix(one, two, rotB);
    changeMatrix(two, one, rotB->inverse3DMatrix());
    
    sendLog(LogLevelDebug);
}

void ImageCluster::findSpotPairs()
{
    commonSpots = std::vector<std::pair<SpotPtr, SpotPtr> >();
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images.size(); j++)
        {
            if (i == j)
                continue;
            
            CommonCirclePair circle = images[i]->getCommonCircleWithImage(images[j]);
            
            if (!(circle.first && circle.second))
                continue;
            
            std::vector<std::pair<SpotPtr, SpotPtr> > orderedSpots;
            circle.first->orderSpotsInPairs(circle.second, &orderedSpots);
            
            for (int i = 0; i < orderedSpots.size(); i++)
            {
                commonSpots.push_back(orderedSpots[i]);
            }
        }
    }
}

void ImageCluster::setSharedSpots(Image *one, Image *two)
{
    int shared = 0;
    
    CommonCirclePair circlePair = one->getCommonCircleWithImage(two);
    
    shared += circlePair.first->makeSuccessful();
    shared += circlePair.second->makeSuccessful();
    
    logged << "Found " << shared << " shared spots." << std::endl;
    sendLog(LogLevelDetailed);

}

void ImageCluster::setSeed()
{
    assert(images.size() == 2);
    
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->setSeeded();
    }
    
    setSharedSpots(images[0], images[1]);
    calculateMatrices(images[0], images[1]);
}

void ImageCluster::process()
{
    double nextPercentage = 5;
    int step = 5;
    unsigned long totalImages = images.size();
    
    logged << "N: Finding common circles with images." << std::endl << "0%... ";
    sendLog();
    
    for (int i = 0; i < totalImages; i++)
    {
        images[i]->commonCirclesWithImages(images);
        
        double percentage = (double)i / (double)totalImages * 100;
        
        if (percentage > nextPercentage)
        {
            logged << nextPercentage << "%... ";
            sendLog();
            nextPercentage += step;
        }
    }
    
    std::cout << std::endl;
    
    logged << "N: Determined common circles between images." << std::endl;
    sendLog();
    
    imageStatistics();
    
    int cycles = FileParser::getKey("CLUSTER_CYCLES", 3);
    
    for (int i = 0; i <= cycles; i++)
    {
        logged << "****** BEGINNING MACROCYCLE " << i << " ******" << std::endl;
        sendLog();
        
        findSeedingClusters();
    
        unsigned long clusterCount = clusters.size();
        unsigned long difference = -1;
        
        while (difference != 0)
        {
            clusterCount = clusters.size();
            std::random_shuffle(clusters.begin(), clusters.end());
            growSeeds();
            findMergableClusters(10);
            findMergableClusters();
            
            logged << "Finished with " << clusters.size() << " clusters." << std::endl;
            sendLog(LogLevelNormal);
            
            difference = clusterCount - clusters.size();
        }
        
        if (difference == 0)
        {
            shakeUp();
            clusterStatistics();
        }
    }
    
    for (int i = 0; i < clusters.size(); i++)
    {
        if (clusters[i]->images.size() > 2)
        {
            clusters[i]->writeRotations("cluster_" + i_to_str(i) + ".rot");
            clusters[i]->writeSpots("cluster_" + i_to_str(i) + ".pdb");
        }
    }
}

bool ImageCluster::clusterMoreTrustworthy(ImageClusterPtr cluster1, ImageClusterPtr cluster2)
{
    return (cluster1->trust > cluster2->trust);
}

bool ImageCluster::clusterMoreUnexpected(ImageClusterPtr cluster1, ImageClusterPtr cluster2)
{
    return (cluster1->images.size() > cluster2->images.size());
}

void ImageCluster::clusterStatistics()
{
    logged << "****** CLUSTER STATISTICS ******" << std::endl;
    sendLog();
    
    std::sort(clusters.begin(), clusters.end(), clusterMoreUnexpected);
    
    for (int i = 0; i < clusters.size(); i++)
    {
        clusters[i]->findSpotPairs();
        
        logged << "Cluster with " << clusters[i]->images.size() << " images and " << commonSpots.size() << " common spots has self-consistency score of " << clusters[i]->lastScore << " (" << lastImageUse << ")" << std::endl;
        sendLog(LogLevelNormal);
    }
    
    std::random_shuffle(clusters.begin(), clusters.end());
}

void ImageCluster::shakeUp()
{
    logged << "****** SHAKING IT UP ******" << std::endl;
    sendLog();
    
    double threshold = FileParser::getKey("CLUSTER_MERGE_THRESHOLD", 0.05);
    
    std::sort(clusters.begin(), clusters.end(), clusterMoreTrustworthy);
    
    for (int i = 0; i < clusters.size(); i++)
    {
        clusters[i]->minimizeWobble();
        
        if (clusters[i]->lastScore > threshold)
            breakUp(clusters[i]);
    }
    
    std::random_shuffle(clusters.begin(), clusters.end());
    randomiseImages();
}

void ImageCluster::findSeedingClusters()
{
    logged << "Starting to find seed clusters of three images." << std::endl;
    sendLog();
    
    int clusterCount = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        Image *oneImage = images[i];
        
        for (int j = i + 1; j < images.size(); j++)
        {
            Image *twoImage = images[j];
            
            if (oneImage->hasSeeded() || twoImage->hasSeeded())
                continue;
            
            if (oneImage->hasCommonCircleWithImage(twoImage))
            {
                std::vector<Image *>twoImages;
                twoImages.push_back(oneImage);
                twoImages.push_back(twoImage);
                
                ImageClusterPtr newCluster = ImageClusterPtr(new ImageCluster(twoImages, this));
                newCluster->setSeed();
                clusters.push_back(newCluster);
                
                clusterCount++;
            }
        }
    }
    
    logged << "N: Found " << clusterCount << " two-image seeds." << std::endl;
    sendLog();
    
    /*for (int i = 0; i < clusters.size(); i++)
    {
        clusters[i]->writeRotations("cluster_" + i_to_str(i) + ".rot");
        clusters[i]->writeSpots("cluster_" + i_to_str(i) + ".pdb");
    }*/
}

void ImageCluster::imageStatistics()
{
    unsigned long totalImages = images.size();
    int hopefulImages = 0;
    int commonCirclesTotal = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        unsigned long pairCount = images[i]->circlePairCount();
        
        if (pairCount > 0)
        {
            hopefulImages++;
            commonCirclesTotal += pairCount;
            
            images[i]->circlePair(0).first->printSpots();
            images[i]->circlePair(0).second->printSpots();
            logged << "***" << std::endl;
            sendLog();
        }
    }
    
    double commonCirclesPerHopefulImage = (double)commonCirclesTotal / (double)hopefulImages;
    
    logged << "N: Total images: " << totalImages << std::endl;
    logged << "N: Hopeful images (common circles with at least one other image): " << hopefulImages << std::endl;
    logged << "N: Common circles per hopeful image: " << commonCirclesPerHopefulImage << std::endl;

    sendLog();
}

bool ImageCluster::containsImage(Image *image)
{
    for (int i = 0; i < images.size(); i++)
    {
        if (images[i] == image)
            return true;
    }
    
    return false;
}

void ImageCluster::growSeeds()
{
    for (int i = 0; i < clusters.size(); i++)
    {
        bool added = false;
        
        for (int j = 0; j < images.size(); j++)
        {
            bool inSeed = images[j]->hasSeeded();
            
            if (inSeed)
                continue;
            
            for (int k = 0; k < clusters[i]->images.size(); k++)
            {
                if (clusters[i]->containsImage(images[j]))
                    continue;

                Image *test = clusters[i]->images[k];
                
                if (images[j]->hasCommonCircleWithImage(test))
                {
                    clusters[i]->images.push_back(images[j]);
                    images[j]->setSeeded();
                    
                    clusters[i]->setSharedSpots(test, images[j]);
                    clusters[i]->expandMatrices(images[j]);
                    clusters[i]->calculateMatrices(test, images[j]);
                    clusters[i]->fillMatrixGaps();
                    
                    added = true;
                }
            }
        }
        
        if (added)
        {
            double score = clusters[i]->selfConsistencyScore();
            logged << "Grown seed to " << clusters[i]->images.size() << " images, self-consistency of " << score << std::endl;
            sendLog(LogLevelDetailed);
        }
    }
}

void ImageCluster::findMergableClusters(int minUnexpected)
{
    logged << "****** BEGINNING CLUSTER MERGE ******" << std::endl;
    logged << "*** MINIMUM " << minUnexpected << " UNEXPECTED MATCHES *** " << std::endl;
    sendLog();
    
    int total = 0;
    std::map<std::pair<ImageClusterPtr, ImageClusterPtr>, std::pair<Image *, Image *>> clusterPairs;
    
    for (int i = 0; i < clusters.size(); i++)
    {
        ImageClusterPtr aCluster = clusters[i];
        
        if (aCluster->images.size() > minUnexpected)
            continue;
        
        ImageClusterPtr chosenCluster = ImageClusterPtr();
        std::pair<Image *, Image *> chosenPair;
        
        for (int j = i + 1; j < clusters.size() && !chosenCluster; j++)
        {
            ImageClusterPtr bCluster = clusters[j];
            
            if (bCluster->images.size() > minUnexpected)
                continue;

            for (int k = 0; k < clusters[i]->images.size() && !chosenCluster; k++)
            {
                for (int l = 0; l < clusters[j]->images.size() && !chosenCluster; l++)
                {
                    if (aCluster->images[k]->hasCommonCircleWithImage(bCluster->images[l]))
                    {
                        if ((aCluster->flaggedForMerging || bCluster->flaggedForMerging))
                            continue;
                        
                        aCluster->flaggedForMerging = true;
                        bCluster->flaggedForMerging = true;
                        chosenCluster = bCluster;
                        chosenPair = std::make_pair(aCluster->images[k], bCluster->images[l]);
                    }
                }
            }
        }
        
        if (chosenCluster)
        {
            clusterPairs[std::make_pair(aCluster, chosenCluster)] = chosenPair;
            total++;
        }
    }
    
    logged << "Found " << total << " clusters to be merged." << std::endl;
    sendLog();
    
    for (std::map<std::pair<ImageClusterPtr, ImageClusterPtr>, std::pair<Image *, Image *>>::iterator i = clusterPairs.begin();
         i != clusterPairs.end(); i++)
    {
        ImageClusterPtr aCluster = i->first.first;
        ImageClusterPtr bCluster = i->first.second;
        std::pair<Image *, Image *> chosenImages = clusterPairs[i->first];
        
        if (aCluster && bCluster)
            merge(aCluster, bCluster, chosenImages);
    }
}

void ImageCluster::minimizeWobble()
{
    resetWobbles();
    std::vector<std::vector<double> > wobbleSteps(images.size());
    double beforeScore = selfConsistencyScore(false);
    
    double averageStep = 0.1;
    double latestScore = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        wobbleSteps[i] = std::vector<double>(3, 0.1);
    }
    
    logged << "Score before wobble: " << beforeScore << std::endl;
    sendLog(LogLevelDetailed);
    
    int cycles = 0;
    
    while (averageStep > 0.01 && cycles < 50)
    {
        for (int i = 0; i < wobbles.size(); i++)
        {
            minimizeParameter(wobbleSteps[i][0], wobbles[i][0], scoreWrapper, this);
            minimizeParameter(wobbleSteps[i][1], wobbles[i][1], scoreWrapper, this);
            minimizeParameter(wobbleSteps[i][2], wobbles[i][2], scoreWrapper, this);
            
            logged << "Corrected wobble by " << wobbles[i][0] << ", " << wobbles[i][1] << ", " << wobbles[i][2] << std::endl;
            sendLog(LogLevelDebug);
        }
        
        latestScore = selfConsistencyScore(true);
        logged << "Wobble cycle score: " << latestScore << std::endl;
        sendLog(LogLevelDetailed);
        
        for (int i = 0; i < wobbleSteps.size(); i++)
        {
            for (int j = 0; j < 3; j++)
                averageStep += wobbleSteps[i][j];
        }
        
        averageStep /= wobbleSteps.size() * 3;
        cycles++;
    }
    
    logged << "After wobbling, consistency went from " << beforeScore << " to " << latestScore << std::endl;
    sendLog();
    
    applyWobbles();
}

MatrixPtr ImageCluster::getImageToImageMatrix(Image *one, Image *two, bool addWobble)
{
    MatrixPtr mat = imageToImageMatrices[one][two]->copy();
    
    int imageTwo = -1;
    int imageOne = -1;
    
    if (addWobble)
    {
        for (int i = 0; i < images.size(); i++)
        {
            if (images[i] == two)
            {
                imageTwo = i;
            }
            
            if (images[i] == one)
            {
                imageOne = i;
            }
        }
        
        {
            double alphaOne = wobbles[imageOne][0];
            double betaOne = wobbles[imageOne][1];
            double gammaOne = wobbles[imageOne][2];
            
            double alphaTwo = wobbles[imageTwo][0];
            double betaTwo = wobbles[imageTwo][1];
            double gammaTwo = wobbles[imageTwo][2];
            
            mat->rotate(-alphaOne, -betaOne, -gammaOne);
            mat->rotate(alphaTwo, betaTwo, gammaTwo);
        }
    }
    
    return mat;
}

double ImageCluster::scoreWrapper(void *object)
{
    double score = static_cast<ImageCluster *>(object)->selfConsistencyScore(true);
    
    return score;
}

int ImageCluster::uniqueUnexpectedMatches()
{
    int total = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        total += images[i]->hasAnyUnexpectedMatch();
    }
    
    return total;
}

bool ImageCluster::nextPermutation(std::vector<Image *> **testImagePointer, Image *image, int count)
{
    if (goodPermutations.size() == 0)
    {
        return std::next_permutation((*testImagePointer)->begin(), (*testImagePointer)->end());
    }
    else if (count > goodPermutations.size())
    {
        return false;
    }
    else
    {
        *testImagePointer = &goodPermutations[image][count];
        return true;
    }
}

double ImageCluster::distanceBetweenSpotPair(std::pair<SpotPtr, SpotPtr> pair, bool addWobble)
{
    Image *one = pair.first->getParentImage();
    Image *two = pair.second->getParentImage();
    
    MatrixPtr matOne = getImageToImageMatrix(images[0], one, addWobble);
    MatrixPtr matTwo = getImageToImageMatrix(images[0], two, addWobble);
    
    vec vecOne = pair.first->estimatedVector();
    vec vecTwo = pair.second->estimatedVector();
    
    matOne->multiplyVector(&vecOne);
    matTwo->multiplyVector(&vecTwo);
    
    vec between = vector_between_vectors(vecOne, vecTwo);
    
    return length_of_vector(between);
}

double ImageCluster::selfConsistencyScore(bool addWobble, int maxSuccesses, int maxTrials)
{
    if (commonSpots.size() == 0)
    {
        findSpotPairs();
    }
    
    double totalDistance = 0;
    
    for (int i = 0; i < commonSpots.size(); i++)
    {
        totalDistance += distanceBetweenSpotPair(commonSpots[i], addWobble);
    }
    
    totalDistance /= commonSpots.size();
    
    lastScore = totalDistance;
    
    return totalDistance;
}
/*
double ImageCluster::selfConsistencyScore(bool addWobble, int maxSuccesses, int maxTrials)
{
    double totalShift = 0;
    int shiftCount = 0;
    int successfulMatches = 0;
    
    bool findingPermutations = (goodPermutations.size() == 0);
    
    for (int i = 0; i < images.size(); i++)
    {
        if (!images[i]->isUnexpectedMatch(this))
            continue;
        
        MatrixPtr compounded = MatrixPtr(new Matrix());
        vec absolute = new_vector(1, 0, 0);
        
        logged << "Checking image " << i << std::endl;
        sendLog(LogLevelDetailed);
        
        std::vector<Image *>testImages = images;
        testImages.erase(testImages.begin() + i);
        
        std::vector<Image *> *testImagePointer = &testImages;
        
        int successfulTrials = 0;
        int trials = 0;
        
        while (nextPermutation(&testImagePointer, images[i], trials) && successfulTrials < maxSuccesses && trials < maxTrials)
        {
            trials++;
            
            vec test = new_vector(1, 0, 0);
        
            MatrixPtr nextMat = getImageToImageMatrix(images[i], testImages[0], addWobble);
            nextMat->multiplyVector(&test);
            
            int prevImage = 0;
            
            for (int j = 1; j < testImages.size(); j++)
            {
                MatrixPtr nextMat = getImageToImageMatrix(testImages[prevImage], testImages[j], addWobble);
                nextMat->multiplyVector(&test);
            
                prevImage = j;
            }
            
            nextMat = getImageToImageMatrix(testImages[prevImage], images[i], addWobble);
            nextMat->multiplyVector(&test);

            vec shiftVec = vector_between_vectors(absolute, test);
            double shift = length_of_vector(shiftVec);
        
            if (shift > 0.00001)
            {
                totalShift += shift;
                shiftCount++;
                successfulTrials++;
                
                if (findingPermutations)
                {
                    goodPermutations[images[i]].push_back(testImages);
                }
            }
        }
        
        if (successfulTrials > 0)
            successfulMatches++;
    }
    
    totalShift /= shiftCount;
    
    if (shiftCount == 0)
    {
        logged << "Failed to calculate a self-consistency score." << std::endl;
        sendLog(LogLevelDetailed);
        return 10;
    }
    
    lastScore = totalShift;
    
    sendLog(LogLevelDebug);
    
    logged << "Cluster with " << images.size() << " images has " << uniqueUnexpectedMatches() << " unexpected matches with average shift " << totalShift << " (" << successfulMatches << ")" << std::endl;
    
    if (uniqueUnexpectedMatches() > 5 && !addWobble)
        sendLog();

    sendLog(LogLevelDetailed);
    
    return totalShift;
}*/

void ImageCluster::breakUp(ImageClusterPtr toBreak, bool recursive)
{
    toBreak->breakUp(recursive);
    
    for (int i = 0; i < clusters.size(); i++)
    {
        if (clusters[i] == toBreak)
        {
            clusters.erase(clusters.begin() + i);
        }
    }
}

void ImageCluster::breakUp(bool recursive)
{
    if (recursive)
    {
        for (int i = 0; i < clusters.size(); i++)
        {
            breakUp(clusters[i], recursive);
        }
    }
    
    if (clusters.size() >= 2)
    {
        ImageClusterPtr subA = clusters[0];
        ImageClusterPtr subB = clusters[1];
        subA->flaggedForMerging = false;
        subB->flaggedForMerging = false;
        
        subA->setParentCluster(parentCluster);
        subB->setParentCluster(parentCluster);
        
        parentCluster->clusters.push_back(subA);
        parentCluster->clusters.push_back(subB);
        
        for (int i = 0; i < this->images.size(); i++)
        {
            images[i]->removeUnexpectedMatch(this);
        }
    }
    
    if (clusters.size() == 0)
    {
        for (int i = 0; i < images.size(); i++)
        {
            images[i]->setSeeded(false);
        }
    }
    
    if (clusters.size() >= 2)
        logged << "Breaking up " << images.size() << " image cluster into " << clusters[0]->images.size() << " and " << clusters[1]->images.size() << " image clusters." << std::endl;
    
    if (clusters.size() == 0)
    {
        logged << "Breaking up " << images.size() << " image seed cluster." << std::endl;
    }
    
    sendLog(LogLevelDetailed);
}

void ImageCluster::addTrust()
{
    trust += uniqueUnexpectedMatches();
    
    for (int i = 0; i < clusters.size(); i++)
    {
        clusters[i]->addTrust();
    }
}

void ImageCluster::merge(ImageClusterPtr subA, ImageClusterPtr subB, std::pair<Image *, Image *> sharedImages)
{
    int findCount = 0;
    subA->flaggedForMerging = false;
    subB->flaggedForMerging = false;

    for (int i = 0; i < clusters.size(); i++)
    {
        if (clusters[i] == subA || clusters[i] == subB)
            findCount++;
    }
    
    if (findCount < 2)
    {
        return;
    }
    
    std::vector<Image *> allImages;
    allImages.reserve(subA->images.size() + subB->images.size());
    allImages.insert(allImages.begin(), subA->images.begin(), subA->images.end());
    allImages.insert(allImages.begin() + subA->images.size(), subB->images.begin(), subB->images.end());
    std::vector<Image *> allUnexpectedMatches;
    
    for (int i = 1; i < allImages.size(); i++)
    {
        bool found = false;
        
        for (int j = 0; j < i && !found; j++)
        {
            if (allImages[i] == allImages[j])
            {
                allImages.erase(allImages.begin() + i);
                i--;
                found = true;
            }
        }
    }
    
    unsigned long moreUnexpected = subA->images.size() + subB->images.size() - allImages.size();
    
    ImageClusterPtr subAB = ImageClusterPtr(new ImageCluster(allImages, this));

    for (int i = 0; i < allUnexpectedMatches.size(); i++)
    {
        allUnexpectedMatches[i]->addUnexpectedMatch(&*subAB);
    }
    
    subA->setParentCluster(&*subAB);
    subB->setParentCluster(&*subAB);
    subAB->clusters.push_back(subA);
    subAB->clusters.push_back(subB);
    subAB->setSharedSpots(sharedImages.first, sharedImages.second);
    subAB->calculateMatrices(sharedImages.first, sharedImages.second);
    subAB->fillMatrixGaps();
    subAB->unexpectedMatches = subA->unexpectedMatches + subB->unexpectedMatches + (int)moreUnexpected;
    clusters.push_back(subAB);
    
    unsigned long clusterSize = clusters.size();
    
    for (int i = 0; i < clusters.size(); i++)
    {
        if (clusters[i] == subA)
        {
            clusters.erase(clusters.begin() + i);
            i--;
        }
        
        if (clusters[i] == subB)
        {
            clusters.erase(clusters.begin() + i);
            i--;
        }
    }
    
    double score = subAB->selfConsistencyScore();

    logged << "Clusters with " << subA->images.size() << " and " << subB->images.size() << " images merge to give " << allImages.size() << " images in supercluster (score = " << score << ")." << std::endl;

    double clusterMergeThreshold = FileParser::getKey("CLUSTER_MERGE_THRESHOLD", 0.15);
    
    if (score > clusterMergeThreshold)
    {
        subAB->minimizeWobble();
        
        if (subAB->lastScore < clusterMergeThreshold)
        {
            logged << "Further reduced to " << subAB->lastScore << " and thus accepted." << std::endl;
            sendLog(LogLevelNormal);
            return;
        }
        
        if (subAB->images.size() > 40)
            sendLog();
        else
            sendLog(LogLevelDetailed);
        breakUp(subAB);
        
        return;
    }
    else
    {
        sendLog(LogLevelNormal);
        addTrust();
        
        return;
    }
}

Image *ImageCluster::sharesCommonImagesWithCluster(ImageClusterPtr otherCluster)
{
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < otherCluster->images.size(); j++)
        {
            if (images[i] == otherCluster->images[j])
                return images[i];
        }
    }
    
    return NULL;
}

void ImageCluster::setParentCluster(ImageCluster *newParent)
{
    parentCluster = newParent;
    copyMatricesToParent();
}

void ImageCluster::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

void ImageCluster::writeRotations(std::string filename)
{
    std::ofstream rotFile;
    rotFile.open(filename);
    int num = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        vec test = new_vector(1, 0, 0);
        MatrixPtr relativeMat = imageToImageMatrices[images[i]][images[0]];

        relativeMat->multiplyVector(&test);
        
        rotFile << images[i]->getFilename() << "\t" << test.h << "\t" << test.k << "\t" << test.l << std::endl;
    }
    
    rotFile.close();
}

void ImageCluster::writeSpots(std::string filename)
{
    std::ofstream spotFile;
    spotFile.open(filename);
    int num = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        MatrixPtr relativeMat = imageToImageMatrices[images[0]][images[i]];

        std::vector<std::string> elements;
        std::vector<vec> positions;
        images[i]->rotatedSpotPositions(relativeMat, &positions, &elements);
        
        for (int j = 0; j < positions.size(); j++)
        {
            num++;
            
            char chainNum = i + 65;
            
            spotFile << "ATOM      4  OW  WAT " << chainNum;
            spotFile << std::setw(4) << j << "    ";
            spotFile << std::setprecision(3) << std::setw(8) << positions[j].h * 50;
            spotFile << std::setprecision(3) << std::setw(8) << positions[j].k * 50;
            spotFile << std::setprecision(3) << std::setw(8) << positions[j].l * 50;
            spotFile << "                       " << elements[j] << "  " << std::endl;
        }
    }
    
    logged << filename << "\t" << images.size() << /*"\t" << uniqueUnexpectedMatches() << */"\t" << lastScore << std::endl;
    
    sendLog();
    
    spotFile.close();
}
