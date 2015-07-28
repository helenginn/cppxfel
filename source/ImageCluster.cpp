//
//  ImageCluster.cpp
//  cppxfel
//
//  Created by Helen Ginn on 28/07/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "ImageCluster.h"
#include "Image.h"
#include "Logger.h"

ImageCluster::ImageCluster(std::vector<Image *>theImages)
{
    images = theImages;
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images.size(); j++)
        {
            imageToImageMatrices[images[i]][images[j]] = MatrixPtr(new Matrix());
        }
    }
}

void ImageCluster::setSeed()
{
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->setSeeded();
    }
}

void ImageCluster::process()
{
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->commonLinesWithImages(images);
    }
    
    logged << "Determined common lines between images." << std::endl;
    sendLog();
    
    imageStatistics();
    
    findSeedingClusters();
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
            
            if (oneImage->hasCommonLinesWithImage(twoImage))
            {
                for (int k = j + 1; k < images.size(); k++)
                {
                    Image *threeImage = images[k];
                    
                    if (threeImage->hasCommonLinesWithImage(oneImage) &&
                        threeImage->hasCommonLinesWithImage(twoImage))
                    {
                        clusterCount++;
                        
                        std::vector<Image *>threeImages;
                        threeImages.push_back(oneImage);
                        threeImages.push_back(twoImage);
                        threeImages.push_back(threeImage);
                        
                        ImageClusterPtr newCluster = ImageClusterPtr(new ImageCluster(threeImages));
                        newCluster->setSeed();
                        clusters.push_back(newCluster);
                    }
                }
            }
        }
    }
    
    logged << "Found " << clusterCount << " three-image seeds." << std::endl;
    sendLog();
}

void ImageCluster::imageStatistics()
{
    unsigned long totalImages = images.size();
    int hopefulImages = 0;
    int commonLinesTotal = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        unsigned long pairCount = images[i]->linePairCount();
        
        if (pairCount > 0)
        {
            hopefulImages++;
            commonLinesTotal += pairCount;
        }
    }
    
    double commonLinesPerHopefulImage = (double)commonLinesTotal / (double)hopefulImages;
    
    logged << "Total images: " << totalImages << std::endl;
    logged << "Hopeful images (common line with at least one other image): " << hopefulImages << std::endl;
    logged << "Common lines per hopeful image: " << commonLinesPerHopefulImage << std::endl;

    sendLog();
}

void ImageCluster::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}
