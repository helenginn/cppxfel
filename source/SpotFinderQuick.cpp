//
//  SpotFinderQuick.cpp
//  cppxfel
//
//  Created by Helen Ginn on 29/04/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotFinderQuick.h"
#include "Image.h"

template<class Value>
void SpotFinderQuick::findSignalToNoise(Value *data, size_t position, int xDim, int yDim, float *signalToNoiseRatio, float *background, float *backgroundVariance)
{
    int backgroundSum = 0;
    int backgroundSumSquared = 0;
    
    for (int k = 0; k < backgroundShifts.size(); k++)
    {
        int shift = backgroundShifts[k];
        size_t absolutePos = position + shift;
        
        if (absolutePos >= xDim * yDim)
        {
            continue;
        }
        
        backgroundSum += data[absolutePos];
        backgroundSumSquared += pow(data[absolutePos], 2);
    }
    
    backgroundSum /= backgroundShifts.size();
    *background = backgroundSum;
    backgroundSumSquared /= backgroundShifts.size();
    
    *backgroundVariance = backgroundSumSquared - pow(backgroundSum, 2);
    *signalToNoiseRatio = data[position] / *backgroundVariance;
}

void SpotFinderQuick::findSpecificSpots(std::vector<SpotPtr> *spots)
{    
    int xDim = image->getXDim();
    int yDim = image->getYDim();
    float minSeparationSquared = minSeparation * minSeparation;
    int peakToEnter = 0;
    
    bool useShort = true;
    int *data = NULL;
    short int *shortData = image->getShortDataPtr();
    
    size_t *pixelTracker = (size_t *)malloc(maxPixels * sizeof(size_t));
    peaks = (Peak *)malloc(sizeof(Peak) * maxHits);

    // calloc sets mem to 0... 1 will be used as masked
    int *trackerMask = (int *)calloc(xDim * yDim, sizeof(int));
    
    int reachedSNRThreshold = 0;
    int reachedThreshold = 0;
    int tooBig = 0;
    int tooSmall = 0;
    
    if (shortData == NULL)
    {
        useShort = false;
        data = image->getDataPtr();
    }
    
    int shifts[] = { - xDim - 1, - xDim, - xDim + 1,
        -1, 1,
        + xDim + 1, + xDim, +xDim + 1 };
    
    for (int i = 0; i < yDim; i++)
    {
        for (int j = 0; j < xDim; j++)
        {
            if (peakToEnter + 1 > maxHits)
            {
                break;
            }

            bool mustContinue = false;
            size_t position = i * xDim + j;
            
            short int value = (useShort ? shortData[position] : data[position]);
            
            if (value < threshold)
            {
                continue;
            }
            
            reachedThreshold++;
            
            for (int k = 0; k < 8; k++)
            {
                size_t otherPosition = position + shifts[k];
                
                if (otherPosition >= xDim * yDim)
                    continue;
                
                short int otherValue = (useShort ? shortData[otherPosition] : data[otherPosition]);
                
                if (value <= otherValue)
                {
                    mustContinue = true;
                    break;
                }
            }
            
            if (mustContinue)
            {
                continue;
            }
            
            float signalToNoiseRatio = 0;
            float background = 0;
            float backgroundVariance = 0;
            
            if (shortData)
            {
                findSignalToNoise(shortData, position, xDim, yDim, &signalToNoiseRatio, &background, &backgroundVariance);
            }
            else
            {
                findSignalToNoise(data, position, xDim, yDim, &signalToNoiseRatio, &background, &backgroundVariance);
            }
            
      //      std::cout << signalToNoiseRatio << std::endl;
            
            if (signalToNoiseRatio < signalToNoiseThreshold)
            {
                continue;
            }
            
            reachedSNRThreshold++;
            
            float backgroundSigma = sqrt(backgroundVariance);
            
            float centreOfMassX = i;
            float centreOfMassY = j;
            pixelTracker[0] = position;
            int totalPixelsToCheck = 1;
            int lastCheckedIndex = 0;
            double totalSignal = 0;
            
            do
            {
                bool mustBreak = false;
                size_t currentPixelToCheck = pixelTracker[lastCheckedIndex];
                
                for (int k = 0; k < 8; k++)
                {
                    size_t relativeToCurrentPixel = currentPixelToCheck + shifts[k];
                    
                    // don't need to check for 0 because it is unsigned!
                    if (relativeToCurrentPixel >= xDim * yDim)
                    {
                        continue;
                    }
                    
                    if (trackerMask[relativeToCurrentPixel] == 1)
                    {
                        continue;
                    }
                    
                    int currentPixelValue = (useShort ? shortData[relativeToCurrentPixel] : data[relativeToCurrentPixel]);
                    float currentSignalToNoise = (currentPixelValue - background) / backgroundSigma;
                    
                    if (currentSignalToNoise > signalToNoiseThreshold)
                    {
                        trackerMask[relativeToCurrentPixel] = 1;
                        
                        if (totalPixelsToCheck == maxPixels)
                        {
                            mustBreak = true;
                            break;
                        }

                        pixelTracker[totalPixelsToCheck] = relativeToCurrentPixel;
                        totalPixelsToCheck++;
                        
                        float signal = currentPixelValue - background;
                        totalSignal += signal;
                        size_t currentY = relativeToCurrentPixel / xDim;
                        
                        if (xDim == 0)
                        {
                            logged << "Warning! xDim is zero - something is very wrong" << std::endl;
                            sendLogAndExit();
                        }
                        
                        size_t currentX = relativeToCurrentPixel % xDim;
                        
                        centreOfMassX += currentX * signal;
                        centreOfMassY += currentY * signal;
                    }
                    
                    if (mustBreak)
                    {
                        break;
                    }
                }
                
                if (mustBreak)
                {
                    break;
                }
                
                lastCheckedIndex++;
            }
            while (lastCheckedIndex != totalPixelsToCheck);
            
        //    std::cout << totalPixelsToCheck << std::endl;
            
            if (totalPixelsToCheck >= maxPixels)
            {
                tooBig++;
                continue;
            }
                
                
            if (totalPixelsToCheck < minPixels)
            {
                tooSmall++;
                continue;
            }
            
            centreOfMassX /= totalSignal;
            centreOfMassY /= totalSignal;
            
            for (int k = 0; k < peakToEnter; k++)
            {
                float distanceSquared = pow(centreOfMassX - peaks[k].centreX, 2) + pow(centreOfMassY - peaks[k].centreY, 2);
                
                if (distanceSquared < minSeparationSquared)
                {
                    mustContinue = true;
                }
            }
            
            if (mustContinue)
            {
                continue;
            }
            
            peaks[peakToEnter].centreX = centreOfMassX;
            peaks[peakToEnter].centreY = centreOfMassY;
            peaks[peakToEnter].totalSignal = totalSignal;
            
            peakToEnter++;
        }
    }
    
    std::ostringstream logged;
    logged << "Pixels reaching intensity threshold: " << reachedThreshold << std::endl;
    logged << "Pixels reaching SNR threshold: " << reachedSNRThreshold << std::endl;
    logged << "Spot had too many pixels: " << tooBig << std::endl;
    logged << "Spot had too few pixels: " << tooSmall << std::endl;
    
    Logger::log(logged);
    
    totalPeaks = peakToEnter + 1;
    
    free(pixelTracker);
    free(trackerMask);
}

void SpotFinderQuick::calculateBackgroundShifts()
{
    int xDim = image->getXDim();
    
    for (int i = -maxRadius; i <= maxRadius; i++)
    {
        for (int j = -maxRadius; j <= maxRadius; j++)
        {
            if (i > -minRadius && i <= minRadius && j > -minRadius && j <= minRadius)
            {
                continue;
            }
            
            int relativePos = i * xDim + j;
            
            backgroundShifts.push_back(relativePos);
        }
    }
}