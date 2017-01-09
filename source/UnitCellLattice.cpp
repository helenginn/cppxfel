//
//  UnitCellLattice.cpp
//  cppxfel
//
//  Created by Helen Ginn on 10/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "UnitCellLattice.h"
#include "csymlib.h"
#include "FileParser.h"
#include "SpotVector.h"
#include "Logger.h"
#include "Detector.h"
#include "CSV.h"
#include "Vector.h"
#include <algorithm>

std::vector<MatrixPtr> UnitCellLattice::symOperators;

void UnitCellLattice::getMaxMillerIndicesForResolution(double resolution, int *hMax, int *kMax, int *lMax)
{
    double *lengths = new double[3];
    
    unitCellMatrix->unitCellLengths(&lengths);
    
    *hMax = lengths[0] / resolution;
    *kMax = lengths[1] / resolution;
    *lMax = lengths[2] / resolution;
    
    delete [] lengths;
}

void UnitCellLattice::addConvolutedPeak(CSVPtr csv, double mean, double stdev, double weight)
{
    double totalIntervals = 300;
    double step = (stdev * 6) / totalIntervals;
    
    for (double x = -stdev * 3; x < stdev * 3; x += step)
    {
        double y = normal_distribution(x, 0, stdev) / totalIntervals;
        csv->addOneToFrequency(x + mean, "Perfect frequency", y * weight);
    }
}

void UnitCellLattice::weightUnitCell()
{
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
    double rlpSize = FileParser::getKey("INITIAL_RLP_SIZE", 0.0001);
    CSVPtr csv = CSVPtr(new CSV(0));
    csv->setupHistogram(0, maxDistance, step, "Distance", 1, "Perfect frequency");
    
    for (int i = 0; i < uniqueSymVectorCount(); i++)
    {
        double distance = uniqueSymVector(i)->distance();
        
        if (distance <= 0)
            continue;
        
        double inverse = 1 / distance;
        inverse *= maxDistance;
        double stdev = rlpSize / 2;
        
        addConvolutedPeak(csv, distance, stdev, inverse);
    }
    
    csv->writeToFile("test_perfect.csv");
    weightedUnitCell = csv;
}

void UnitCellLattice::setup(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum, double resolution)
{
    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);
    
    if (!symOperators.size())
    {
        Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup, a, b, c, alpha, beta, gamma);
    }
    
    unitCellOnly = Matrix::matrixFromUnitCell(a, b, c, alpha, beta, gamma);
    
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 0);
    maxDistance = 0;
    
    int maxMillerIndexTrialH, maxMillerIndexTrialK, maxMillerIndexTrialL;
    
    if (resolution == 0 && maxMillerIndexTrial == 0)
    {
        resolution = 1 / FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    }
    
    if (resolution == 0)
    {
        maxMillerIndexTrialH = maxMillerIndexTrial;
        maxMillerIndexTrialK = maxMillerIndexTrial;
        maxMillerIndexTrialL = maxMillerIndexTrial;
    }
    else
    {
        getMaxMillerIndicesForResolution(resolution, &maxMillerIndexTrialH, &maxMillerIndexTrialK, &maxMillerIndexTrialL);
    }
    
    int count = 0;
    
    orderedDistances.push_back(0);
    
    for (int i = -maxMillerIndexTrialH; i <= maxMillerIndexTrialH; i++)
    {
        for (int j = -maxMillerIndexTrialK; j <= maxMillerIndexTrialK; j++)
        {
            for (int k = -maxMillerIndexTrialL; k <= maxMillerIndexTrialL; k++)
            {
                if (ccp4spg_is_sysabs(spaceGroup, i, j, k))
                {
                    continue;
                }
                
                vec hkl = new_vector(i, j, k);
                vec hkl_transformed = copy_vector(hkl);
                
                unitCellMatrix->multiplyVector(&hkl_transformed);
                
                double distance = length_of_vector(hkl_transformed);
                
                if (distance > 1 / resolution)
                    continue;
                
                if (distance > maxDistance)
                    maxDistance = distance;
                
                SpotVectorPtr newStandardVector = SpotVectorPtr(new SpotVector(hkl_transformed, hkl));
                
                vec3<int> integer = vec3<int>(i, j, k);
                integerVectors.push_back(integer);
                int asym = CSym::ccp4spg_is_in_asu(spaceGroup, i, j, k);
                
                spotVectors.push_back(newStandardVector);
                
                if (asym)
                {
                    uniqueSymVectors.push_back(newStandardVector);
                    orderedDistances.push_back(newStandardVector->distance());
                }
                
                count++;
            }
        }
    }
    
    weightUnitCell();
    
    std::sort(orderedDistances.begin(), orderedDistances.end(), std::less<double>());
    
    minDistance = FLT_MAX;
    
    for (int i = 0; i < 3; i++)
    {
        vec hkl = new_vector((i == 0), (i == 1), (i == 2));
        unitCellMatrix->multiplyVector(&hkl);
        
        if (length_of_vector(hkl) < minDistance)
            minDistance = length_of_vector(hkl);
    }
}

void UnitCellLattice::refineUnitCell(PowderHistogram aHistogram)
{
    histogram = aHistogram;
}

double UnitCellLattice::weightForDistance(double distance)
{
    return weightedUnitCell->valueForHistogramEntry("Perfect frequency", distance);
}