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
UnitCellLatticePtr UnitCellLattice::mainLattice;

void UnitCellLattice::getMaxMillerIndicesForResolution(double resolution, int *hMax, int *kMax, int *lMax)
{
    double *lengths = new double[3];
    
    unitCellMatrix->unitCellLengths(&lengths);
    
    *hMax = lengths[0] / resolution;
    *kMax = lengths[1] / resolution;
    *lMax = lengths[2] / resolution;
    
    delete [] lengths;
}

void UnitCellLattice::weightUnitCell()
{
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
    double rlpSize = FileParser::getKey("INITIAL_RLP_SIZE", 0.0001);
    CSVPtr distCSV = CSVPtr(new CSV(0));
    double maxDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    double maxAngleDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.04);
    distCSV->setupHistogram(0, maxDistance, step, "Distance", 1, "Perfect frequency");
    CSVPtr angleCSV = CSVPtr(new CSV(0));
    angleCSV->setupHistogram(0, 90, 0.1, "Angle", 2, "Perfect frequency", "Cosine");
    
    for (int i = 0; i < angleCSV->entryCount(); i++)
    {
        double angle = angleCSV->valueForEntry("Angle", i);
        double cosine = cos(angle * M_PI / 180);
        angleCSV->setValueForEntry(i, "Cosine", cosine);
    }
    
    double distTotal = 0;
    double angleTotal = 0;
    
    for (int i = 0; i < uniqueSymVectorCount(); i++)
    {
        double distance = uniqueSymVector(i)->distance();
        
        if (distance <= 0)
            continue;
        
        double inverse = 1 / distance;
        inverse *= maxDistance;
        double stdev = rlpSize / 2;
        
        distCSV->addConvolutedPeak(1, distance, stdev, inverse);
        distTotal += inverse;
        
        if (distance > maxAngleDistance)
            continue;
        
        inverse /= uniqueSymVectorCount() * standardVectorCount();
        
        for (int j = 0; j < standardVectorCount(); j++)
        {
            if (j == i)
                continue;
            
            double jDistance = standardVector(j)->distance();
            
            if (jDistance > maxAngleDistance)
                continue;
            
            double angle = standardVector(j)->angleWithVector(uniqueSymVector(i));
            
            angle *= 180 / M_PI;
            angle = (angle > 90) ? 180 - angle : angle;
            
            angleCSV->addConvolutedPeak(1, angle, 0.3, inverse);
            angleTotal += inverse;
        }
    }
    
    angleCSV->writeToFile("perfect_angles.csv");
    distCSV->writeToFile("perfect_unit_cell.csv");
    weightedUnitCell = distCSV;
    weightedAngles = angleCSV;
    distanceToAngleRatio = distTotal / angleTotal;
}

void UnitCellLattice::updateUnitCellData()
{

    lockUnitCellDimensions(&_aDim, &_bDim, &_cDim, &_alpha, &_beta, &_gamma);
    unitCellOnly = Matrix::matrixFromUnitCell(unitCell);
    
    for (int i = 0; i < standardVectorCount(); i++)
    {
        SpotVectorPtr spotVec = standardVector(i);
        vec hkl = spotVec->getHKL();
        unitCellOnly->multiplyVector(&hkl);
        spotVec->setSpotDiff(hkl);
        spotVec->setUpdate();
    }
    
    orderedDistances.clear();

    for (int i = 0; i < uniqueSymVectorCount(); i++)
    {
        SpotVectorPtr uniqueVector = uniqueSymVector(i);
        orderedDistances.push_back(uniqueVector->distance());
    }
    
    std::sort(orderedDistances.begin(), orderedDistances.end(), std::less<double>());
    weightUnitCell();
}

void UnitCellLattice::setup(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum, double resolution)
{
    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);
    
    _aDim = a;
    _bDim = b;
    _cDim = c;
    _alpha = alpha;
    _beta = beta;
    _gamma = gamma;
    
    unitCell[0] = _aDim;
    unitCell[1] = _bDim;
    unitCell[2] = _cDim;
    unitCell[3] = _alpha;
    unitCell[4] = _beta;
    unitCell[5] = _gamma;
    
    if (!symOperators.size())
    {
        Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup, a, b, c, alpha, beta, gamma);
    }
    
    unitCellOnly = Matrix::matrixFromUnitCell(unitCell);
    
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 0);
    maxDistance = 0;
    
    int maxMillerIndexTrialH, maxMillerIndexTrialK, maxMillerIndexTrialL;
    
    if (resolution == 0 && maxMillerIndexTrial == 0)
    {
        resolution = 1 / (FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15) + 0.05);
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
    return weightedUnitCell->valueForHistogramEntry(1, distance);
}

double UnitCellLattice::weightForAngle(double angle)
{
    return weightedAngles->valueForHistogramEntry("Perfect frequency", angle) * distanceToAngleRatio;
}

double UnitCellLattice::weightForCosine(double cosine)
{
    return weightedAngles->valueForHistogramEntry("Perfect frequency", cosine, "Cosine") * distanceToAngleRatio;
}


void UnitCellLattice::lockUnitCellDimensions(double *a, double *b, double *c, double *alpha, double *beta, double *gamma)
{
    double spgNum = spaceGroup->spg_num;
    if (spgNum >= 75 && spgNum <= 194)
    {
        *b = *a;
    }
    if (spgNum >= 195)
    {
        *b = *a;
        *c = *a;
        *alpha = 90;
        *beta = 90;
        *gamma = 90;
    }
}

UnitCellLatticePtr UnitCellLattice::getMainLattice()
{
    if (!mainLattice)
    {
        double spaceGroupNum = FileParser::getKey("SPACE_GROUP", 0);
        
        std::vector<double> unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());
        
        mainLattice = UnitCellLatticePtr(new UnitCellLattice(unitCell[0], unitCell[1], unitCell[2],
                                                             unitCell[3], unitCell[4], unitCell[5], spaceGroupNum));
    }
    
    return mainLattice;
}