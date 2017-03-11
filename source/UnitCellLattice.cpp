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

bool UnitCellLattice::setupLattice = false;
std::vector<MatrixPtr> UnitCellLattice::symOperators;
UnitCellLatticePtr UnitCellLattice::mainLattice;

void UnitCellLattice::getMaxMillerIndicesForResolution(double resolution, int *hMax, int *kMax, int *lMax)
{
    double lengths[3];
    unitCellMatrix->unitCellLengths(lengths);
    
    *hMax = lengths[0] / resolution;
    *kMax = lengths[1] / resolution;
    *lMax = lengths[2] / resolution;
}

void UnitCellLattice::weightUnitCell()
{
    CSVPtr distCSV = CSVPtr(new CSV(0));
    double maxDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15) + DISTANCE_BUFFER;
    double maxAngleDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.04);
    distCSV->setupHistogram(0, maxDistance, powderStep, "Distance", 2, "Perfect frequency", "Correction");
    CSVPtr angleCSV = CSVPtr(new CSV(0));
    angleCSV->setupHistogram(0, 90, 0.1, "Angle", 2, "Perfect frequency", "Cosine");
    
    double km = FileParser::getKey("INITIAL_BANDWIDTH", 0.0013);
    double n = 2;
    km = pow(km, n);
    
    for (int i = 0; i < angleCSV->entryCount(); i++)
    {
        double angle = angleCSV->valueForEntry("Angle", i);
        double cosine = cos(angle * M_PI / 180);
        angleCSV->setValueForEntry(i, "Cosine", cosine);
    }
    
    double distTotal = 0;
    double angleTotal = 0;
    int csvPos = 0;
    
    double normDist = 0;
    
    for (int i = 0; i < orderedDistanceCount() - 1; i++)
    {
        double firstDistance = orderedDistance(i);
        double secondDistance = orderedDistance(i + 1);
        double separation = secondDistance - firstDistance;
        double vmax = 100000 * separation;
        
        if (separation < 0.000001)
        {
            continue;
        }
        
        if (normDist == 0)
        {
            normDist = sqrt(secondDistance);
        }
        
        while (true)
        {
            if (csvPos >= distCSV->entryCount())
                break;
            
            double csvDist = distCSV->valueForEntry("Distance", csvPos);
            
            if (csvDist > secondDistance)
                break;
            
            double x = std::min(csvDist - firstDistance, secondDistance - csvDist);
            
            double y = vmax * pow(x, n) / (km + pow(x, n));
            
            distCSV->setValueForEntry(csvPos, "Perfect frequency", y);
            
            csvPos++;
        }
    }
    
    double averageableDist = maxDistance / 8;
    int stepsPerAverage = (averageableDist / powderStep + 0.5) / 2;
    
    csvPos = stepsPerAverage;
    
    while (true)
    {
        if (csvPos >= distCSV->entryCount() - stepsPerAverage)
            break;
        
        double integration = 0;
        
        for (int i = -stepsPerAverage; i <= stepsPerAverage; i++)
        {
            int csvPosAdd = csvPos + i;
            
            double value = distCSV->valueForEntry("Perfect frequency", csvPosAdd);
            
            integration += value * powderStep;
        }
        
        integration = 1 / integration;
        
        distCSV->setValueForEntry(csvPos, "Correction", integration);

        csvPos++;
    }
    
    csvPos = 0;
    
    while (true)
    {
        int checkPos = csvPos;
        
        if (checkPos < stepsPerAverage)
        {
            checkPos = stepsPerAverage;
        }
        
        if (checkPos > distCSV->entryCount() - stepsPerAverage)
        {
            checkPos = distCSV->entryCount() - stepsPerAverage;
        }
        
        if (csvPos >= distCSV->entryCount())
            break;
        
        double correction = distCSV->valueForEntry("Correction", checkPos);
        double weight = distCSV->valueForEntry("Perfect frequency", csvPos);
        
        weight *= correction;
        distCSV->setValueForEntry(csvPos, "Perfect frequency", weight);
        
        csvPos++;
    }
    
    for (int i = 0; i < uniqueSymVectorCount(); i++)
    {
        double distance = uniqueSymVector(i)->distance();
        
        if (distance <= 0)
            continue;
        
        double inverse = 1 / distance;
        inverse *= maxDistance;
        
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
    setupLock.lock();
    
    if (setupLattice)
    {
        setupLock.unlock();
        return;
    }
    
    powderStep = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
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
    
    maxDistance = 0;
    
    int maxMillerIndexTrialH, maxMillerIndexTrialK, maxMillerIndexTrialL;
    
    if (resolution == 0)
    {
        resolution = 1 / (FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15) + DISTANCE_BUFFER);
    }
    
    getMaxMillerIndicesForResolution(resolution, &maxMillerIndexTrialH, &maxMillerIndexTrialK, &maxMillerIndexTrialL);
    
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
                
                int asym = CSym::ccp4spg_is_in_asu(spaceGroup, i, j, k);
                
                spotVectors.push_back(newStandardVector);
                orderedDistances.push_back(newStandardVector->distance());
                
                if (asym)
                {
                    uniqueSymVectors.push_back(newStandardVector);
                }
                
                count++;
            }
        }
    }
    
    std::sort(orderedDistances.begin(), orderedDistances.end(), std::less<double>());
    
    weightUnitCell();
    
    minDistance = FLT_MAX;
    
    for (int i = 0; i < 3; i++)
    {
        vec hkl = new_vector((i == 0), (i == 1), (i == 2));
        unitCellMatrix->multiplyVector(&hkl);
        
        if (length_of_vector(hkl) < minDistance)
            minDistance = length_of_vector(hkl);
    }
    
    setupLock.unlock();
    setupLattice = true;
}

void UnitCellLattice::refineUnitCell(PowderHistogram aHistogram)
{
    histogram = aHistogram;
}

double UnitCellLattice::weightForDistance(double distance)
{
    double correctedDistance = distance + powderStep / 2;
    return correctedDistance * weightedUnitCell->valueForHistogramEntry(1, correctedDistance);
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