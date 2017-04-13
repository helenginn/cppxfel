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

#define ANGLE_FUNNEL_START 1.5
#define ANGLE_DISTANCE_BUFFER 0.02

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
    distCSV->setupHistogram(0, maxDistance, powderStep, "Distance", 2, "Perfect frequency", "Correction");

    double km = FileParser::getKey("INITIAL_BANDWIDTH", 0.0013);
    double n = 2;
    km = pow(km, n);
    
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

	CSVPtr angleCSV = CSVPtr(new CSV(3, "angle", "aLength", "bLength"));

    for (int i = 0; i < standardVectorCount(); i++)
    {
        double distance = standardVector(i)->distance();

		if (distance <= 0)
			continue;

        if (distance > maxAngleDistance)
            continue;

        for (int j = 0; j < standardVectorCount(); j++)
        {
            double jDistance = standardVector(j)->distance();

			if (jDistance > distance)
				continue;

            if (jDistance > maxAngleDistance)
                continue;
            
            double angle = standardVector(j)->angleWithVector(standardVector(i));

            angle *= 180 / M_PI;
            angle = (angle > 90) ? 180 - angle : angle;

			if (angle != angle)
			{
				continue;
			}

			angleCSV->addEntry(3, angle, distance, jDistance);
        }
    }

	angleCSV->plotPDB("perfect_angles.pdb", "aLength", "bLength", "angle");

    angleCSV->writeToFile("perfect_angles.csv");
    distCSV->writeToFile("perfect_unit_cell.csv");
    weightedUnitCell = distCSV;
    weightedAngles = angleCSV;
    distanceToAngleRatio = distTotal / angleTotal;

	double maxAngle = 90.;
	double intervals = LOOKUP_INTERVALS;
	double lengthStep = maxAngleDistance / intervals;
	double angleStep = maxAngle / intervals;

	double angleTolerance = ANGLE_FUNNEL_START * 1.1;
	double lengthTolerance = 1.1 * ANGLE_FUNNEL_START * maxAngleDistance / 90;
	double total = pow(intervals, 3);
	memset(lookupIntervals, 0, total * sizeof(float));

	time_t startTime;
	time(&startTime);
	logged << "Calculating vector pair scores in advance..." << std::endl;
	sendLog();

	for (int dist1 = 0; dist1 < intervals; dist1++)
	{
		logged << ".";
		sendLog();

		std::vector<std::vector<double> > dist1Filter;
		double chosenDist1 = dist1 * lengthStep;

		for (int l = 0; l < angleCSV->entryCount(); l++)
		{
			double paLength = angleCSV->valueForEntry(1, l);
			double aLengthDiff = fabs(paLength - chosenDist1);
			if (aLengthDiff > lengthTolerance)
			{
				continue;
			}

			dist1Filter.push_back(angleCSV->entry(l));
		}

		for (int dist2 = 0; dist2 < intervals; dist2++)
		{
			double chosenDist2 = dist2 * lengthStep;
			std::vector<std::vector<double> > dist2Filter;

			for (int l = 0; l < dist1Filter.size(); l++)
			{
				double pbLength = dist1Filter[l][2];

				double bLengthDiff = fabs(pbLength - chosenDist2);
				if (bLengthDiff > lengthTolerance)
				{
					continue;
				}

				dist2Filter.push_back(dist1Filter[l]);
			}

			for (int angles = 0; angles < intervals; angles++)
			{
				double bestDist = FLT_MAX;

				double chosenAngle = ((double)angles) * angleStep;

				for (int l = 0; l < dist2Filter.size(); l++)
				{
					double pAngle = dist2Filter[l][0];
					double paLength = dist2Filter[l][1];
					double pbLength = dist2Filter[l][2];
					double aLengthDiff = fabs(paLength - chosenDist1);
					double bLengthDiff = fabs(pbLength - chosenDist2);

					double angleDiff = fabs(pAngle - chosenAngle);
					if (angleDiff > angleTolerance)
					{
						continue;
					}

					double scaleUp = 90. / maxAngleDistance;

					aLengthDiff *= scaleUp;
					bLengthDiff *= scaleUp;

					double distSq = angleDiff * angleDiff +
					aLengthDiff * aLengthDiff +
					bLengthDiff * bLengthDiff;

					if (distSq < bestDist)
					{
						bestDist = distSq;
					}
				}

				int position = dist1 * intervals * intervals +
				dist2 * intervals + angles;
				bestDist = sqrt(bestDist);

				if (bestDist > ANGLE_FUNNEL_START)
				{
					lookupIntervals[position] = -1;
				}
				else
				{
					bestDist *= M_PI / 180.;
					double score = 2 * sin((60 / ANGLE_FUNNEL_START) * bestDist + M_PI / 6) - 1;
					lookupIntervals[position] = -score;
				}


			}
		}
	}

	time_t endTime;
	time(&endTime);

	double seconds = endTime - startTime;
	logged << "Pre-calculation took " << seconds << " seconds." << std::endl;
	sendLog();
}

void UnitCellLattice::updateUnitCellData()
{
    lockUnitCellDimensions();
	std::vector<double> unitCell = getUnitCell();
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

void UnitCellLattice::setup()
{
    setupLock.lock();
    
    if (setupLattice)
    {
        setupLock.unlock();
        return;
    }
    
    powderStep = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
	maxAngleDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);

	if (maxAngleDistance > 0)
	{
		maxAngleDistance += ANGLE_DISTANCE_BUFFER;
	}

	std::vector<double> unitCell = getUnitCell();

	if (!symOperators.size())
    {
        Matrix::symmetryOperatorsForSpaceGroup(&symOperators, getSpaceGroup(), unitCell);
    }

	unitCellOnly = Matrix::matrixFromUnitCell(unitCell);
    
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    maxDistance = 0;
    
    int maxMillerIndexTrialH, maxMillerIndexTrialK, maxMillerIndexTrialL;
	double resolution = 1 / (FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15) + DISTANCE_BUFFER);

    getMaxMillerIndicesForResolution(resolution, &maxMillerIndexTrialH, &maxMillerIndexTrialK, &maxMillerIndexTrialL);
    
    int count = 0;
    
    orderedDistances.push_back(0);
    
    for (int i = -maxMillerIndexTrialH; i <= maxMillerIndexTrialH; i++)
    {
        for (int j = -maxMillerIndexTrialK; j <= maxMillerIndexTrialK; j++)
        {
            for (int k = -maxMillerIndexTrialL; k <= maxMillerIndexTrialL; k++)
            {
                if (ccp4spg_is_sysabs(getSpaceGroup(), i, j, k))
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
                
                int asym = CSym::ccp4spg_is_in_asu(getSpaceGroup(), i, j, k);
                
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

double UnitCellLattice::weightForDistance(double distance)
{
    double correctedDistance = distance + powderStep / 2;
    return correctedDistance * weightedUnitCell->valueForHistogramEntry(1, correctedDistance);
}

double UnitCellLattice::weightForPair(double dist1, double dist2, double angle)
{
	if (dist1 > maxAngleDistance || dist2 > maxAngleDistance)
	{
		return -1;
	}

	if (dist1 != dist1 || dist2 != dist2 || angle != angle)
	{
		return -1;
	}

	int dist1Step = dist1 * LOOKUP_INTERVALS / (maxAngleDistance);
	int dist2Step = dist2 * LOOKUP_INTERVALS / (maxAngleDistance);
	int angleStep = angle * LOOKUP_INTERVALS / 90.;

	int position = dist1Step * LOOKUP_INTERVALS * LOOKUP_INTERVALS +
	dist2Step * LOOKUP_INTERVALS + angleStep;

	return lookupIntervals[position];
}

UnitCellLatticePtr UnitCellLattice::getMainLattice()
{
    if (!mainLattice)
    {
        mainLattice = UnitCellLatticePtr(new UnitCellLattice());
    }
    
    return mainLattice;
}