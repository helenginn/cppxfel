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
#include "Miller.h"
#include "RefinementStrategy.h"

#define ANGLE_FUNNEL_START 1.5
#define ANGLE_DISTANCE_BUFFER 0.005

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

void UnitCellLattice::weightUnitCellThread(void *object, int offset)
{
	UnitCellLattice *me = static_cast<UnitCellLattice *>(object);

	double indexingRlp = FileParser::getKey("INDEXING_RLP_SIZE", 0.001);
	int maxThreads = FileParser::getMaxThreads();

	double maxAngle = 90.;
	double intervals = LOOKUP_INTERVALS;

	int threadsPerInterval = intervals / (double)maxThreads + 1;
	int start = threadsPerInterval * offset;
	int end = threadsPerInterval * (offset + 1);

	double lengthStep = me->maxAngleDistance / intervals;
	double angleStep = maxAngle / intervals;

	double angleTolerance = ANGLE_FUNNEL_START * 1.0;
	double lengthTolerance = indexingRlp;

	for (int dist1 = start; dist1 < intervals && dist1 < end; dist1++)
	{
		std::ostringstream logged;
		logged << ".";
		Logger::log(logged);

		std::vector<std::vector<double> > dist1Filter;
		double chosenDist1 = ((double)dist1 + 0.5) * lengthStep;

		for (int l = 0; l < me->angleCSV->entryCount(); l++)
		{
			double paLength = me->angleCSV->valueForEntry(1, l);
			double aLengthDiff = fabs(paLength - chosenDist1);
			if (aLengthDiff > lengthTolerance)
			{
				continue;
			}

			dist1Filter.push_back(me->angleCSV->entry(l));
		}

		for (int dist2 = 0; dist2 < intervals; dist2++)
		{
			double chosenDist2 = ((double)dist2 + 0.5) * lengthStep;
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

				double chosenAngle = ((double)angles + 0.5) * angleStep;

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

					double scaleUp = ANGLE_FUNNEL_START / indexingRlp;

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
				double score = 0;
				double minAngleFunnelRad = atan(indexingRlp / chosenDist1)
				+ atan(indexingRlp / chosenDist2);

				double minAngleFunnelDeg = minAngleFunnelRad * 180 / M_PI;

				if (bestDist <= minAngleFunnelDeg)
				{
				//	bestDist /= minAngleFunnelDeg / 10;
				//	score = Detector::lookupCache(bestDist);
					score = 1 - bestDist / minAngleFunnelDeg;
				}

				me->lookupIntervals[position] = score;
			}
		}
	}
}

void UnitCellLattice::weightUnitCell()
{
	angleCSV = CSVPtr(new CSV(3, "angle", "aLength", "bLength"));

    for (int i = 0; i < uniqueSymVectorCount(); i++)
    {
        double distance = uniqueSymVector(i)->distance();

		if (distance <= 0)
			continue;

        if (distance > maxAngleDistance)
            continue;

        for (int j = 0; j < standardVectorCount(); j++)
        {
            double jDistance = standardVector(j)->distance();

		//	if (jDistance > distance)
		//		continue;

            if (jDistance > maxAngleDistance)
                continue;

			double minDistance = (distance > jDistance) ? jDistance : distance;
			double maxDistance = (distance <= jDistance) ? jDistance : distance;

            double angle = standardVector(j)->angleWithVector(standardVector(i));

            angle *= 180 / M_PI;
            angle = (angle > 90) ? 180 - angle : angle;

			if (angle != angle)
			{
				continue;
			}

			angleCSV->addEntry(3, angle, maxDistance, minDistance);
        }
    }

	angleCSV->plotPDB("perfect_angles.pdb", "aLength", "bLength", "angle");

    angleCSV->writeToFile("perfect_angles.csv");
    weightedAngles = angleCSV;

	double intervals = LOOKUP_INTERVALS;
	double total = pow(intervals, 3);

	memset(lookupIntervals, 0, total * sizeof(float));
	lookupIntervalPtr = &lookupIntervals[0];

	int maxThreads = FileParser::getMaxThreads();

	time_t startTime;
	time(&startTime);
	logged << "Calculating vector pair scores in advance..." << std::endl;
	sendLog();

	boost::thread_group threads;

	for (int i = 0; i < maxThreads; i++)
	{
		boost::thread *thr = new boost::thread(weightUnitCellThread, this, i);
		threads.add_thread(thr);
	}

	threads.join_all();

	time_t endTime;
	time(&endTime);

	double seconds = endTime - startTime;
	logged << "Pre-calculation took " << seconds << " seconds." << std::endl;
	sendLog();
}

void UnitCellLattice::updateUnitCellData()
{
    std::vector<double> newUnitCell;
	newUnitCell.push_back(_aDim);
	newUnitCell.push_back(_bDim);
	newUnitCell.push_back(_cDim);
	newUnitCell.push_back(_alpha);
	newUnitCell.push_back(_beta);
	newUnitCell.push_back(_gamma);
	setUnitCell(newUnitCell);
	lockUnitCellDimensions();
	newUnitCell = getUnitCell();

    MatrixPtr newUnitCellOnly = Matrix::matrixFromUnitCell(newUnitCell);

	unitCellOnly->copyComponents(newUnitCellOnly);

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
  //  weightUnitCell();
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

	if (unitCell.size() < 6)
	{
		logged << "Sorry, need unit cell dimensions specified by UNIT_CELL to do this." << std::endl;
		sendLogAndExit();
	}

	_aDim = unitCell[0];
	_bDim = unitCell[1];
	_cDim = unitCell[2];
	_alpha = unitCell[3];
	_beta = unitCell[4];
	_gamma = unitCell[5];

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

double UnitCellLattice::weightForPair(double dist1, double dist2, double angle)
{
	if (!weightedAngles)
	{
		std::lock_guard<std::mutex> lg(setupLock);

		if (!weightedAngles)
		{
			weightUnitCell();
		}
	}

	if (dist1 > maxAngleDistance || dist2 > maxAngleDistance)
	{
		return 0;
	}

	if (dist1 != dist1 || dist2 != dist2 || angle != angle)
	{
		return 0;
	}

	int dist1Step = dist1 * LOOKUP_INTERVALS / (maxAngleDistance);
	int dist2Step = dist2 * LOOKUP_INTERVALS / (maxAngleDistance);
	int angleStep = angle * LOOKUP_INTERVALS / 90.;

	int position = dist1Step * LOOKUP_INTERVALS * LOOKUP_INTERVALS +
	dist2Step * LOOKUP_INTERVALS + angleStep;

	return lookupIntervalPtr[position];
}

UnitCellLatticePtr UnitCellLattice::getMainLattice()
{
    if (!mainLattice)
    {
        mainLattice = UnitCellLatticePtr(new UnitCellLattice());
    }
    
    return mainLattice;
}

double UnitCellLattice::refineUnitCellScore(void *object)
{
	UnitCellLattice *lat = static_cast<UnitCellLattice *>(object);
	double score = 0;

	for (int i = 0; i < lat->allMillers.size(); i++)
	{
		MillerPtr miller = lat->allMillers[i];
		miller->getMatrix()->recalculateOrientationMatrix();
		miller->recalculateWavelength();

		double wave = miller->getWavelength();
		double mean = miller->getMtzParent()->getWavelength();

		double diff = fabs(mean - wave);
		diff *= 1 / 0.005;

		double contrib = Detector::lookupCache(diff);

		score -= contrib;
	}

	return score;
}

void UnitCellLattice::refineMtzs(std::vector<MtzPtr> newMtzs)
{
	for (int i = 0; i < newMtzs.size(); i++)
	{
		// collect Millers
		std::vector<MillerPtr> millers = newMtzs[i]->strongMillers();
		newMtzs[i]->getMatrix()->setUnitCell(unitCellOnly);

		for (int j = 0; j < millers.size(); j++)
		{
			millers[j]->getMatrix()->setUnitCell(unitCellOnly);
		}

		allMillers.reserve(allMillers.size() + millers.size());
		allMillers.insert(allMillers.end(), millers.begin(), millers.end());
	}

	logged << "Unit cell dimension before: " << printUnitCell() << std::endl;
	sendLog();

	RefinementStrategyPtr strategy = RefinementStrategy::userChosenStrategy();
	strategy->setEvaluationFunction(refineUnitCellScore, this);
	strategy->addParameter(this, getUnitCellA, setUnitCellA, 0.1, 0.0001);
	strategy->addParameter(this, getUnitCellB, setUnitCellB, 0.1, 0.0001);
	strategy->addParameter(this, getUnitCellC, setUnitCellC, 0.1, 0.0001);
	strategy->refine();
	strategy->setVerbose(true);

	logged << "Unit cell dimension after: " << printUnitCell() << std::endl;
	sendLog();


}