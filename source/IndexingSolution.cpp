//
//  IndexingSolution.cpp
//  cppxfel
//
//  Created by Helen Ginn on 18/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include <float.h>
#include "IndexingSolution.h"
#include "FileParser.h"
#include "Logger.h"
#include <string>
#include "misc.h"
#include <iomanip>
#include "UnitCellLattice.h"
#include "Reflection.h"

int IndexingSolution::spaceGroupNum = 0;
CSym::CCP4SPG *IndexingSolution::spaceGroup = NULL;
std::vector<double> IndexingSolution::unitCell;
int IndexingSolution::maxMillerIndexTrial = 4;
MatrixPtr IndexingSolution::unitCellOnly;
MatrixPtr IndexingSolution::unitCellMatrix;
MatrixPtr IndexingSolution::unitCellMatrixInverse;
double IndexingSolution::distanceTolerance = 4000;
double IndexingSolution::distanceToleranceReciprocal = 0.0025;
double IndexingSolution::angleTolerance = 1.0 * M_PI / 180;
double IndexingSolution::solutionAngleSpread = 8.0 * M_PI / 180;
double IndexingSolution::approximateCosineDelta = 0.017;
UnitCellLatticePtr IndexingSolution::lattice;
Reflection *IndexingSolution::newReflection;
bool IndexingSolution::notSetup = true;
bool IndexingSolution::finishedSetup = false;
bool IndexingSolution::checkingCommonSpots = true;
std::mutex IndexingSolution::setupMutex;


void IndexingSolution::setupStandardVectors()
{
    if (finishedSetup) return;

    setupMutex.lock();

    if (finishedSetup)
    {
        setupMutex.unlock();
        return;
    }

    notSetup = false;

    distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.0);
    distanceToleranceReciprocal = 1 / distanceTolerance;
    angleTolerance = FileParser::getKey("MINIMUM_TRUST_ANGLE", 1.0) * M_PI / 180;
    solutionAngleSpread = FileParser::getKey("SOLUTION_ANGLE_SPREAD", 8.0) * M_PI / 180;
    checkingCommonSpots = FileParser::getKey("CHECKING_COMMON_SPOTS", true);

    double fortyFiveAngle = M_PI / 4;
    double slightlySmallerAngle = fortyFiveAngle - angleTolerance;

    approximateCosineDelta = fabs(cos(fortyFiveAngle) - cos(slightlySmallerAngle));

    spaceGroupNum = FileParser::getKey("SPACE_GROUP", 0);

    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);

    unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());

    if (unitCell.size() < 6)
    {
        std::cout << "Please supply target unit cell in keyword UNIT_CELL." << std::endl;
        exit(1);
    }

    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 4);

    lattice = UnitCellLattice::getMainLattice();

    newReflection = new Reflection();
    newReflection->setUnitCell(unitCell);
    newReflection->setSpaceGroup(spaceGroupNum);

    finishedSetup = true;

    setupMutex.unlock();
}

void IndexingSolution::calculateSimilarStandardVectorsForImageVectors(std::vector<SpotVectorPtr> vectors)
{
    for (int i = 0; i < vectors.size(); i++)
    {
        vectors[i]->addSimilarLengthStandardVectors(lattice->getStandardVectors(), distanceTolerance);
    }
}

bool IndexingSolution::matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2, bool force)
{
    double minTrace = FLT_MAX;

    std::ostringstream logged;


    for (int k = 0; k < symOperatorCount(); k++)
    {
        MatrixPtr symOp = symOperator(k);

        for (int l = 0; l < newReflection->ambiguityCount(); l++)
        {
            MatrixPtr mat3 = mat2->copy();
            MatrixPtr ambiguity = newReflection->matrixForAmbiguity(l);
            mat3->getRotation()->preMultiply(*symOp);
            mat3->getRotation()->preMultiply(*ambiguity);

            MatrixPtr subtractedMat = mat1->getRotation()->copy();
            subtractedMat->subtract(mat3->getRotation());
            MatrixPtr transposedMat = subtractedMat->copy()->transpose();
            transposedMat->multiply(*subtractedMat);

            double trace = transposedMat->trace();

            if (trace < minTrace)
            {
                minTrace = trace;
             //   logged << "Trace with existing solution: " << trace << std::endl;
             //   Logger::log(logged);
            }
        }
    }

    double badMinTrace = sqrt(4 * (1 - cos(solutionAngleSpread)));

    return (minTrace < badMinTrace);
}

bool IndexingSolution::vectorPairLooksLikePair(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr standard1, SpotVectorPtr standard2)
{
    double standardCos = standard1->cosineWithVector(standard2);
    double realCos = firstObserved->cosineWithVector(secondObserved);

    double difference = fabs(realCos - standardCos);

    return (difference < approximateCosineDelta);
}

typedef std::map<double, std::pair<SpotVectorPtr, SpotVectorPtr> > MatchPair;

bool IndexingSolution::vectorMatchesVector(SpotVectorPtr firstVector, SpotVectorPtr secondVector, std::vector<SpotVectorPtr> *firstMatches, std::vector<SpotVectorPtr> *secondMatches)
{
    MatchPair matchPairs;
    double realAngle = firstVector->angleWithVector(secondVector);

    for (int i = 0; i < uniqueSymVectorCount(); i++) // FIXME: standardVectorCount
    {
        double firstVectorTrust = firstVector->trustComparedToStandardVector(uniqueSymVector(i));
        double firstTolerance = firstVector->getMinDistanceTolerance();

        if (firstVectorTrust > firstTolerance)
        {
            for (int j = 0; j < standardVectorCount(); j++)
            {
                double secondVectorTrust = secondVector->trustComparedToStandardVector(standardVector(j));
                double secondTolerance = secondVector->getMinDistanceTolerance();

                if (secondVectorTrust > secondTolerance)
                {
                    double expectedAngle = uniqueSymVector(i)->angleWithVector(standardVector(j));

                                        double degrees = expectedAngle * 180 / M_PI;

                                        if (degrees < 5.0 || degrees > 180 - 5.0)
                                        {
                                                continue;
                                        }

                    double difference = fabs(realAngle - expectedAngle);

                    if (difference < angleTolerance)
                    {
                        std::pair<SpotVectorPtr, SpotVectorPtr> aPair = std::make_pair(uniqueSymVector(i), standardVector(j));
                        double score = 1 / firstVectorTrust + 1 / secondVectorTrust;

                        matchPairs[score] = aPair;
                        firstMatches->push_back(uniqueSymVector(i));
                        secondMatches->push_back(standardVector(j));
//                        return true;
                    }
                }
            }
        }
    }

    for (MatchPair::iterator it = matchPairs.begin(); it != matchPairs.end(); it++)
    {
        firstMatches->push_back(it->second.first);
        secondMatches->push_back(it->second.second);
    }

    return true;
}

bool match_greater_than_match(std::pair<MatrixPtr, double> a, std::pair<MatrixPtr, double> b)
{
    return a.second > b.second;
}

MatrixPtr IndexingSolution::createSolution()
{
    std::vector<SpotVectorPtr> allSpotVectors;
    std::vector<MatrixPtr> matrices;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        allSpotVectors.push_back(it->first);
    }

    for (int i = 0; i < allSpotVectors.size() - 1 && i < 1; i++)
    {
        for (int j = i + 1; j < allSpotVectors.size(); j++)
        {
            MatrixPtr mat = createSolution(allSpotVectors[i], allSpotVectors[j]);
            matrices.push_back(mat);
        }
    }

    std::vector<std::pair<MatrixPtr, double> > scoredSolutions;
    for (int j = 0; j < matrices.size(); j++)
    {
        MatrixPtr aMat = matrices[j];
        double score = 0;
        int count = 0;

        for (int k = 0; k < matrices.size(); k++)
        {
            MatrixPtr bMat = matrices[k]->copy();

            double angle = bMat->similarityToRotationMatrix(aMat, solutionAngleSpread);

            if (angle < solutionAngleSpread && angle != -1)
            {
                double addition = pow(solutionAngleSpread - angle, 2);
                score += addition;
                count++;
            }
        }

#ifndef __OPTIMIZE__
        logged << matrices[j]->summary() << "\t" << score << "\t" << count << std::endl;
#endif

        scoredSolutions.push_back(std::make_pair(aMat, score));
    }

#ifndef __OPTIMIZE__
    sendLog(LogLevelDebug);
#endif

    std::sort(scoredSolutions.begin(), scoredSolutions.end(), match_greater_than_match);


    MatrixPtr chosenMat = MatrixPtr(new Matrix());

    if (scoredSolutions.size())
    {
        chosenMat = scoredSolutions[0].first;
#ifndef __OPTIMIZE__
        logged << "Chosen solution:\t" << chosenMat->summary() << "\tfrom " << spotVectors.size() << " vectors"<< std::endl;
        sendLog(LogLevelDebug);

        logged << chosenMat->description() << std::endl;
        sendLog(LogLevelDebug);
#endif
    }

    if (modMatrix)
    {
        chosenMat->getRotation()->preMultiply(*modMatrix);
        chosenMat->recalculateOrientationMatrix();
    }

    return chosenMat;
}

MatrixPtr IndexingSolution::createSolution(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr firstStandard)
{
    if (!firstStandard)
    {
        firstStandard = spotVectors[firstObserved];
    }

    SpotVectorPtr secondStandard = spotVectors[secondObserved];

    // Get our important variables.
    vec observedVec1 = firstObserved->getUnitVector();
    vec observedVec2 = secondObserved->getUnitVector();

    vec simulatedVec1 = firstStandard->getUnitVector();
    vec simulatedVec2 = secondStandard->getUnitVector();

    // Rotate reciprocal space so that the first simulated vector lines up with the first observed vector.
    MatrixPtr rotateSpotDiffMatrix = rotation_between_vectors(simulatedVec1, observedVec1);

    // Rotate this vector by first matrix so that we have treated the two simulated vectors equally
    rotateSpotDiffMatrix->multiplyVector(&simulatedVec2); // checked

    double resultantAngle = 0;

    // Now we twirl around the firstAxisUnit until the rotated simulated vector matches the second observed vector
    // as closely as possible.

    MatrixPtr secondTwizzleMatrix = closest_rotmat_analytical(simulatedVec2, observedVec2, observedVec1, &resultantAngle);

    // We want to apply the first matrix and then the second matrix, so we multiply these.
    rotateSpotDiffMatrix->multiply(*secondTwizzleMatrix);

    // Create the goods.
    MatrixPtr fullMat = MatrixPtr(new Matrix());
    fullMat->setComplexMatrix(lattice->getUnitCellOnly()->copy(), rotateSpotDiffMatrix);

    // Send back the goods.
    return fullMat;
}

bool IndexingSolution::spotVectorHasAnAppropriateDistance(SpotVectorPtr observedVector)
{
    double myDistance = observedVector->distance();
    double myTolerance = observedVector->getMinDistanceTolerance();

    for (int i = 0; i < standardVectorCount(); i++)
    {
        double distance = standardVector(i)->distance();

        if (fabs(distance - myDistance) < 1 / myTolerance)
        {
            return true;
        }
    }

    return false;
}

void IndexingSolution::pruneSpotVectors(std::vector<SpotVectorPtr> *spotVectors)
{
    int count = 0;

    for (int i = 0; i < spotVectors->size(); i++)
    {
        bool appropriate = spotVectorHasAnAppropriateDistance((*spotVectors)[i]);

        if (!appropriate)
        {
            count++;
            spotVectors->erase(spotVectors->begin() + i);
            i--;
        }
    }
}

void IndexingSolution::removeSpotVectors(std::vector<SpotVectorPtr> *spotVectors)
{/*
    int count = 0;

    for (int i = 0; i < spotVectors->size(); i++)
    {
        bool exists = this->spotVectors.count((*spotVectors)[i]);

        if (exists)
        {
            count++;
            spotVectors->erase(spotVectors->begin() + i);
            i--;
        }
    }

    logged << "Removed " << count << " spot vectors which led to a bad solution." << std::endl;
    sendLog();*/
}

bool IndexingSolution::spotsAreNotTooClose(SpotVectorPtr observedVector)
{
    SpotPtr spot1 = observedVector->getFirstSpot();
    SpotPtr spot2 = observedVector->getSecondSpot();

    double squareMinDistance = pow(lattice->getMinDistance(), 2) * 0.7;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        SpotPtr spot3 = vector->getFirstSpot();
        SpotPtr spot4 = vector->getSecondSpot();

        if (spot3->closeToSecondSpot(spot1, squareMinDistance) || spot3->closeToSecondSpot(spot2, squareMinDistance)
            || spot4->closeToSecondSpot(spot1, squareMinDistance) || spot4->closeToSecondSpot(spot2, squareMinDistance))
            return false;
    }

    return true;
}

bool IndexingSolution::vectorAgreesWithExistingVectors(SpotVectorPtr observedVector, SpotVectorPtr standardVector)
{
    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        bool similar = vectorPairLooksLikePair(it->first, observedVector, it->second, standardVector);

        if (!similar)
            return false;
    }

    return true;
}

bool IndexingSolution::vectorSolutionsAreCompatible(SpotVectorPtr observedVector, SpotVectorPtr standardVector)
{
    int count = 0;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr myVector = it->first;

        MatrixPtr newSolution = createSolution(observedVector, myVector, standardVector);

        for (SpotVectorMatrixMap2D::iterator jt = matrices.begin(); jt != matrices.end() && count < 5; jt++)
        {
            for (SpotVectorMatrixMap::iterator kt = jt->second.begin(); kt != jt->second.end(); kt++)
            {
                MatrixPtr mat = kt->second;

                if (!IndexingSolution::matrixSimilarToMatrix(newSolution, mat))
                {
                    return false;
                }
            }

            count++;
        }
    }

    return true;
}

void IndexingSolution::addMatrix(SpotVectorPtr observedVector1, SpotVectorPtr observedVector2, MatrixPtr solution)
{
    matrices[observedVector1][observedVector2] = solution;
    matrices[observedVector2][observedVector1] = solution;

    double theta, phi, psi;
    solution->eulerAngles(&theta, &phi, &psi);

    double allThetas = averageTheta * matrixCount;
    double allPhis = averagePhi * matrixCount;
    double allPsis = averagePsi * matrixCount;

    matrixCount++;
    allThetas += theta;
    allPhis += phi;
    allPsis += psi;

    allThetas /= matrixCount;
    allPhis /= matrixCount;
    allPsis /= matrixCount;

    averageTheta = allThetas;
    averagePhi = allPhis;
    averagePsi = allPsis;
}

void IndexingSolution::addVectorToList(SpotVectorPtr observedVector, SpotVectorPtr standardVector)
{
    spotVectors[observedVector] = standardVector;

    if (spotVectors.size() == 1)
        return;

    for (SpotVectorMap::iterator i = spotVectors.begin(); i != spotVectors.end(); i++)
    {
        if (i->first == observedVector)
            continue;

        MatrixPtr newSolution = createSolution(observedVector, i->first);
        addMatrix(i->first, observedVector, newSolution);
    }
}


int IndexingSolution::extendFromSpotVectors(std::vector<SpotVectorPtr> *possibleVectors, int limit)
{
    int added = 0;

    for (int i = lastUsed; i < possibleVectors->size() && i < lastUsed + 15000; i++)
    {
        SpotVectorPtr possibleVector = (*possibleVectors)[i];
        bool commonSpots = false;
        bool duplicates = false;

        for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
        {
            SpotVectorPtr myVector = it->first;

            if (checkingCommonSpots && myVector->hasCommonSpotWithVector(possibleVector))
            {
                commonSpots = true;
            }

            if (myVector == possibleVector)
            {
                duplicates = true;
                break;
            }
        }

        if (checkingCommonSpots && !commonSpots)
            continue;

        if (duplicates)
            continue;

        // Find a standard vector which works with all the other vectors
        std::vector<SpotVectorPtr> standardSelection = possibleVector->standardVectorsOfSameDistance();

        for (int j = 0; j < standardSelection.size(); j++)
        {
            bool agrees = vectorAgreesWithExistingVectors(possibleVector, standardSelection[j]);

            if (!agrees)
                continue;

            agrees = vectorSolutionsAreCompatible(possibleVector, standardSelection[j]);

            if (!agrees)
            {
                continue;
            }

            agrees = spotsAreNotTooClose(possibleVector);

            if (!agrees)
            {
                logged << "Too close!" << std::endl;
                sendLog(LogLevelDetailed);

                continue;
            }

            if (agrees)
            {
                added++;

                addVectorToList(possibleVector, standardSelection[j]);

                possibleVectors->erase(possibleVectors->begin() + i);
                i--;
                lastUsed = i;

                if (added == limit && limit > 0)
                    return added;
            }
        }
    }

    return added;
}

IndexingSolution::IndexingSolution(SpotVectorMap firstMap, SpotVectorMap secondMap, SpotVectorMatrixMap2D matrixMap1, SpotVectorMatrixMap2D matrixMap2, MatrixPtr symOperator)
{
    spotVectors = firstMap;

    for (SpotVectorMap::iterator it = secondMap.begin(); it != secondMap.end(); it++)
    {
        SpotVectorPtr myVector = it->first;
        spotVectors[myVector] = it->second->vectorRotatedByMatrix(symOperator);
    }

    matrices = matrixMap1;
}

std::string IndexingSolution::getNetworkPDB()
{
    std::ostringstream pdbLog;
    int count = 0;

    pdbLog << "HETATM";
    pdbLog << std::fixed;
    pdbLog << std::setw(5) << count << "                   ";
    pdbLog << std::setprecision(2) << std::setw(8)  << 0;
    pdbLog << std::setprecision(2) << std::setw(8) << 0;
    pdbLog << std::setprecision(2) << std::setw(8) << 0;
    pdbLog << "                       N" << std::endl;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        count++;
        SpotVectorPtr standard = it->second;

        pdbLog << "HETATM";
        pdbLog << std::fixed;
        pdbLog << std::setw(5) << count << "                   ";
        pdbLog << std::setprecision(2) << std::setw(8)  << standard->getH();
        pdbLog << std::setprecision(2) << std::setw(8) << standard->getK();
        pdbLog << std::setprecision(2) << std::setw(8) << standard->getL();
        pdbLog << "                       O" << std::endl;
    }

    return pdbLog.str();
}

std::string IndexingSolution::printNetwork()
{
    std::ostringstream printed;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        SpotVectorPtr standard = it->second;

        double x1 = vector->getFirstSpot()->getX();
        double y1 = vector->getFirstSpot()->getY();

        double x2 = vector->getSecondSpot()->getX();
        double y2 = vector->getSecondSpot()->getY();

        double h = standard->getH();
        double k = standard->getK();
        double l = standard->getL();

        printed << "vec\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << "\t" << h << "\t" << k << "\t" << l << std::endl;
    }

    return printed.str();
}

IndexingSolutionPtr IndexingSolution::copy()
{
    IndexingSolutionPtr newPtr = IndexingSolutionPtr(new IndexingSolution());

    newPtr->distanceTolerance = distanceTolerance;
    newPtr->angleTolerance = angleTolerance;
    newPtr->spotVectors = spotVectors;
    newPtr->matrices = matrices;
    newPtr->averageTheta = averageTheta;
    newPtr->averagePhi = averagePhi;
    newPtr->averagePsi = averagePsi;
    newPtr->lastUsed = lastUsed;

    return newPtr;
}

IndexingSolution::IndexingSolution()
{

}

IndexingSolution::IndexingSolution(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr firstMatch, SpotVectorPtr secondMatch)
{
    averageTheta = 0;
    averagePhi = 0;
    averagePsi = 0;
    matrixCount = 0;
    lastUsed = 0;

    addVectorToList(firstVector, firstMatch);
    addVectorToList(secondVector, secondMatch);
}

std::vector<IndexingSolutionPtr> IndexingSolution::startingSolutionsForVectors(SpotVectorPtr firstVector, SpotVectorPtr secondVector)
{
    if (notSetup)
    {
        setupStandardVectors();
    }

    if (!firstVector->hasCommonSpotWithVector(secondVector))
        return std::vector<IndexingSolutionPtr>();

    std::vector<IndexingSolutionPtr> solutions;

    std::vector<SpotVectorPtr> firstMatches, secondMatches;

    vectorMatchesVector(firstVector, secondVector, &firstMatches, &secondMatches);
    std::ostringstream logged;

    if (!firstMatches.size())
    {
        return solutions;
    }
    else
    {
        std::ostringstream logged;
        logged << "Seed matches " << firstMatches.size() << " combinations." << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelDetailed);
    }

    for (int i = 0; i < firstMatches.size(); i++)
    {
        IndexingSolutionPtr newSolution = IndexingSolutionPtr(new IndexingSolution(firstVector, secondVector, firstMatches[i], secondMatches[i]));
        solutions.push_back(newSolution);
    }

    return solutions;
}

IndexingSolution::~IndexingSolution()
{
    spotVectors.clear();
    SpotVectorMap().swap(spotVectors);
}

std::vector<double> IndexingSolution::totalDistances()
{
    std::vector<double> distances;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        //   SpotVectorPtr standard = it->second;

        distances.push_back(vector->distance());
    }

    return distances;
}

std::vector<double> IndexingSolution::totalAngles()
{
    std::vector<double> angles;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector1 = it->first;
        SpotVectorPtr standard1 = it->second;

        for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
        {
            SpotVectorPtr vector2 = it->first;
            SpotVectorPtr standard2 = it->second;

            double angle1 = standard2->angleWithVector(standard1);
            double angle2 = vector2->angleWithVector(vector1);

            double angleDiff = fabs(angle2 - angle1);

            angles.push_back(angleDiff * 180 / M_PI);
        }
    }

    return angles;
}

std::vector<double> IndexingSolution::totalDistanceTrusts()
{
    std::vector<double> trusts;

    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr vector = it->first;
        SpotVectorPtr standard = it->second;

        double distanceDiff = fabs(vector->distance() - standard->distance());
        double trust = 1 / distanceDiff;

        if (trust > distanceTolerance * 10)
            trust = distanceTolerance * 10;

        trusts.push_back(trust);
    }

    return trusts;
}

double IndexingSolution::getMinDistance()
{
    setupStandardVectors();
    return lattice->getMinDistance();
}
