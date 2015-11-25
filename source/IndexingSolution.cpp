//
//  IndexingSolution.cpp
//  cppxfel
//
//  Created by Helen Ginn on 18/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "IndexingSolution.h"
#include "FileParser.h"
#include "Logger.h"
#include <string>

std::vector<SpotVectorPtr> IndexingSolution::standardVectors;
int IndexingSolution::spaceGroupNum = 0;
CSym::CCP4SPG *IndexingSolution::spaceGroup = NULL;
std::vector<double> IndexingSolution::unitCell;
int IndexingSolution::maxMillerIndexTrial = 4;
MatrixPtr IndexingSolution::unitCellOnly;
MatrixPtr IndexingSolution::unitCellMatrix;
MatrixPtr IndexingSolution::unitCellMatrixInverse;
double IndexingSolution::distanceTolerance = 4000;
double IndexingSolution::angleTolerance = 1.0;
double IndexingSolution::solutionAngleSpread = 8.0 * M_PI / 180;
std::vector<MatrixPtr> IndexingSolution::symOperators;
Reflection *IndexingSolution::newReflection;
bool IndexingSolution::notSetup = true;
bool IndexingSolution::finishedSetup = false;

void IndexingSolution::setupStandardVectors()
{
    notSetup = false;
    
    spaceGroupNum = FileParser::getKey("SPACE_GROUP", 0);
    
    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);
    
    Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup);
    
    unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());
    
    unitCellOnly = Matrix::matrixFromUnitCell(unitCell[0], unitCell[1], unitCell[2],
                                              unitCell[3], unitCell[4], unitCell[5]);
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    if (unitCell.size() < 6)
    {
        std::cout << "Please supply target unit cell in keyword UNIT_CELL." << std::endl;
        exit(1);
    }
    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 4);
    double maxDistance = 0;
    
    for (int i = -maxMillerIndexTrial; i <= maxMillerIndexTrial; i++)
    {
        for (int j = -maxMillerIndexTrial; j <= maxMillerIndexTrial; j++)
        {
            for (int k = -maxMillerIndexTrial; k <= maxMillerIndexTrial; k++)
            {
                if (spaceGroupNum != 19)
                    if (ccp4spg_is_sysabs(spaceGroup, i, j, k))
                        continue;
                
                vec hkl = new_vector(i, j, k);
                vec hkl_transformed = copy_vector(hkl);
                
                if (length_of_vector(hkl) > maxMillerIndexTrial)
                    continue;
                
                unitCellMatrix->multiplyVector(&hkl_transformed);
                
                double distance = length_of_vector(hkl_transformed);
                
                if (distance > maxDistance)
                    maxDistance = distance;
                
                SpotVectorPtr newStandardVector = SpotVectorPtr(new SpotVector(hkl_transformed, hkl));
                
                standardVectors.push_back(newStandardVector);
            }
        }
    }
    
    newReflection = new Reflection();
    newReflection->setUnitCellDouble(&unitCell[0]);
    newReflection->setSpaceGroup(spaceGroupNum);
    
    finishedSetup = true;
}

bool IndexingSolution::matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2)
{
    std::ostringstream logged;
    MatrixPtr mat3 = mat2->copy();
    
    for (int k = 0; k < symOperators.size(); k++)
    {
        MatrixPtr symOperator = symOperators[k];
        
        for (int l = 0; l < newReflection->ambiguityCount(); l++)
        {
            MatrixPtr ambiguity = newReflection->matrixForAmbiguity(l);
            mat3->preMultiply(*symOperator);
            mat3->preMultiply(*ambiguity);
            
            double radianSpread = solutionAngleSpread * M_PI / 180;
            double angle = mat1->similarityToRotationMatrix(mat3, radianSpread);
            
            if (angle < radianSpread && angle != -1)
            {
                return true;
            }
        }
    }
    
    return false;
}

bool IndexingSolution::vectorPairLooksLikePair(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr standard1, SpotVectorPtr standard2, MatrixPtr symOperator)
{
    double firstVectorTrust = firstObserved->trustComparedToStandardVector(standard1);
    double secondVectorTrust = secondObserved->trustComparedToStandardVector(standard2);
    
    if (firstVectorTrust < distanceTolerance)
        return false;
    
    if (secondVectorTrust < distanceTolerance)
        return false;
    
    double standardAngle = standard1->angleWithVector(standard2);
    double realAngle = firstObserved->angleWithVector(secondObserved, symOperator);
    
    double difference = fabs(realAngle - standardAngle);
    
    if (difference > angleTolerance)
        return false;
    
    return true;
}

bool IndexingSolution::allVectorMatches(SpotVectorPtr firstVector, SpotVectorPtr secondVector, std::vector<SpotVectorPtr> *firstMatches, std::vector<SpotVectorPtr> *secondMatches)
{
    bool good = false;
    
    for (int i = 0; i < standardVectors.size(); i++)
    {
        double firstVectorTrust = firstVector->trustComparedToStandardVector(standardVectors[i]);
        
        if (firstVectorTrust > distanceTolerance)
        {
            for (int j = 0; j < standardVectors.size(); j++)
            {
                double secondVectorTrust = secondVector->trustComparedToStandardVector(standardVectors[j]);
                
                if (secondVectorTrust > distanceTolerance)
                {
                    double realAngle = firstVector->angleWithVector(secondVector);
                    double expectedAngle = standardVectors[i]->angleWithVector(standardVectors[j]);
                    
                    double difference = fabs(realAngle - expectedAngle);
                    
                    if (difference < angleTolerance)
                    {
                  /*      std::ostringstream logged;
                        logged << "Pair matched trust, trust, angle diff:\t" << firstVectorTrust << "\t" << secondVectorTrust << "\t" << difference << std::endl;
                        Logger::mainLogger->addStream(&logged, LogLevelDebug);
                    */    
                        firstMatches->push_back(standardVectors[i]);
                        secondMatches->push_back(standardVectors[j]);
                        
                        good = true;
                    }
                }
            }
        }
    }
    
    return false;
}

bool IndexingSolution::vectorMatchesVector(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr *firstMatch, SpotVectorPtr *secondMatch)
{
    for (int i = 0; i < standardVectors.size(); i++)
    {
        double firstVectorTrust = firstVector->trustComparedToStandardVector(standardVectors[i]);
        
        if (firstVectorTrust > distanceTolerance)
        {
            for (int j = 0; j < standardVectors.size(); j++)
            {
                double secondVectorTrust = secondVector->trustComparedToStandardVector(standardVectors[j]);
                
                if (secondVectorTrust > distanceTolerance)
                {
                    double realAngle = firstVector->angleWithVector(secondVector);
                    double expectedAngle = standardVectors[i]->angleWithVector(standardVectors[j]);
                    
                    double difference = fabs(realAngle - expectedAngle);
                    
                    if (difference < angleTolerance)
                    {
                        *firstMatch = standardVectors[i];
                        *secondMatch = standardVectors[j];
                        
                        return true;
                    }
                }
            }
        }
    }
    
    return false;
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
    
    for (int i = 0; i < allSpotVectors.size() - 1; i++)
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
            //   double angle = angleBetweenVectors(aMat, bMat);
            
            if (angle < solutionAngleSpread && angle != -1)
            {
                double addition = pow(solutionAngleSpread - angle, 2);
                score += addition;
                count++;
            }
        }
        
        logged << matrices[j]->summary() << "\t" << score << "\t" << count << std::endl;
        
        scoredSolutions.push_back(std::make_pair(aMat, score));
    }

    sendLog(LogLevelDetailed);

    std::sort(scoredSolutions.begin(), scoredSolutions.end(), match_greater_than_match);
    
    
    MatrixPtr chosenMat = MatrixPtr();
    
    if (scoredSolutions.size())
    {
        chosenMat = scoredSolutions[0].first;
        logged << "Chosen solution:\t" << chosenMat->summary() << "\tfrom " << spotVectors.size() << " vectors"<< std::endl;
        sendLog(LogLevelNormal);
    }
    return chosenMat;
}

MatrixPtr IndexingSolution::createSolution(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved)
{
    SpotVectorPtr firstStandard = spotVectors[firstObserved];
    SpotVectorPtr secondStandard = spotVectors[secondObserved];
    
    vec observedVec1 = firstObserved->getVector();
    vec observedVec2 = secondObserved->getVector();
 
    vec simulatedVec1 = firstStandard->getVector();
    vec simulatedVec2 = secondObserved->getVector();
    
    MatrixPtr rotateSpotDiffMatrix = rotation_between_vectors(observedVec1, simulatedVec1);
    
    vec firstSpotVec = firstObserved->getFirstSpot()->estimatedVector();
    vec reverseSpotVec = reverseVector(firstSpotVec);
    rotateSpotDiffMatrix->multiplyVector(&reverseSpotVec);
    
    vec firstAxisUnit = copy_vector(simulatedVec1);
    scale_vector_to_distance(&firstAxisUnit, 1);
    
    vec rotatedObservedVec2 = copy_vector(observedVec2);
    rotateSpotDiffMatrix->multiplyVector(&rotatedObservedVec2);
    
    double resultantAngle = 0;
    
    MatrixPtr secondTwizzleMatrix = closest_rotation_matrix(rotatedObservedVec2, simulatedVec2, firstAxisUnit, &resultantAngle);
    
    MatrixPtr combinedMatrix = rotateSpotDiffMatrix->copy();
    combinedMatrix->multiply(*secondTwizzleMatrix);
    
    MatrixPtr rotateFinalMatrix = combinedMatrix->inverse3DMatrix();
    rotateFinalMatrix->rotate(0, 0, -M_PI/2);
    
    MatrixPtr fullMat = MatrixPtr(new Matrix());
    fullMat->setComplexMatrix(unitCellOnly, rotateFinalMatrix);
    
    logged << fullMat->summary() << std::endl;
    sendLog(LogLevelDetailed);
    
    
    return fullMat;
}

bool IndexingSolution::solutionCompatibleForMerge(IndexingSolutionPtr otherSolution)
{
    if (spotVectors.size() >= 3 && otherSolution->spotVectors.size() >= 3)
        return true;
    
    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        for (SpotVectorMap::iterator it2 = otherSolution->spotVectors.begin(); it2 != otherSolution->spotVectors.end(); it2++)
        {
            SpotVectorPtr myVector = it->first;
            SpotVectorPtr herVector = it2->first;
            
            if (myVector->hasCommonSpotWithVector(herVector))
            {
                return true;
            }
        }
    }
    
    return false;
}

bool IndexingSolution::vectorAgreesWithExistingVectors(SpotVectorPtr observedVector, SpotVectorPtr standardVector)
{
    for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
    {
        SpotVectorPtr myVector = it->first;
        
        bool similar = vectorPairLooksLikePair(it->first, observedVector, it->second, standardVector, MatrixPtr(new Matrix()));
        
        if (!similar)
            return false;
    }
    
    return true;
}

int IndexingSolution::extendFromSpotVectors(std::vector<SpotVectorPtr> *possibleVectors)
{
    int added = 0;
    
    for (int i = 0; i < possibleVectors->size() && i < 3000; i++)
    {
        SpotVectorPtr herVector = (*possibleVectors)[i];
        bool commonSpots = false;
        bool duplicates = false;
        
        for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
        {
            SpotVectorPtr myVector = it->first;
            
            if (myVector->hasCommonSpotWithVector(herVector))
            {
                commonSpots = true;
            }
            
            if (myVector == herVector)
            {
                duplicates = true;
            }
        }
        
        if (!commonSpots || duplicates)
            continue;
        
        sendLog(LogLevelDebug);
        
        for (int j = 0; j < standardVectors.size(); j++)
        {
            bool agrees = vectorAgreesWithExistingVectors(herVector, standardVectors[j]);
            
            if (agrees)
            {
                added++;
                spotVectors[herVector] = standardVectors[j];
                possibleVectors->erase(possibleVectors->begin() + i);
                i--;
            }
        }
    }
    
    return added;
}

IndexingSolutionPtr IndexingSolution::mergeWithSolution(IndexingSolutionPtr otherSolution)
{
 /*   if (!solutionCompatibleForMerge(otherSolution))
    {
        return IndexingSolutionPtr();
    }
   */ 
    bool hadCommonSpot = false;
    
    for (int i = 0; i < newReflection->ambiguityCount(); i++)
    {
        MatrixPtr ambiguity = newReflection->matrixForAmbiguity(i);
        
        for (int j = 0; j < symOperators.size(); j++)
        {
            MatrixPtr symOperator = symOperators[j]->copy();
            
            symOperator->multiply(*ambiguity);
            
            for (SpotVectorMap::iterator it = spotVectors.begin(); it != spotVectors.end(); it++)
            {
                bool allGood = true;
                
                for (SpotVectorMap::iterator it2 = otherSolution->spotVectors.begin(); it2 != otherSolution->spotVectors.end() && allGood; it2++)
                {
                    SpotVectorPtr myVector = it->first;
                    SpotVectorPtr herVector = it2->first;
                    
                    if (myVector->hasCommonSpotWithVector(herVector))
                    {
                        hadCommonSpot = true;
                    }
                    
                    if (!vectorPairLooksLikePair(myVector, herVector, it->second, it2->second, symOperator))
                    {
                        allGood = false;
                        break;
                    }
                }
                
                if (!hadCommonSpot && spotVectors.size() < 3 && otherSolution->spotVectorCount() < 3)
                    continue;
                
                if (allGood)
                {
                    IndexingSolutionPtr newSolution = IndexingSolutionPtr(new IndexingSolution(spotVectors, otherSolution->spotVectors, matrices, otherSolution->matrices, symOperator));
                    
                    logged << "Merging solution (" << spotVectors.size()  << " + " << otherSolution->spotVectors.size() << ") -> " << newSolution->spotVectorCount() << " sym " << j << std::endl;
                    sendLog(LogLevelDebug);
                    
                    return newSolution;
                }
            }
        }
    }
    
    return IndexingSolutionPtr();
}

IndexingSolution::IndexingSolution(SpotVectorMap firstMap, SpotVectorMap secondMap, SpotVectorMatrixMap2D matrixMap1, SpotVectorMatrixMap2D matrixMap2, MatrixPtr symOperator)
{
    distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.0);
    angleTolerance = FileParser::getKey("MINIMUM_TRUST_ANGLE", 1.0) * M_PI / 180;
    solutionAngleSpread = FileParser::getKey("SOLUTION_ANGLE_SPREAD", 8.0) * M_PI / 180;
    
    spotVectors = firstMap;
    
    for (SpotVectorMap::iterator it = secondMap.begin(); it != secondMap.end(); it++)
    {
        SpotVectorPtr myVector = it->first;
        spotVectors[myVector] = it->second->vectorRotatedByMatrix(symOperator);
    }
    
    matrices = matrixMap1;
    /*
    for (SpotVectorMatrixMap2D::iterator it = matrixMap2.begin(); it != matrixMap2.end(); it++)
    {
        SpotVectorPtr firstVec = it->first;
        SpotVectorMatrixMap internalMap = it->second;
        
        for (SpotVectorMatrixMap::iterator it2 = internalMap.begin(); it2 != internalMap.end(); it2++)
        {
            SpotVectorPtr secondVec = it2->first;
            MatrixPtr matrix = it2->second;
            
            matrices[firstVec][secondVec] = matrix;
        }
    }*/
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

IndexingSolution::IndexingSolution(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr firstMatch, SpotVectorPtr secondMatch)
{
    distanceTolerance = FileParser::getKey("MINIMUM_TRUST_DISTANCE", 4000.0);
    angleTolerance = FileParser::getKey("MINIMUM_TRUST_ANGLE", 1.0) * M_PI / 180;
    
    
    spotVectors[firstVector] = firstMatch;
    spotVectors[secondVector] = secondMatch;
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
    
    /*std::vector<SpotVectorPtr> firstMatches, secondMatches;
    
    allVectorMatches(firstVector, secondVector, &firstMatches, &secondMatches);
    
    if (firstMatches.size() == 0)
        return std::vector<IndexingSolutionPtr>();
    */
    
    SpotVectorPtr firstMatch = SpotVectorPtr();
    SpotVectorPtr secondMatch = SpotVectorPtr();
    
    vectorMatchesVector(firstVector, secondVector, &firstMatch, &secondMatch);
    
    if (!(firstMatch && secondMatch))
        return solutions;
    
    IndexingSolutionPtr newSolution = IndexingSolutionPtr(new IndexingSolution(firstVector, secondVector, firstMatch, secondMatch));
    solutions.push_back(newSolution);
    
    return solutions;
}