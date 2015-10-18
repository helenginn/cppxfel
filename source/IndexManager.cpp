//
//  IndexManager.cpp
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "IndexManager.h"
#include "FileParser.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "Logger.h"
#include <algorithm>

IndexManager::IndexManager(std::vector<Image *> newImages)
{
    images = newImages;
    
    spaceGroupNum = FileParser::getKey("SPACE_GROUP", 0);
    
    if (spaceGroupNum == 0)
    {
        std::cout << "Please eprovide space group number in SPACE_GROUP" << std::endl;
        exit(1);
    }
    
    spaceGroup = ccp4spg_load_by_ccp4_num(spaceGroupNum);
    
    unitCell = FileParser::getKey("UNIT_CELL", std::vector<double>());
    
    if (unitCell.size() < 6)
    {
        std::cout << "Please supply target unit cell in keyword UNIT_CELL." << std::endl;
        exit(1);
    }
    
    unitCellOnly = Matrix::matrixFromUnitCell(unitCell[0], unitCell[1], unitCell[2],
                                                    unitCell[3], unitCell[4], unitCell[5]);
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    maxDistance = 0;
    
    for (int i = -5; i <= 5; i++)
    {
        for (int j = -5; j <= 5; j++)
        {
            for (int k = -5; k <= 5; k++)
            {
                if (ccp4spg_is_sysabs(spaceGroup, i, j, k))
                    continue;
                
                vec hkl = new_vector(i, j, k);
                vec hkl_transformed = copy_vector(hkl);
                
                if (length_of_vector(hkl) > 4)
                    continue;
                
                unitCellMatrix->multiplyVector(&hkl_transformed);
                
                double distance = length_of_vector(hkl_transformed);

                if (distance > maxDistance)
                    maxDistance = distance;
                
                vectorDistances.push_back(std::make_pair(hkl, distance));
            }
        }
    }
    
    for (int i = 0; i < vectorDistances.size(); i++)
    {
        logged << vectorDistances[i].first.h << "\t"
         << vectorDistances[i].first.k << "\t"
         << vectorDistances[i].first.l << "\t"
        << vectorDistances[i].second << std::endl;
    }
    
    sendLog();
}

bool match_greater_than_match(Match a, Match b)
{
    return a.second > b.second;
}

bool greater_than_scored_matrix(std::pair<MatrixPtr, double> a, std::pair<MatrixPtr, double> b)
{
    return a.second > b.second;
}

void IndexManager::indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset)
{
    double angleTolerance = 1.0 * M_PI / 180;
    double trustTolerance = 5000;
    int maxThreads = FileParser::getMaxThreads();
    
    std::ostringstream logged;
    
    for (int i = offset; i < indexer->images.size(); i += maxThreads)
    {
        indexer->images[i]->compileDistancesFromSpots(indexer->maxDistance);
       
        std::vector<Match> possibleMatches;
        std::vector<MatrixPtr> possibleSolutions;
        
        for (int j = 0; j < indexer->images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr dataPoint = indexer->images[i]->spotVector(j);
            
            for (int k = 0; k < indexer->vectorDistances.size(); k++)
            {
                VectorDistance trial;
                
                trial.first = copy_vector(indexer->vectorDistances[k].first);
                trial.second = indexer->vectorDistances[k].second;
                
                double trust = 1 / (fabs(dataPoint->distance() - trial.second));
                
                if (trust > trustTolerance)
                {
                    std::pair<SpotVectorPtr, VectorDistance> vectorPair = std::make_pair(dataPoint, trial);
                    
                    possibleMatches.push_back(std::make_pair(vectorPair, trust));
                }
            }
        }
        
        std::sort(possibleMatches.begin(), possibleMatches.end(), match_greater_than_match);

        for (int j = 0; j < possibleMatches.size(); j++)
        {
            double trust = possibleMatches[j].second;
            SpotVectorPtr dataPoint = possibleMatches[j].first.first;
            VectorDistance trial = possibleMatches[j].first.second;
            
            logged << trust << "\t" << trial.second << "\t" << dataPoint->distance() << "\t" << trial.first.h << "\t" <<
            trial.first.k << "\t" <<
            trial.first.l << "\t" << std::endl;
        }
        
        Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
     
        logged << "Total spots for " << indexer->images[i]->getFilename() << ": " << indexer->images[i]->spotCount() << std::endl;
        logged << "Total vectors: " << indexer->images[i]->spotVectorCount() << std::endl;
        Logger::mainLogger->addStream(&logged); logged.str("");
        logged << "Total possible vector matches: " << possibleMatches.size() << std::endl;
        Logger::mainLogger->addStream(&logged); logged.str("");
        
        
        for (int j = 0; j < possibleMatches.size(); j++)
        {
            for (int k = j + 1; k < possibleMatches.size(); k++)
            {
                std::pair<SpotVectorPtr, VectorDistance> vectorPair1 = possibleMatches[j].first;
                std::pair<SpotVectorPtr, VectorDistance> vectorPair2 = possibleMatches[k].first;
                
                SpotVectorPtr observed1 = vectorPair1.first;
                SpotVectorPtr observed2 = vectorPair2.first;
                
                vec observedVec1 = observed1->getVector();
                vec observedVec2 = observed2->getVector();
                
                if (!vectors_are_equal(observedVec1, observedVec2))
                {
                    vec simulatedVec1 = possibleMatches[j].first.second.first;
                    vec simulatedVec2 = possibleMatches[k].first.second.first;
                    
                    double angle1 = fabs(angleBetweenVectors(observedVec1, observedVec2));
                    double angle2 = fabs(angleBetweenVectors(simulatedVec1, simulatedVec2));
                    
                    double angleDiff = fabs(angle1 - angle2);
                    
                    if (angleDiff < angleTolerance)
                    {
                        // rotation to get from observedvec1 to observedvec2.
                        MatrixPtr rotateSpotDiffMatrix = rotation_between_vectors(observedVec1, simulatedVec1);
                     
                        vec firstSpotVec = observed1->getFirstSpot()->estimatedVector();
                        vec secondSpotVec = observed1->getSecondSpot()->estimatedVector();
                        vec reverseSpotVec = reverseVector(firstSpotVec);
                        rotateSpotDiffMatrix->multiplyVector(&reverseSpotVec);
                //        vec rereversed = reverseVector(reverseSpotVec);

                        // we have now found one axis responsible for the first observation
                        // now we must rotate around that axis (rereversed) until
                        // observedVec2 matches simulatedVec2
                        // first we make a unit vector to rotate around
                        
                        vec firstAxisUnit = copy_vector(simulatedVec1);
                        scale_vector_to_distance(&firstAxisUnit, 1);
                        
                        // we need to rotate observedVec2 by the first rotation matrix
                        vec rotatedObservedVec2 = copy_vector(observedVec2);
                        rotateSpotDiffMatrix->multiplyVector(&rotatedObservedVec2);
                        
                        // and now we can find the appropriate rotation matrix given our rotation axis firstAxisUnit
                        
                        double resultantAngle = 0;
                        
                        MatrixPtr secondTwizzleMatrix = closest_rotation_matrix(rotatedObservedVec2, simulatedVec2, firstAxisUnit, &resultantAngle);
                        
                        if (resultantAngle > angleTolerance)
                            continue;
                        
                        MatrixPtr combinedMatrix = rotateSpotDiffMatrix->copy();
                        combinedMatrix->multiply(*secondTwizzleMatrix);
                        
                        vec possibleMiller = copy_vector(firstSpotVec);
                        combinedMatrix->multiplyVector(&possibleMiller);
                        indexer->unitCellMatrixInverse->multiplyVector(&possibleMiller);
                        
                        logged << "Possible pair of matches " << j << ", " << k << std::endl;
                        logged << "Relative Millers (" << simulatedVec1.h << ", " << simulatedVec1.k << ", " << simulatedVec1.l << ") and (" << simulatedVec2.h << ", " << simulatedVec2.k << ", " << simulatedVec2.l << ") from trusts " << possibleMatches[j].second << " and " << possibleMatches[k].second << std::endl;

                        Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
                    
                        
                        SpotPtr pair1Spot1 = observed1->getFirstSpot();
                        SpotPtr pair1Spot2 = observed1->getSecondSpot();
                        
                        SpotPtr pair2Spot1 = observed2->getFirstSpot();
                        SpotPtr pair2Spot2 = observed2->getSecondSpot();
                        
                        
                        MatrixPtr rotateFinalMatrix = combinedMatrix->inverse3DMatrix();
                        rotateFinalMatrix->rotate(0, 0, -M_PI/2);
                        
                        MatrixPtr fullMat = MatrixPtr(new Matrix());
                        fullMat->setComplexMatrix(indexer->unitCellOnly, rotateFinalMatrix);
                        
                        possibleSolutions.push_back(fullMat);
                      
                        j++;
                    }
                }
            }
        }
        
        if (possibleSolutions.size() == 0)
            continue;
        
        // dummy Reflection to give us indexing ambiguity info
        
        MillerPtr miller = MillerPtr(new Miller(NULL, 0, 0, 0));
        
        Reflection *newReflection = new Reflection();
        newReflection->setUnitCellDouble(&indexer->unitCell[0]);
        newReflection->setSpaceGroup(indexer->spaceGroupNum);
        
        logged << "Indexing ambiguities for this space group: " << newReflection->ambiguityCount() << std::endl;
        logged << possibleSolutions.size() << " total solutions for image " << indexer->images[i]->getFilename() << std::endl;
        
        std::ostringstream dummyVecStr;
        
        std::vector<std::pair<MatrixPtr, double> > scoredSolutions;
        
        Logger::mainLogger->addStream(&logged); logged.str("");
        
        for (int j = 0; j < possibleSolutions.size(); j++)
        {
            MatrixPtr aMat = possibleSolutions[j];
            vec aDummyVec = new_vector(1, 0, 0);
            aMat->multiplyVector(&aDummyVec);
            double score = 0;
            int count = 0;
            
            for (int k = 0; k < possibleSolutions.size(); k++)
            {
                MatrixPtr bMat = possibleSolutions[k];
                
                for (int l = 0; l < newReflection->ambiguityCount(); l++)
                {
                    MatrixPtr ambiguity = newReflection->matrixForAmbiguity(l);
                    
                    vec bDummyVec = new_vector(1, 0, 0);
                    
                    ambiguity->multiplyVector(&bDummyVec);
                    bMat->multiplyVector(&bDummyVec);
                    
                    double angle = fabs(angleBetweenVectors(aDummyVec, bDummyVec));
                    double tolerance = 8 * M_PI / 180;
                    
                    if (angle < tolerance)
                    {
                        score += tolerance - angle;
                        count++;
                    }
                }
            }
            
            logged << j << "\t" << score << "\t" << count << std::endl;
            
            scoredSolutions.push_back(std::make_pair(aMat, score));
        }
        
        Logger::mainLogger->addStream(&logged, LogLevelDetailed); logged.str("");
        
        std::sort(scoredSolutions.begin(), scoredSolutions.end(), greater_than_scored_matrix);
        
        logged << "Highest scored solution: " << scoredSolutions[0].second << std::endl;
        Logger::mainLogger->addStream(&logged); logged.str("");
        
        indexer->images[i]->setUpIOMRefiner(scoredSolutions[0].first);
        IOMRefinerPtr refiner = indexer->images[i]->getIOMRefiner(0);
        refiner->refineOrientationMatrix();
        mtzSubset->push_back(refiner->newMtz(0));
    }
}

void IndexManager::index()
{
    int maxThreads = FileParser::getMaxThreads();
    
    boost::thread_group threads;
    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(indexThread, this, &managerSubsets[i], i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    int total = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }
    
    mtzs.reserve(total);
    int lastPos = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        mtzs.insert(mtzs.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }
}

void IndexManager::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}
