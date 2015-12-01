//
//  IndexingSolution.h
//  cppxfel
//
//  Created by Helen Ginn on 18/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__IndexingSolution__
#define __cppxfel__IndexingSolution__

#include <stdio.h>
#include "parameters.h"
#include "SpotVector.h"
#include "Matrix.h"
#include "csymlib.h"
#include "LoggableObject.h"
#include "Holder.h"

typedef std::map<SpotVectorPtr, SpotVectorPtr> SpotVectorMap;
typedef std::map<SpotVectorPtr, MatrixPtr> SpotVectorMatrixMap;
typedef std::map<SpotVectorPtr, SpotVectorMatrixMap > SpotVectorMatrixMap2D;

class IndexingSolution : LoggableObject
{
private:
    static std::vector<SpotVectorPtr> standardVectors;
    static std::vector<MatrixPtr> symOperators;
    static Reflection *newReflection;
    SpotVectorMap spotVectors;
    SpotVectorMatrixMap2D matrices;
    
    bool vectorAgreesWithExistingVectors(SpotVectorPtr observedVector, SpotVectorPtr standardVector);
    static bool vectorMatchesVector(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr *firstMatch, SpotVectorPtr *secondMatch);
    MatrixPtr createSolution(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr firstStandard = SpotVectorPtr());
    bool solutionCompatibleForMerge(IndexingSolutionPtr otherSolution);
    bool vectorPairLooksLikePair(SpotVectorPtr firstObserved, SpotVectorPtr secondObserved, SpotVectorPtr standard1, SpotVectorPtr standard2, MatrixPtr symOperator);
    static bool allVectorMatches(SpotVectorPtr firstVector, SpotVectorPtr secondVector, std::vector<SpotVectorPtr> *firstMatches, std::vector<SpotVectorPtr> *secondMatches);
    void addVectorToList(SpotVectorPtr observedVector, SpotVectorPtr standardVector);
    void addMatrix(SpotVectorPtr observedVector1, SpotVectorPtr observedVector2, MatrixPtr solution);
    bool vectorSolutionsAreCompatible(SpotVectorPtr observedVector, SpotVectorPtr standardVector);
    
    static double distanceTolerance;
    static double angleTolerance;
    static double solutionAngleSpread;
    static int spaceGroupNum;
    static CSym::CCP4SPG *spaceGroup;
    static std::vector<double> unitCell;
    static int maxMillerIndexTrial;
    static MatrixPtr unitCellOnly;
    static MatrixPtr unitCellMatrix;
    static MatrixPtr unitCellMatrixInverse;
    static bool notSetup;
    static bool finishedSetup;
    
    IndexingSolution();
    IndexingSolution(SpotVectorPtr firstVector, SpotVectorPtr secondVector, SpotVectorPtr firstMatch, SpotVectorPtr secondMatch);
    IndexingSolution(SpotVectorMap firstMap, SpotVectorMap secondMap, SpotVectorMatrixMap2D matrixMap1, SpotVectorMatrixMap2D matrixMap2, MatrixPtr symOperator);
    
public:
    IndexingSolutionPtr mergeWithSolution(IndexingSolutionPtr otherSolution);
    static std::vector<IndexingSolutionPtr> startingSolutionsForVectors(SpotVectorPtr firstVector, SpotVectorPtr secondVector);
    int extendFromSpotVectors(std::vector<SpotVectorPtr> *possibleVectors, int limit = 0);
    MatrixPtr createSolution();
    static bool matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2);
    static void setupStandardVectors();
    
    IndexingSolutionPtr copy();
    
    int spotVectorCount()
    {
        return (int)spotVectors.size();
    }
    
    bool wasSuccessful()
    {
        return spotVectors.size();
    }
    
    
    std::string printNetwork();
};

#endif /* defined(__cppxfel__IndexingSolution__) */
