//
//  IndexManager.h
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__IndexManager__
#define __cppxfel__IndexManager__

#include "cmtzlib.h"
#include "csymlib.h"
#include <stdio.h>
#include <vector>
#include <map>
#include "Image.h"
#include "Vector.h"
#include "parameters.h"
#include "Logger.h"
#include "Reflection.h"
#include "LoggableObject.h"
#include <mutex>
#include "Detector.h"

class IndexManager : LoggableObject
{
protected:
    UnitCellLatticePtr lattice;
    std::vector<ImagePtr> images;
    std::vector<ImagePtr> mergeImages;
    std::vector<double> unitCell;
    std::vector<MatrixPtr> symOperators;
    std::vector<SpotVectorPtr> goodVectors;
    MatrixPtr unitCellOnly;
    MatrixPtr unitCellMatrix;
    MatrixPtr unitCellMatrixInverse;
    Reflection *newReflection;
    CSym::CCP4SPG *spaceGroup;
    DetectorWeakPtr _activeDetector;
    int spaceGroupNum;
    std::vector<MtzPtr> mtzs;
    double minimumTrustDistance;
    double minimumTrustAngle;
    double solutionAngleSpread;
    double lastTime;
    ImagePtr getNextImage();
    int nextImage;
    std::mutex indexMutex;
    bool modifyParameters();
    PseudoScoreType scoreType;
    double proportionDistance;
    static double pseudoDistanceScore(void *object);
    static double pseudoAngleScore(void *object);
    CSVPtr angleCSV;
    CSVPtr angleConsistencyCSV;
    bool _canLockVectors;
    
    void updateAllSpots();
    bool matrixSimilarToMatrix(MatrixPtr mat1, MatrixPtr mat2);
    int indexOneImage(ImagePtr image, std::vector<MtzPtr> *mtzSubset);
    double maxMillerIndexTrial;
    double maxDistance;
    double smallestDistance;
    double minReciprocalDistance;
    PowderHistogram generatePowderHistogram(int intraPanel = -1, int perfectPadding = 0);
    std::vector<VectorDistance> vectorDistances;
    std::vector<IOMRefinerPtr> consolidateOrientations(ImagePtr image1, ImagePtr image2, int *oneHand, int *otherHand, int *both);
    PseudoScoreType checkVectors(SpotVectorPtr vec1, SpotVectorPtr vec2);
    bool checkVector(SpotVectorPtr spotVector, bool permissive = false);
    void processVector(SpotVectorPtr vec, double *score, double *count, bool lock, bool skipCheck);
    
    DetectorPtr getActiveDetector()
    {
        DetectorPtr det = _activeDetector.lock();
        if (!det)
        {
            det = Detector::getMaster();
        }
        
        return det;
    }
    
public:
    ImagePtr getImage(int i)
    {
        return images[i];
    }
    
    std::vector<MtzPtr> getMtzs()
    {
        return mtzs;
    }
    
    void setMergeImages(std::vector<ImagePtr> otherImages)
    {
        mergeImages = otherImages;
    }
    
    void setActiveDetector(DetectorPtr detector)
    {
        _activeDetector = detector;
    }
    
    void setPseudoScoreType(PseudoScoreType type)
    {
        scoreType = type;
    }
    
    PseudoScoreType getPseudoScoreType()
    {
        return scoreType;
    }
    
    UnitCellLatticePtr getLattice()
    {
        return lattice;
    }
    
    void setProportionDistance(double newProp)
    {
        proportionDistance = newProp;
    }
    
    void combineLists();
    void indexingParameterAnalysis();
    static void indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset);
    void index();
    void pseudoAngleCSV();
    void powderPattern(std::string csvName = "powder.csv", bool force = true);
    void refineUnitCell();
    
    void lockVectors()
    {
        _canLockVectors = true;
    }
    
    static double pseudoScore(void *object);
    static double debugPseudoScore(void *object);
    IndexManager(std::vector<ImagePtr>images);
};


#endif /* defined(__cppxfel__IndexManager__) */
