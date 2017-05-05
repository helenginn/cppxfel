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

class IndexManager : LoggableObject, public boost::enable_shared_from_this<IndexManager>
{
protected:
	struct SpotVectorPair
	{
		SpotVectorPtr vec1;
		SpotVectorPtr vec2;
	};

    UnitCellLatticePtr lattice;
    std::vector<ImagePtr> images;
    std::vector<ImagePtr> mergeImages;
	std::vector<SpotPtr> goodSpots;
    std::vector<SpotVectorPtr> goodVectors;
	std::vector<SpotVectorPair> goodVectorPairs;
    DetectorWeakPtr _activeDetector;
    double interPanelDistance;
    double intraPanelDistance;
    int spaceGroupNum;
    std::vector<MtzPtr> mtzs;
    double _maxFrequency;
    ImagePtr getNextImage();
    int nextImage;
    std::mutex indexMutex;
    PseudoScoreType scoreType;
    double proportionDistance;
    CSVPtr angleCSV;
    bool _canLockVectors;
    static int _cycleNum;
    PseudoScoreWeightingAxis _axisWeighting;
    
    double maxMillerIndexTrial;
    double maxDistance;
    double smallestDistance;
    double minReciprocalDistance;
    PowderHistogram generatePowderHistogram(int intraPanel = -1, int perfectPadding = 0);
    std::vector<VectorDistance> vectorDistances;
    std::vector<MtzPtr> consolidateOrientations(ImagePtr image1, ImagePtr image2, int *oneHand, int *otherHand, int *both);
    PseudoScoreType checkVectors(SpotVectorPtr vec1, SpotVectorPtr vec2);
    bool checkVector(SpotVectorPtr spotVector, bool permissive = false);
    bool processVector(SpotVectorPtr vec, double *score, double *count, bool lock, bool skipCheck);
    
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
    
    void setActiveDetector(DetectorPtr detector, GeometryScoreType type)
    {
        _activeDetector = detector;
        detector->setIndexManager(shared_from_this(), type);
		goodVectorPairs.clear();
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
    
	static double pseudoAngleScore(void *object);
	static double pseudoDistanceScore(void *object, bool writeToCSV = false, std::string tag = "");
    void combineLists();
    static void indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset);
    void index();
	bool writtenPDB;
	static double staticPseudoAnglePDB(void *object);
	CSVPtr pseudoAnglePDB(bool writePost = false);
    void powderPattern(std::string csvName = "powder.csv", bool force = true);
    
    void lockVectors()
    {
        _canLockVectors = true;
    }
    
    void setAxisWeighting(PseudoScoreWeightingAxis weighting)
    {
        _axisWeighting = weighting;
    }
    
    static void setCycleNum(int num)
    {
        _cycleNum = num;
    }
    
    static int getCycleNum()
    {
        return _cycleNum;
    }
    
    void resetMaxFrequency()
    {
        _maxFrequency = -1;
    }
    
    void clearGoodVectors()
    {
        goodVectors.clear();
		goodVectorPairs.clear();
    }
    
    void plotGoodVectors();
    static double pseudoScore(void *object);
    static double debugPseudoScore(void *object);
    IndexManager(std::vector<ImagePtr>images);
};


#endif /* defined(__cppxfel__IndexManager__) */
