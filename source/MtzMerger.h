//
//  MtzMerger.h
//  cppxfel
//
//  Created by Helen Ginn on 09/05/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__MtzMerger__
#define __cppxfel__MtzMerger__

#include <stdio.h>
#include "parameters.h"
#include "LoggableObject.h"
#include <mutex>

typedef enum
{
    MtzRejectionNotRejected,
    MtzRejectionCorrelation,
    MtzRejectionRFactor,
    MtzRejectionPartCorrel,
    MtzRejectionMinRefl,
    MtzRejectionOther,
} MtzRejectionReason;


class MtzMerger : public LoggableObject
{
private:
    std::mutex *reflCountMutex;
    std::vector<MtzPtr> allMtzs;
    std::vector<MtzPtr> someMtzs; // those we wish to merge
    MtzPtr mergedMtz;
    int cycle;
    int rejectedReflections;
    bool lowMemoryMode;
    double partCorrelThreshold;
    double correlationThreshold;
    double rFactorThreshold;
    int minReflectionCounts;
    ScalingType scalingType;
    bool excludeWorst;
    double rejectSigma;
    std::string filename;
    bool silent;
    int friedel;
    bool freeOnly;
    bool needToScale;
    
    void splitAllMtzs(std::vector<MtzPtr> &firstHalfMtzs, std::vector<MtzPtr> &secondHalfMtzs);
    MtzRejectionReason isMtzAccepted(MtzPtr mtz);
    std::map<MtzRejectionReason, int> rejectNums;
    std::mutex *rejectMutex;
    bool mtzIsPruned(MtzPtr mtz);
    void summary();
    void writeParameterCSV();
    void groupMillerThread(int offset);
    void groupMillers();
    void addMtzMillers(MtzPtr mtz);
    void makeEmptyReflectionShells(MtzPtr whichMtz);
    double maxResolution();
    static void groupMillerThreadWrapper(MtzMerger *object, int offset);
    std::string makeFilename(std::string prefix);
    
    void scaleIndividual(MtzPtr mtz);
    void scale();
    void fixSigmas();
    void removeReflections();
    void mergeMillersThread(int offset);
    void mergeMillers();
    int totalObservations();
    static void mergeMillersThreadWrapper(MtzMerger *object, int offset);
    static void writeAnomalousMtz(MtzPtr negative, MtzPtr positive, MtzPtr mean, std::string filename);
    void createAnomalousDiffMtz(MtzPtr negative, MtzPtr positive);

    void incrementRejectedReflections();
    
protected:

    
public:
    MtzMerger();
    void merge();
    void mergeFull(bool anomalous = false);
    void mergeAnomalous();
    
    void setCycle(int num);
    
    void setScalingType(ScalingType type)
    {
        scalingType = type;
    }
    
    void setAllMtzs(std::vector<MtzPtr> mtzs)
    {
        allMtzs = mtzs;
    }
    
    void setExcludeWorst(bool worst)
    {
        excludeWorst = worst;
    }

    void setFilename(std::string newName)
    {
        filename = newName;
    }
    
    std::string getFilename()
    {
        return filename;
    }
    
    MtzPtr getMergedMtz()
    {
        return mergedMtz;
    }
    
    void setSilent(bool newSilent)
    {
        silent = newSilent;
    }
    
    void setFriedel(int newFriedel)
    {
        friedel = newFriedel;
    }
    
    void setFreeOnly(bool free)
    {
        freeOnly = free;
    }
    
    void setNeedToScale(bool need)
    {
        needToScale = need;
    }
};

#endif /* defined(__cppxfel__MtzMerger__) */
