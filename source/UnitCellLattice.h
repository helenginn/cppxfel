//
//  UnitCellLattice.h
//  cppxfel
//
//  Created by Helen Ginn on 10/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__UnitCellLattice__
#define __cppxfel__UnitCellLattice__

#include <stdio.h>
#include "FreeLattice.h"
#include "csymlib.h"
#include "Matrix.h"
#include "parameters.h"

class UnitCellLattice : public FreeLattice
{
private:
    static UnitCellLatticePtr mainLattice;
    CSym::CCP4SPG *spaceGroup;
    MatrixPtr unitCellMatrix;
    MatrixPtr unitCellOnly;
    MatrixPtr unitCellMatrixInverse;
    static std::vector<MatrixPtr> symOperators;
    int maxMillerIndexTrial;
    double maxDistance;
    double minDistance;
    double powderStep;
    double unitCell[6];
    PowderHistogram histogram;
    std::vector<SpotVectorPtr> uniqueSymVectors;
    CSVPtr weightedUnitCell;
    CSVPtr weightedAngles;
    void updateUnitCellData();
    void lockUnitCellDimensions(double *a, double *b, double *c, double *alpha, double *beta, double *gamma);
    double distanceToAngleRatio;
    std::mutex setupLock;
    static bool setupLattice;
    
    double _aDim;
    double _bDim;
    double _cDim;
    double _alpha;
    double _beta;
    double _gamma;
    
    UnitCellLattice(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum, double resolution = 0) : FreeLattice()
    {
        setup(a, b, c, alpha, beta, gamma, spaceGroupNum, resolution);
    }

public:
    void setup(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum, double resolution = 0);
    void refineUnitCell(PowderHistogram histogram);
    void getMaxMillerIndicesForResolution(double resolution, int *hMax, int *kMax, int *lMax);
    void weightUnitCell();
    
    static UnitCellLatticePtr getMainLattice();
    
    int standardVectorCount()
    {
        return (int)spotVectors.size();
    }
    
    int uniqueSymVectorCount()
    {
        return (int)uniqueSymVectors.size();
    }
    
    SpotVectorPtr uniqueSymVector(int i)
    {
        return uniqueSymVectors[i];
    }
    
    SpotVectorPtr standardVector(int i)
    {
        return spotVectors[i];
    }
    
    int symOperatorCount()
    {
        return (int)symOperators.size();
    }
    
    MatrixPtr symOperator(int i)
    {
        return symOperators[i];
    }
    
    double getMaxDistance()
    {
        return maxDistance;
    }
    
    MatrixPtr getUnitCellOnly()
    {
        return unitCellOnly;
    }
    
    double getMinDistance()
    {
        return minDistance;
    }
    
    double weightForDistance(double distance);
    double weightForAngle(double angle);
    double weightForCosine(double cosine);
    
    std::vector<SpotVectorPtr> getStandardVectors()
    {
        return spotVectors;
    }
    
    
    CSVPtr getWeightedUnitCell()
    {
        return weightedUnitCell;
    }
    
    static double getUnitCellA(void *object)
    {
        return static_cast<UnitCellLattice *>(object)->_aDim;
    }
    
    static void setUnitCellA(void *object, double newDim)
    {
        static_cast<UnitCellLattice *>(object)->_aDim = newDim;
        static_cast<UnitCellLattice *>(object)->updateUnitCellData();
    }
    
    static double getUnitCellB(void *object)
    {
        return static_cast<UnitCellLattice *>(object)->_bDim;
    }
    
    static void setUnitCellB(void *object, double newDim)
    {
        static_cast<UnitCellLattice *>(object)->_bDim = newDim;
        static_cast<UnitCellLattice *>(object)->updateUnitCellData();
    }
    
    static double getUnitCellC(void *object)
    {
        return static_cast<UnitCellLattice *>(object)->_cDim;
    }
    
    static void setUnitCellC(void *object, double newDim)
    {
        static_cast<UnitCellLattice *>(object)->_cDim = newDim;
        static_cast<UnitCellLattice *>(object)->updateUnitCellData();
    }
    
    static double getUnitCellAlpha(void *object)
    {
        return static_cast<UnitCellLattice *>(object)->_alpha;
    }
    
    static void setUnitCellAlpha(void *object, double newDim)
    {
        static_cast<UnitCellLattice *>(object)->_alpha = newDim;
        static_cast<UnitCellLattice *>(object)->updateUnitCellData();
    }
    
    static double getUnitCellBeta(void *object)
    {
        return static_cast<UnitCellLattice *>(object)->_beta;
    }
    
    static void setUnitCellBeta(void *object, double newDim)
    {
        static_cast<UnitCellLattice *>(object)->_beta = newDim;
        static_cast<UnitCellLattice *>(object)->updateUnitCellData();
    }
    
    static double getUnitCellGamma(void *object)
    {
        return static_cast<UnitCellLattice *>(object)->_gamma;
    }
    
    static void setUnitCellGamma(void *object, double newDim)
    {
        static_cast<UnitCellLattice *>(object)->_gamma = newDim;
        static_cast<UnitCellLattice *>(object)->updateUnitCellData();
    }
};

#endif /* defined(__cppxfel__UnitCellLattice__) */
