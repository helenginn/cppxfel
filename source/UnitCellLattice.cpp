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

std::vector<MatrixPtr> UnitCellLattice::symOperators;

void UnitCellLattice::setup(double a, double b, double c, double alpha, double beta, double gamma, int spaceGroupNum)
{
    spaceGroup = CSym::ccp4spg_load_by_ccp4_num(spaceGroupNum);
    
    Matrix::symmetryOperatorsForSpaceGroup(&symOperators, spaceGroup);
    
    unitCellOnly = Matrix::matrixFromUnitCell(a, b, c, alpha, beta, gamma);
    
    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    unitCellMatrix = MatrixPtr(new Matrix());
    unitCellMatrix->setComplexMatrix(unitCellOnly, rotationMat);
    
    unitCellMatrixInverse = unitCellMatrix->inverse3DMatrix();
    
    maxMillerIndexTrial = FileParser::getKey("MAX_MILLER_INDEX_TRIAL", 8);
    maxDistance = 0;
    
    for (int i = -maxMillerIndexTrial; i <= maxMillerIndexTrial; i++)
    {
        for (int j = -maxMillerIndexTrial; j <= maxMillerIndexTrial; j++)
        {
            for (int k = -maxMillerIndexTrial; k <= maxMillerIndexTrial; k++)
            {
                if (spaceGroupNum != 19 && spaceGroupNum != 178)
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
                
                spotVectors.push_back(newStandardVector);
            }
        }
    }
}