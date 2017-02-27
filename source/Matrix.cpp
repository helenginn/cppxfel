//
//  Matrix.cpp
//  RaddoseViewer
//
//  Created by Helen Ginn on 19/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#include "Matrix.h"
#include <cmath>
#include <cstdlib>
#include "csymlib.h"
#include <cstring>
#include "Vector.h"
#include "misc.h"
#include <sstream>
#include <iomanip>
#include "parameters.h"
#include <cctbx/sgtbx/space_group.h>
#include "Logger.h"
//#include <boost/python.hpp>
#include "FileParser.h"
#include "Hdf5ManagerProcessing.h"

using cctbx::sgtbx::rt_mx;
using cctbx::sgtbx::rot_mx;

MatrixPtr Matrix::identityMatrix = MatrixPtr(new Matrix());

double *Matrix::array()
{
    return components;
}

void Matrix::maxMillers(int (&millers)[3], double maxResolution)
{
    for (int i = 0; i < 3; i++)
    {
        vec testHKL = new_vector(i == 0, i == 1, i == 2);
        
        multiplyVector(&testHKL);
        
        double maxD = 1 / maxResolution;
        double hklLength = length_of_vector(testHKL);
        
        millers[i] = fabs(maxD / hklLength);
    }
}

void Matrix::setIdentity()
{
    for (int i = 0; i < 16; i++)
        components[i] = 0;
    
    for (int i = 0; i < 16; i += 5)
        components[i] = 1;
}

Matrix::Matrix(void)
{
    setIdentity();
    
    unitCell = MatrixPtr();
    rotation = MatrixPtr();
    
    eulerA = 0;
    eulerB = 0;
    eulerC = 0;
}

std::string Matrix::summary()
{
    vec bDummyVec = new_vector(1, 0, 0);
    
    this->multiplyVector(&bDummyVec);
    scale_vector_to_distance(&bDummyVec, 1);
    
    double alpha = acos(bDummyVec.l / length_of_vector(bDummyVec));
    double beta = atan(bDummyVec.k / bDummyVec.h);
    
    double theta, phi, psi;
    eulerAngles(&theta, &phi, &psi);
    
    std::ostringstream dummyVecStr;
    dummyVecStr << bDummyVec.h << "\t" << bDummyVec.k << "\t" << bDummyVec.l << "\t" << alpha << "\t" << beta << "\t" << theta << "\t" << phi << "\t" << psi;
    
    return dummyVecStr.str();
}

std::string Matrix::description(bool detailed, bool submatrix)
{
    std::ostringstream description;
    
    if (detailed == true && unitCell)
    {
        description << "unitcell ";
        description << unitCell->description(false, true) << std::endl;
        description << "rotation ";
        description << rotation->description(false, true);
        
        return description.str();
    }
    
    description << std::setprecision(14);
    
    if (!submatrix)
        description << "matrix ";
    
    description << components[0] << " ";
    description << components[4] << " ";
    description << components[8] << " ";
    description << components[1] << " ";
    description << components[5] << " ";
    description << components[9] << " ";
    description << components[2] << " ";
    description << components[6] << " ";
    description << components[10] << " ";
    
    return description.str();
}

void Matrix::sensibleComponents(double *returnedComponents[9])
{
    (*returnedComponents)[0] = components[0];
    (*returnedComponents)[1] = components[4];
    (*returnedComponents)[2] = components[8];
    (*returnedComponents)[3] = components[1];
    (*returnedComponents)[4] = components[5];
    (*returnedComponents)[5] = components[9];
    (*returnedComponents)[6] = components[2];
    (*returnedComponents)[7] = components[6];
    (*returnedComponents)[8] = components[10];
}

void Matrix::eulerAngles(double *theta, double *phi, double *psi, bool force)
{
    Matrix *chosenMat = rotation ? &*rotation : this;
    
    if (!(eulerA == 0 && eulerB == 0 && eulerC == 0) && !force)
    {
        *theta = eulerA;
        *phi = eulerB;
        *psi = eulerC;
    }
    
    double sinTheta = chosenMat->components[2];
    *theta = asin(sinTheta);
    double cosTheta = cos(*theta);
    
    *psi = atan2((chosenMat->components[6] / cosTheta), (chosenMat->components[10] / cosTheta));
    *phi = atan2((chosenMat->components[1] / cosTheta), (chosenMat->components[0] / cosTheta));
    
    eulerA = *theta;
    eulerB = *phi;
    eulerC = *psi;
}

double Matrix::similarityToRotationMatrix(MatrixPtr mat2, double tolerance, bool force)
{
    double theta1, phi1, psi1;
    eulerAngles(&theta1, &phi1, &psi1, force);
    
    double theta2, phi2, psi2;
    mat2->eulerAngles(&theta2, &phi2, &psi2, force);
    
    double thetaDiff = fabs(theta2 - theta1);
    
    if (thetaDiff > tolerance)
        return -1;
    
    double psiDiff = fabs(psi2 - psi1);
    
    if (psiDiff > tolerance)
        return -1;
    
    double phiDiff = fabs(phi2 - phi1);
    
    if (phiDiff > tolerance)
        return -1;
    
    double sumSqr = pow(theta1 - theta2, 2) + pow(psi1 - psi2, 2) + pow(phi1 - phi2, 2);
    
    return sqrt(sumSqr);
}

MatrixPtr Matrix::getNegativeCopy()
{
    MatrixPtr newMat = this->copy();
    
    for (int i = 0; i < 16; i++)
    {
        newMat->components[i] = - newMat->components[i];
    }
    
    return newMat;
}

// CCTBX_REWRITE: here

void Matrix::symmetryOperatorsForSpaceGroup(std::vector<MatrixPtr> *matrices, CSym::CCP4SPG *spaceGroup, double a, double b, double c, double alpha, double beta, double gamma, bool orthogonal)
{
    cctbx::sgtbx::space_group_symbols spaceGroupSymbol = cctbx::sgtbx::space_group_symbols(spaceGroup->spg_num);
    std::string hallSymbol = spaceGroupSymbol.hall();
    
    cctbx::sgtbx::space_group spaceGroupCctbx = cctbx::sgtbx::space_group(hallSymbol);
    
    scitbx::af::shared<rt_mx> allOps = spaceGroupCctbx.all_ops();
    scitbx::af::double6 params = scitbx::af::double6(a, b, c, alpha, beta, gamma);
    cctbx::uctbx::unit_cell uc = cctbx::uctbx::unit_cell(params);
    
    std::ostringstream logged;
    
    logged << "Number of symmetry operators for spg " << spaceGroup->spg_num << ": " << allOps.size() << std::endl;
    
    for (int j = 0; j < allOps.size(); j++)
    {
        rot_mx rotation = allOps[j].r();
        cctbx::uc_mat3 orthogonal_rotation = uc.matrix_cart(rotation);
        
        for (int i = 0; i < 1; i++)
        {
            MatrixPtr newMat = MatrixPtr(new Matrix());
            
            newMat->components[0] = orthogonal ? orthogonal_rotation[0] : rotation[0];
            newMat->components[4] = orthogonal ? orthogonal_rotation[1] : rotation[1];
            newMat->components[8] = orthogonal ? orthogonal_rotation[2] : rotation[2];
            
            newMat->components[1] = orthogonal ? orthogonal_rotation[3] : rotation[3];
            newMat->components[5] = orthogonal ? orthogonal_rotation[4] : rotation[4];
            newMat->components[9] = orthogonal ? orthogonal_rotation[5] : rotation[5];
            
            newMat->components[2] = orthogonal ? orthogonal_rotation[6] : rotation[6];
            newMat->components[6] = orthogonal ? orthogonal_rotation[7] : rotation[7];
            newMat->components[10] = orthogonal ? orthogonal_rotation[8] : rotation[8];
            
            MatrixPtr transposeMat = newMat->copy();
            matrices->push_back(transposeMat->transpose());
            MatrixPtr transNegMat = transposeMat->getNegativeCopy();
            matrices->push_back(transNegMat); // why
        }
    }
    
    Logger::mainLogger->addStream(&logged, LogLevelDetailed);
}

void Matrix::printDescription(bool detailed)
{
    Logger::mainLogger->addString(description(detailed));
}

void Matrix::prettyPrint()
{
    int width = 10;
    std::ostringstream logged;
    logged << "Pretty: " << std::endl;
    logged << std::setw(width) << std::setprecision(4) << components[0];
    logged << std::setw(width) << std::setprecision(4) << components[4];
    logged << std::setw(width) << std::setprecision(4) << components[8] << std::endl;
    logged << std::setw(width) << std::setprecision(4) << components[1];
    logged << std::setw(width) << std::setprecision(4) << components[5];
    logged << std::setw(width) << std::setprecision(4) << components[9] << std::endl;
    logged << std::setw(width) << std::setprecision(4) << components[2];
    logged << std::setw(width) << std::setprecision(4) << components[6];
    logged << std::setw(width) << std::setprecision(4) << components[10] << std::endl;
    Logger::mainLogger->addStream(&logged);
}

std::vector<double> Matrix::unitCellFromReciprocalUnitCell(double a, double b, double c, double alpha, double beta, double gamma)
{
    scitbx::af::double6 params = scitbx::af::double6(a, b, c, alpha, beta, gamma);
    cctbx::uctbx::unit_cell uc = cctbx::uctbx::unit_cell(params);
    
    cctbx::uctbx::unit_cell rc = uc.reciprocal();
    
    scitbx::af::double6 newParams = rc.parameters();
    
    std::vector<double> unitCell;

    for (int i = 0; i < newParams.size(); i++)
    {
        unitCell.push_back(newParams[i]);
    }
    
    return unitCell;
}

MatrixPtr Matrix::matrixFromUnitCell(double *unitCell)
{
    MatrixPtr aMatrix = MatrixPtr(new Matrix());
    aMatrix->rotation = MatrixPtr(new Matrix());
    aMatrix->unitCell = MatrixPtr(new Matrix());
    aMatrix->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
    
    return aMatrix->unitCell;
}

// might be totally crap
void Matrix::orientationMatrixUnitCell(double *a, double *b, double *c)
{
    MatrixPtr orientationTranspose = this->transpose();
    Matrix transposeTimesMatrix = *orientationTranspose * *this;
    MatrixPtr inverted = transposeTimesMatrix.inverse3DMatrix();
    
    double aSquared = (*inverted)[0];
    double bSquared = (*inverted)[5];
    double cSquared = (*inverted)[10];
    
    *a = sqrt(aSquared);
    *b = sqrt(bSquared);
    *c = sqrt(cSquared);
}

void Matrix::threeDimComponents(double **componentArray)
{
    (*componentArray)[0] = components[0];
    (*componentArray)[1] = components[4];
    (*componentArray)[2] = components[8];
    
    (*componentArray)[3] = components[1];
    (*componentArray)[4] = components[5];
    (*componentArray)[5] = components[9];
    
    (*componentArray)[6] = components[2];
    (*componentArray)[7] = components[6];
    (*componentArray)[8] = components[10];
}

void Matrix::add(MatrixPtr secondMatrix)
{
    for (int i = 0; i < 15; i++)
    {
        components[i] += secondMatrix->components[i];
    }
}

void Matrix::subtract(MatrixPtr secondMatrix)
{
    for (int i = 0; i < 15; i++)
    {
        components[i] -= secondMatrix->components[i];
    }
}


void Matrix::unitCellLengths(double **lengths)
{
    Matrix *mat = unitCell != NULL ? &*unitCell : this;
    
    double aLengthSqr = pow(mat->components[0], 2) + pow(mat->components[1], 2) + pow(mat->components[2], 2);
    double bLengthSqr = pow(mat->components[4], 2) + pow(mat->components[5], 2) + pow(mat->components[6], 2);
    double cLengthSqr = pow(mat->components[8], 2) + pow(mat->components[9], 2) + pow(mat->components[10], 2);
    
    (*lengths)[0] = 1 / sqrt(aLengthSqr);
    (*lengths)[1] = 1 / sqrt(bLengthSqr);
    (*lengths)[2] = 1 / sqrt(cLengthSqr);
}

/** Warning: only use for unit cells with 90 degree angles! */
void Matrix::scaleUnitCellAxes(double aScale, double bScale, double cScale)
{
    double *oldLengths = new double[3];
    unitCellLengths(&oldLengths);
    
    double *newLengths = new double[3];
    newLengths[0] = oldLengths[0] * aScale;
    newLengths[1] = oldLengths[1] * bScale;
    newLengths[2] = oldLengths[2] * cScale;
    
    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    
    changeOrientationMatrixDimensions(newLengths[0], newLengths[1], newLengths[2], unitCell[3], unitCell[4], unitCell[5]);
    
    delete [] oldLengths;
    delete [] newLengths;
}

void Matrix::recalculateOrientationMatrix()
{
    MatrixPtr newMat = unitCell->copy();
    
    *newMat *= *rotation;
    memcpy(components, newMat->components, 16 * sizeof(double));
}

void Matrix::setComplexMatrix(MatrixPtr newUnitCell, MatrixPtr newRotation)
{
    unitCell = newUnitCell;
    rotation = newRotation;
    
    recalculateOrientationMatrix();
}

void Matrix::changeOrientationMatrixDimensions(double newA, double newB, double newC, double alpha, double beta, double gamma)
{
    std::ostringstream logged;
    
    double *lengths = new double[3];
    unitCellLengths(&lengths);
    
    if (!unitCell || !rotation)
    {
        double aScale = lengths[0] / newA;
        double bScale = lengths[1] / newB;
        double cScale = lengths[2] / newC;
        
        components[0] *= aScale;
        components[1] *= aScale;
        components[2] *= aScale;

        components[4] *= bScale;
        components[5] *= bScale;
        components[6] *= bScale;

        components[8] *= cScale;
        components[9] *= cScale;
        components[10] *= cScale;

        delete [] lengths;
        
        return;
    }
    
    logged << "Original cell axes: " << lengths[0] << ", " << lengths[1] << ", " << lengths[2];
    
    double cosAlpha = cos(alpha * M_PI / 180);
    double cosBeta = cos(beta * M_PI / 180);
    double cosGamma = cos(gamma * M_PI / 180);
    double sinBeta = sin(beta * M_PI / 180);
    double sinGamma = sin(gamma * M_PI / 180);
    
    double denom = sinBeta * sinGamma;
    double rCosAlpha = cosBeta * cosGamma - cosAlpha;
    rCosAlpha /= denom;
    
    double s1rca2 = sqrt(1. - rCosAlpha * rCosAlpha);
    
    MatrixPtr ortho = MatrixPtr(new Matrix());
    ortho->components[0] = newA;
    ortho->components[4] = cosGamma * newB;
    ortho->components[8] = cosBeta * newC;
    ortho->components[1] = 0.;
    ortho->components[5] = sinGamma * newB;
    ortho->components[9] = -sinBeta * rCosAlpha * newC;
    ortho->components[2] = 0.;
    ortho->components[6] = 0.;
    ortho->components[10] = sinBeta * newC * s1rca2;
    
    this->unitCell = ortho->inverse3DMatrix();
    
    recalculateOrientationMatrix();
    
    unitCellLengths(&lengths);
    logged << "; new cell axes: " << lengths[0] << ", " << lengths[1] << ", " << lengths[2] << std::endl;
    
    delete [] lengths;
    
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
}

// Order = b * this
Matrix Matrix::operator*=(Matrix &b)
{
    Matrix newMat;
    
    (newMat)[0] = b[0] * components[0] + b[4] * components[1]
    + b[8] * components[2];//; + b[12] * components[3];
    (newMat)[1] = b[1] * components[0] + b[5] * components[1]
    + b[9] * components[2];// + b[13] * components[3];
    (newMat)[2] = b[2] * components[0] + b[6] * components[1]
    + b[10] * components[2];// + b[14] * components[3];
    /*	(newMat)[3] = b[3] * components[0] + b[7] * components[1]
     + b[11] * components[2] + b[15] * components[3];
     */
    
    (newMat)[4] = b[0] * components[4] + b[4] * components[5]
    + b[8] * components[6];// + b[12] * components[7];
    (newMat)[5] = b[1] * components[4] + b[5] * components[5]
    + b[9] * components[6];// + b[13] * components[7];
    (newMat)[6] = b[2] * components[4] + b[6] * components[5]
    + b[10] * components[6];// + b[14] * components[7];
    /*	(newMat)[7] = b[3] * components[4] + b[7] * components[5]
     + b[11] * components[6] + b[15] * components[7];
     */
    
    (newMat)[8] = b[0] * components[8] + b[4] * components[9]
    + b[8] * components[10];// + b[12] * components[11];
    (newMat)[9] = b[1] * components[8] + b[5] * components[9]
    + b[9] * components[10];// + b[13] * components[11];
    (newMat)[10] = b[2] * components[8] + b[6] * components[9]
    + b[10] * components[10] ;//+ b[14] * components[11];
    /*	(newMat)[11] = b[3] * components[8] + b[7] * components[9]
     + b[11] * components[10] + b[15] * components[11];
     */
    
    (newMat)[12] = b[0] * components[12] + b[4] * components[13]
    + b[8] * components[14];// + b[12] * components[15];
    (newMat)[13] = b[1] * components[12] + b[5] * components[13]
    + b[9] * components[14];// + b[13] * components[15];
    (newMat)[14] = b[2] * components[12] + b[6] * components[13]
    + b[10] * components[14];// + b[14] * components[15];
    /*	(newMat)[15] = b[3] * components[12] + b[7] * components[13]
     + b[11] * components[14] + b[15] * components[15];
     */
    memcpy(this->components, newMat.components, 16 * sizeof(double));
    
    return *this;
}

bool Matrix::isIdentity()
{
    for (int i = 0; i < 16; i++)
    {
        if (i == 0 || i == 5 || i == 10 || i == 15)
        {
            if (fabs(components[i] - 1) > 0.00001)
                return false;
        }
        else if (components[i] > 0.00001)
        {
            return false;
        }
    }
    
    return true;
}

void Matrix::multiply(double scalar)
{
    for (int i = 0; i < 16; i++)
    {
        if (i == 3 || i == 7 || i > 11)
            continue;
        
        components[i] *= scalar;
    }
}

void Matrix::multiply(Matrix &b)
{
    (*this) *= b;
}

// Order = this * b (I think)
void Matrix::preMultiply(Matrix &b)
{
    Matrix newMat = b * (*this);
    memcpy(this->components, newMat.components, sizeof(double) * 16);
}

Matrix Matrix::operator*(Matrix &b)
{
    Matrix newMat;
    memcpy(newMat.components, this->components, 16 * sizeof(double));
    
    newMat *= b;
    
    return newMat;
}

MatrixPtr Matrix::copy(void)
{
    MatrixPtr newMat = MatrixPtr(new Matrix());
    memcpy(newMat->components, this->components, 16 * sizeof(double));
    
    if (isComplex())
    {
        newMat->unitCell = unitCell->copy();
        newMat->rotation = rotation->copy();
    }
    
    newMat->eulerA = eulerA;
    newMat->eulerB = eulerB;
    newMat->eulerC = eulerC;
    
    return newMat;
}

Matrix::Matrix(double *newComponents)
{
    components[0] = newComponents[0];
    components[1] = newComponents[3];
    components[2] = newComponents[6];
    components[3] = 0;
    
    components[4] = newComponents[1];
    components[5] = newComponents[4];
    components[6] = newComponents[7];
    components[7] = 0;
    
    components[8] = newComponents[2];
    components[9] = newComponents[5];
    components[10] = newComponents[8];
    components[11] = 0;
    
    components[12] = 0;
    components[13] = 0;
    components[14] = 0;
    components[15] = 1;
    
    eulerA = 0;
    eulerB = 0;
    eulerC = 0;
}

void Matrix::translate(double x, double y, double z)
{
    components[12] += x;
    components[13] += y;
    components[14] += z;
    
}

void Matrix::scale(double a)
{
    scale(a, a, a);
}

void Matrix::scale(double a, double b, double c)
{
    MatrixPtr scaleMat = MatrixPtr(new Matrix());
    scaleMat->components[0] *= a;
    scaleMat->components[5] *= b;
    scaleMat->components[10] *= c;
    
    this->multiply(*scaleMat);
}

void Matrix::rotateHK(double hRot, double kRot)
{
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    
    this->rotate(hRad, kRad, 0);
}

void Matrix::rotateABC(MatrixPtr oldMatrix, double aRot, double bRot, double cRot)
{
    double aRad = aRot * M_PI / 180;
    double bRad = bRot * M_PI / 180;
    double cRad = cRot * M_PI / 180;
    
    vec vecs[3];
    
    for (int i = 0; i < 3; i++)
    {
        vecs[i] = new_vector(i == 0, i == 1, i == 2);
        oldMatrix->multiplyVector(&vecs[i]);
    }
    
    this->rotateRoundUnitVector(vecs[0], aRad);
    this->rotateRoundUnitVector(vecs[1], bRad);
    this->rotateRoundUnitVector(vecs[2], cRad);
}

void Matrix::rotate(double alpha, double beta, double gamma)
{
    if (alpha != 0)
    {
        double xAxis[] = { 1, 0, 0 };
        rotateRoundUnitVector(xAxis, alpha);
    }
    
    if (beta != 0)
    {
        double yAxis[] = { 0, 1, 0 };
        rotateRoundUnitVector(yAxis, beta);
    }
    
    if (gamma != 0)
    {
        double zAxis[] = { 0, 0, 1 };
        rotateRoundUnitVector(zAxis, gamma);
    }
}

void Matrix::rotateModelAxes(double alpha, double beta, double gamma)
{
    Matrix *alphaMat = new Matrix();
    (*alphaMat)[5] = cos(-alpha);
    (*alphaMat)[9] = -sin(-alpha);
    (*alphaMat)[6] = sin(-alpha);
    (*alphaMat)[10] = cos(-alpha);
    
    Matrix *betaMat = new Matrix();
    (*betaMat)[0] = cos(-beta);
    (*betaMat)[8] = sin(-beta);
    (*betaMat)[2] = -sin(-beta);
    (*betaMat)[10] = cos(-beta);
    
    Matrix *gammaMat = new Matrix();
    (*gammaMat)[0] = cos(-gamma);
    (*gammaMat)[4] = -sin(-gamma);
    (*gammaMat)[1] = sin(-gamma);
    (*gammaMat)[5] = cos(-gamma);
    
    Matrix *chosenMatrix = this;
    
    Matrix alphaThis = *alphaMat * *chosenMatrix;
    Matrix betaAlphaThis = *betaMat * alphaThis;
    Matrix allThis = *gammaMat * betaAlphaThis;
    
    memcpy(chosenMatrix->components, allThis.components, 16 * sizeof(double));
    
}

void Matrix::rotateRoundUnitVector(vec unitVector, double radians)
{
    double *axis = new double[3];
    
    axis[0] = unitVector.h;
    axis[1] = unitVector.k;
    axis[2] = unitVector.l;
    
    rotateRoundUnitVector(axis, radians);
    
    delete [] axis;
}

void Matrix::rotateRoundUnitVector(double *unitVector, double radians)
{
    Matrix matrix = Matrix();
    
    (matrix)[0] = cos(radians) + pow(unitVector[0], 2) * (1 - cos(radians));
    (matrix)[4] = unitVector[0] * unitVector[1] * (1 - cos(radians))
    - unitVector[2] * sin(radians);
    (matrix)[8] = unitVector[0] * unitVector[2] * (1 - cos(radians))
    + unitVector[1] * sin(radians);
    
    (matrix)[1] = unitVector[1] * unitVector[0] * (1 - cos(radians))
    + unitVector[2] * sin(radians);
    (matrix)[5] = cos(radians) + pow(unitVector[1], 2) * (1 - cos(radians));
    (matrix)[9] = unitVector[1] * unitVector[2] * (1 - cos(radians))
    - unitVector[0] * sin(radians);
    
    (matrix)[2] = unitVector[2] * unitVector[0] * (1 - cos(radians))
    - unitVector[1] * sin(radians);
    (matrix)[6] = unitVector[2] * unitVector[1] * (1 - cos(radians))
    + unitVector[0] * sin(radians);
    (matrix)[10] = cos(radians) + pow(unitVector[2], 2) * (1 - cos(radians));
    
    Matrix *chosenMatrix = this;
    
    if (rotation)
    {
        chosenMatrix = &*rotation;
    }
    
    chosenMatrix->multiply(matrix);
    
    if (rotation)
    {
        recalculateOrientationMatrix();
    }
}

void Matrix::multiplyVector(vec *vector)
{
    vec oldVec = copy_vector(*vector);
    
    vector->h = components[0] * oldVec.h + components[4] * oldVec.k
    + components[8] * oldVec.l;
    vector->k = components[1] * oldVec.h + components[5] * oldVec.k
    + components[9] * oldVec.l;
    vector->l = components[2] * oldVec.h + components[6] * oldVec.k
    + components[10] * oldVec.l;
}

void Matrix::newMultiplyVector(double *vector[])
{
    double *oldVector = (double *) malloc(sizeof(double) * 4);
    memcpy(oldVector, *vector, sizeof(double) * 4);
    
    (*vector)[0] = components[0] * (*vector)[0] + components[4] * (*vector)[1]
    + components[8] * (*vector)[2] + components[12] * (*vector)[3];
    (*vector)[1] = components[1] * (*vector)[0] + components[5] * (*vector)[1]
    + components[9] * (*vector)[2] + components[13] * (*vector)[3];
    (*vector)[2] = components[2] * (*vector)[0] + components[6] * (*vector)[1]
    + components[10] * (*vector)[2] + components[14] * (*vector)[3];
    (*vector)[3] = components[3] * (*vector)[0] + components[7] * (*vector)[1]
    + components[11] * (*vector)[2] + components[15] * (*vector)[3];
    
    free(oldVector);
    
}

void Matrix::identity(void)
{
    Matrix *newMatrix = new Matrix();
    
    (*this) = (*newMatrix);
}

double Matrix::getEwaldSphere(vec *vector)
{
    vec index = new_vector(vector->h, vector->k, vector->l);
    this->multiplyVector(&index);
    
    double ewald_radius = index.h * index.h + index.k * index.k
    + index.l * index.l;
    if (index.l == 0)
        return 0;
    
    ewald_radius /= (0 - 2 * index.l);
    double ewald_wavelength = 1 / ewald_radius;
    
    return ewald_wavelength;
}

MatrixPtr Matrix::matrixFromEulerAngles(double theta, double phi, double psi)
{
    MatrixPtr matrix = MatrixPtr(new Matrix());
    matrix->rotate(psi, theta, phi);
    
    return matrix;
}

Matrix Matrix::inverse2DMatrix()
{
    Matrix newMat = Matrix();
    
    double a = components[0];
    double b = components[4];
    double c = components[1];
    double d = components[5];
    
    newMat.components[0] = d;
    newMat.components[5] = a;
    newMat.components[1] = -b;
    newMat.components[4] = -c;
    
    double scale = 1 / (a * d - b * c);
    
    for (int i = 0; i < 16; i++)
    {
        newMat[i] *= scale;
    }
    
    return newMat;
}

void Matrix::rotate2D(double angle)
{
    Matrix matrix = Matrix();
    
    matrix[0] = cos(angle);
    matrix[1] = sin(angle);
    matrix[4] = -sin(angle);
    matrix[5] = cos(angle);
    
    this->preMultiply(matrix);
}

double invertValue(double topLeft, double bottomRight, double topRight, double bottomLeft)
{
    return topLeft * bottomRight - bottomLeft * topRight;
}

// CCTBX_REWRITE: here
MatrixPtr Matrix::inverse3DMatrix()
{
    double a = components[0];
    double b = components[4];
    double c = components[8];
    double d = components[1];
    double e = components[5];
    double f = components[9];
    double g = components[2];
    double h = components[6];
    double i = components[10];
    
    double det = determinant();
    MatrixPtr newMat = MatrixPtr(new Matrix());
    
    newMat->components[0] = (e * i - f * h) / det;
    newMat->components[4] = -(d * i - f * g) / det;
    newMat->components[8] = (d * h - e * g) / det;
    newMat->components[1] = - (b * i - c * h) / det;
    newMat->components[5] = (a * i - c * g) / det;
    newMat->components[9] = - (a * h - b * g) / det;
    newMat->components[2] = (b * f - c * e) / det;
    newMat->components[6] = - (a * f - c * d) / det;
    newMat->components[10] = (a * e - b * d) / det;

    return newMat;
}

double Matrix::determinant()
{
    double a = components[0];
    double b = components[4];
    double c = components[8];
    double d = components[1];
    double e = components[5];
    double f = components[9];
    double g = components[2];
    double h = components[6];
    double i = components[10];
    
    double det = a*e*i + b*f*g + c*d*h + c*e*g - b*d*i - a*f*h;

    return det;
}

double Matrix::trace()
{
    double trace = 0;
    trace += components[0];
    trace += components[5];
    trace += components[10];
    
    return trace;
}

MatrixPtr Matrix::transpose()
{
    MatrixPtr transpose = MatrixPtr(new Matrix());
    
    (*transpose)[0] = components[0];
    (*transpose)[1] = components[4];
    (*transpose)[2] = components[8];
    (*transpose)[3] = components[12];
    (*transpose)[4] = components[1];
    (*transpose)[5] = components[5];
    (*transpose)[6] = components[9];
    (*transpose)[8] = components[2];
    (*transpose)[9] = components[6];
    (*transpose)[10] = components[10];
    (*transpose)[11] = components[14];
    (*transpose)[12] = components[3];
    (*transpose)[13] = components[7];
    (*transpose)[14] = components[11];
    (*transpose)[15] = components[15];
    
    return transpose;
}

void Matrix::print(void)
{
    std::cout << components[0] << "\t" << components[4] << "\t" << components[8]
    << std::endl;
    std::cout << components[1] << "\t" << components[5] << "\t" << components[9]
    << std::endl;
    std::cout << components[2] << "\t" << components[6] << "\t"
    << components[10] << std::endl;
}

void Matrix::translation(double **vector)
{
    memcpy(*vector, &components[12], sizeof(double) * 3);
    
}

// CCTBX_REWRITE: here
cctbx::miller::index<double> Matrix::multiplyIndex(cctbx::miller::index<> *index)
{
    int h = index->as_tiny()[0];
    int k = index->as_tiny()[1];
    int l = index->as_tiny()[2];
    
    vec hkl = new_vector(h, k, l);
    
    this->multiplyVector(&hkl);
    
    cctbx::miller::index<double> newIndex = cctbx::miller::index<double>(hkl.h, hkl.k, hkl.l);
    
    return newIndex;
}

bool Matrix::writeToHdf5(std::string address)
{
    Hdf5ManagerProcessingPtr processingManager = Hdf5ManagerProcessing::getProcessingManager();
    
    hsize_t dim = 16;
    
    bool success = processingManager->createDataset(address, 1, &dim, H5T_NATIVE_DOUBLE);
    
    if (!success)
    {
        return success;
    }
    
    return processingManager->writeDataset(address, (void **)&components, H5T_NATIVE_DOUBLE);    
}

void Matrix::copyComponents(MatrixPtr mat2)
{
    memcpy(this->components, mat2->components, sizeof(double) * 16);
}
