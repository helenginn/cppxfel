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
#include "Logger.h"
#include "UnitCellLattice.h"
#include "FileParser.h"
#include "Hdf5ManagerProcessing.h"

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

MatrixPtr Matrix::matFromCCP4(CSym::ccp4_symop *symop)
{
    MatrixPtr mat = MatrixPtr(new Matrix());
    mat->components[0] = symop->rot[0][0];
    mat->components[4] = symop->rot[0][1];
    mat->components[8] = symop->rot[0][2];
    
    mat->components[1] = symop->rot[1][0];
    mat->components[5] = symop->rot[1][1];
    mat->components[9] = symop->rot[1][2];
    
    mat->components[2] = symop->rot[2][0];
    mat->components[6] = symop->rot[2][1];
    mat->components[10] = symop->rot[2][2];
    
    return mat;
}

void Matrix::symmetryOperatorsForSpaceGroup(std::vector<MatrixPtr> *matrices, CSym::CCP4SPG *spaceGroup, std::vector<double> cell)
{
    MatrixPtr unitCell = Matrix::matrixFromUnitCell(cell);
    MatrixPtr reverse = unitCell->inverse3DMatrix();

	if (!spaceGroup)
	{
		std::ostringstream logged;
		logged << "Need SPACE_GROUP to load symmetry operators." << std::endl;
		LoggableObject::staticLogAndExit(logged);
	}

    for (int j = 0; j < spaceGroup->nsymop; j++)
    {
        MatrixPtr newMat = matFromCCP4(&spaceGroup->symop[j]);
        MatrixPtr ortho = newMat;
        ortho->preMultiply(*unitCell);
        ortho->multiply(*reverse);
        
        matrices->push_back(ortho);
    }
}

void Matrix::printDescription(bool detailed)
{
    Logger::mainLogger->addString(description(detailed));
}

MatrixPtr Matrix::matrixFromUnitCell(std::vector<double> &unitCell)
{
    MatrixPtr aMatrix = MatrixPtr(new Matrix());
    aMatrix->rotation = MatrixPtr(new Matrix());
    aMatrix->unitCell = MatrixPtr(new Matrix());
    aMatrix->changeOrientationMatrixDimensions(unitCell);
    
    return aMatrix->unitCell;
}

void Matrix::subtract(MatrixPtr secondMatrix)
{
    for (int i = 0; i < 15; i++)
    {
        components[i] -= secondMatrix->components[i];
    }
}


void Matrix::unitCellLengths(double *lengths)
{
    Matrix *mat = unitCell != NULL ? &*unitCell : this;
    
    double aLengthSqr = pow(mat->components[0], 2) + pow(mat->components[1], 2) + pow(mat->components[2], 2);
    double bLengthSqr = pow(mat->components[4], 2) + pow(mat->components[5], 2) + pow(mat->components[6], 2);
    double cLengthSqr = pow(mat->components[8], 2) + pow(mat->components[9], 2) + pow(mat->components[10], 2);
    
    (lengths)[0] = 1 / sqrt(aLengthSqr);
    (lengths)[1] = 1 / sqrt(bLengthSqr);
    (lengths)[2] = 1 / sqrt(cLengthSqr);
}

void Matrix::recalculateOrientationMatrix()
{
    copyComponents(unitCell);
    multiply(*rotation);
}

void Matrix::setComplexMatrix(MatrixPtr newUnitCell, MatrixPtr newRotation)
{
    unitCell = newUnitCell;
    rotation = newRotation;
    
    recalculateOrientationMatrix();
}

void Matrix::changeOrientationMatrixDimensions(std::vector<double> cell)
{
    if (!unitCell || !rotation)
    {
        return;
    }
    
    double cosAlpha = cos(cell[3] * M_PI / 180);
    double cosBeta = cos(cell[4] * M_PI / 180);
    double cosGamma = cos(cell[5] * M_PI / 180);
    double sinBeta = sin(cell[4] * M_PI / 180);
    double sinGamma = sin(cell[5] * M_PI / 180);
    
    double denom = sinBeta * sinGamma;
    double rCosAlpha = cosBeta * cosGamma - cosAlpha;
    rCosAlpha /= denom;
    
    double s1rca2 = sqrt(1. - rCosAlpha * rCosAlpha);
    
    MatrixPtr ortho = MatrixPtr(new Matrix());
    ortho->components[0] = cell[0];
    ortho->components[1] = cosGamma * cell[1];
    ortho->components[2] = cosBeta * cell[2];
    ortho->components[4] = 0.;
    ortho->components[5] = sinGamma * cell[1];
    ortho->components[6] = -sinBeta * rCosAlpha * cell[2];
    ortho->components[8] = 0.;
    ortho->components[9] = 0.;
    ortho->components[10] = sinBeta * cell[2] * s1rca2;
    
    this->unitCell = ortho->inverse3DMatrix();

    recalculateOrientationMatrix();
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

void Matrix::rotate(double alpha, double beta, double gamma)
{
    if (alpha != 0)
    {
        vec xAxis = new_vector(1, 0, 0);
        rotateRoundUnitVector(xAxis, alpha);
    }
    
    if (beta != 0)
    {
        vec yAxis = new_vector(0, 1, 0);
        rotateRoundUnitVector(yAxis, beta);
    }
    
    if (gamma != 0)
    {
        vec zAxis = new_vector(0, 0, 1);
        rotateRoundUnitVector(zAxis, gamma);
    }
}


void Matrix::rotateRoundUnitVector(vec axis, double radians)
{
    Matrix matrix = Matrix();
    
    double x = axis.h;
    double x2 = axis.h * axis.h;

    double y = axis.k;
    double y2 = axis.k * axis.k;
    
    double z = axis.l;
    double z2 = axis.l * axis.l;
    
    double cosa = cos(radians);
    double sina = sin(radians);
    
    (matrix)[0] = cosa + x2 * (1 - cosa);
    (matrix)[4] = x * y * (1 - cosa) - z * sina;
    (matrix)[8] = x * z * (1 - cosa) + y * sina;
    
    (matrix)[1] = y * x * (1 - cosa) + z * sina;
    (matrix)[5] = cosa + y2 * (1 - cosa);
    (matrix)[9] = z * y * (1 - cosa) - x * sina;
    
    (matrix)[2] = z * x * (1 - cosa) - y * sina;
    (matrix)[6] = z * y * (1 - cosa) + x * sina;
    (matrix)[10] = cosa + z2 * (1 - cosa);
    
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

void Matrix::identity(void)
{
    Matrix *newMatrix = new Matrix();
    
    (*this) = (*newMatrix);
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

    return newMat->transpose();
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
    
    double det = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;

    
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

MatrixPtr Matrix::randomOrientation()
{
    MatrixPtr rotation = MatrixPtr(new Matrix());
    
    for (int i = 0; i < 8; i++)
    {
        double alpha = (double)rand() / RAND_MAX * 2 * M_PI;
        double beta = (double)rand() / RAND_MAX * 2 * M_PI;
        double gamma = (double)rand() / RAND_MAX * 2 * M_PI;
        
        MatrixPtr single = MatrixPtr(new Matrix());
        single->rotate(alpha, beta, gamma);
        rotation->multiply(*single);
    }
    
    return rotation;
}

MatrixPtr Matrix::randomOrientationMatrix()
{
    MatrixPtr rotation = randomOrientation();
    MatrixPtr unitCell = UnitCellLattice::getMainLattice()->getUnitCellOnly()->copy();
    
    MatrixPtr mat = MatrixPtr(new Matrix());
    mat->setComplexMatrix(unitCell, rotation);
    
    return mat;
}