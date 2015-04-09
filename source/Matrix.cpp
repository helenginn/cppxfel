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
#include <cstring>
#include "Vector.h"
#include "misc.h"
#include <sstream>
#include <iomanip>
#include "parameters.h"
#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <scitbx/vec3.h>
#include <cctbx/crystal_orientation.h>
#include "Logger.h"

/*+(void) populate: (GLdouble*) aGLMatrix
 fromFrustumLeft: (GLdouble) left
 andRight: (GLdouble) right
 andBottom: (GLdouble) bottom
 andTop: (GLdouble) top
 andNear: (GLdouble) near
 andFar: (GLdouble) far */

double *Matrix::array()
{
	return components;
}

Matrix::Matrix(void)
{
	for (int i = 0; i < 16; i++)
		components[i] = 0;

	for (int i = 0; i < 16; i += 5)
		components[i] = 1;
}

std::string Matrix::description()
{
	std::ostringstream description;

	description << std::setprecision(14);

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

void Matrix::printDescription()
{
	std::cout << description() << std::endl;
/*
	for (int i = 0; i < 9; i += 3)
	{
		for (int j = i; j < i + 3; j++)
		{
			std::cout << components[j] << "\t";
		}
		std::cout << std::endl;
	}*/
}

MatrixPtr Matrix::matrixFromUnitCell(double a, double b, double c, double alpha, double beta, double gamma)
{
    scitbx::af::double6 params = scitbx::af::double6(a, b, c, alpha, beta, gamma);
    cctbx::uctbx::unit_cell uc = cctbx::uctbx::unit_cell(params);
    
    scitbx::mat3<double> mat = uc.orthogonalization_matrix().inverse();
    
    MatrixPtr aMatrix = MatrixPtr(new Matrix());
    aMatrix->assignFromCctbxMatrix(mat);
    
    return aMatrix;
}

void Matrix::orientationMatrixUnitCell(double *a, double *b, double *c)
{
    Matrix orientationTranspose = this->transpose();
    Matrix transposeTimesMatrix = orientationTranspose * *this;
    Matrix inverted = transposeTimesMatrix.inverse3DMatrix();
    
    double aSquared = inverted[0];
    double bSquared = inverted[5];
    double cSquared = inverted[10];
    
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

scitbx::mat3<double> Matrix::cctbxMatrix()
{
    scitbx::mat3<double> cctbxMat = scitbx::mat3<double>();
    
    cctbxMat(0, 0) = components[0];
    cctbxMat(1, 0) = components[4];
    cctbxMat(2, 0) = components[8];
    
    cctbxMat(0, 1) = components[1];
    cctbxMat(1, 1) = components[5];
    cctbxMat(2, 1) = components[9];
    
    cctbxMat(0, 2) = components[2];
    cctbxMat(1, 2) = components[6];
    cctbxMat(2, 2) = components[10];
    
    return cctbxMat;
}

void Matrix::unitCellLengths(double **lengths)
{
    scitbx::mat3<double> cctbxMat = cctbxMatrix();
    scitbx::mat3<double> inverseMat = cctbxMat.error_minimizing_inverse(10);
    
    for (int i = 0; i < 3; i++)
    {
        scitbx::vec3<double> axis = scitbx::vec3<double>(i == 0, i == 1, i == 2);
        scitbx::vec3<double> bigAxis = inverseMat * axis;
        (*lengths)[i] = bigAxis.length();
    }
    

}

void Matrix::assignFromCctbxMatrix(scitbx::mat3<double> newMat)
{
    components[0] = newMat(0, 0);
    components[4] = newMat(1, 0);
    components[8] = newMat(2, 0);

    components[1] = newMat(0, 1);
    components[5] = newMat(1, 1);
    components[9] = newMat(2, 1);

    components[2] = newMat(0, 2);
    components[6] = newMat(1, 2);
    components[10] = newMat(2, 2);
}

void Matrix::changeOrientationMatrixDimensions(double newA, double newB, double newC, double alpha, double beta, double gamma)
{
    ostringstream logged;
    
    scitbx::mat3<double> cctbxMat = cctbxMatrix();
    
    double *lengths = new double[3];
    unitCellLengths(&lengths);
    
    logged << "Original cell axes: " << lengths[0] << ", " << lengths[1] << ", " << lengths[2];
    
    scitbx::af::double6 params = scitbx::af::double6(lengths[0], lengths[1], lengths[2], alpha, beta, gamma);
    cctbx::uctbx::unit_cell uc = cctbx::uctbx::unit_cell(params);
    
    scitbx::mat3<double> ortho = uc.orthogonalization_matrix();
    
    scitbx::mat3<double> rotation = cctbxMat * ortho;
    
    scitbx::af::double6 newParams = scitbx::af::double6(newA, newB, newC, alpha, beta, gamma);
    cctbx::uctbx::unit_cell newUnitCell = cctbx::uctbx::unit_cell(newParams);
    scitbx::mat3<double> newOrtho = newUnitCell.orthogonalization_matrix().error_minimizing_inverse(10);
    scitbx::mat3<double> newMat = rotation * newOrtho;
    
    assignFromCctbxMatrix(newMat);
    unitCellLengths(&lengths);
    logged << "; new cell axes: " << lengths[0] << ", " << lengths[1] << ", " << lengths[2] << std::endl;
    
    delete [] lengths;
    
    Logger::mainLogger->addStream(&logged);
}


Matrix Matrix::operator=(Matrix &b)
{
    memcpy(b.components, this->components, 16 * sizeof(double));

    return *this;
}

Matrix Matrix::operator*=(Matrix &b)
{
    Matrix newMat;

	(newMat)[0] = b[0] * components[0] + b[4] * components[1]
			+ b[8] * components[2] + b[12] * components[3];
	(newMat)[1] = b[1] * components[0] + b[5] * components[1]
			+ b[9] * components[2] + b[13] * components[3];
	(newMat)[2] = b[2] * components[0] + b[6] * components[1]
			+ b[10] * components[2] + b[14] * components[3];
	(newMat)[3] = b[3] * components[0] + b[7] * components[1]
			+ b[11] * components[2] + b[15] * components[3];

	(newMat)[4] = b[0] * components[4] + b[4] * components[5]
			+ b[8] * components[6] + b[12] * components[7];
	(newMat)[5] = b[1] * components[4] + b[5] * components[5]
			+ b[9] * components[6] + b[13] * components[7];
	(newMat)[6] = b[2] * components[4] + b[6] * components[5]
			+ b[10] * components[6] + b[14] * components[7];
	(newMat)[7] = b[3] * components[4] + b[7] * components[5]
			+ b[11] * components[6] + b[15] * components[7];

	(newMat)[8] = b[0] * components[8] + b[4] * components[9]
			+ b[8] * components[10] + b[12] * components[11];
	(newMat)[9] = b[1] * components[8] + b[5] * components[9]
			+ b[9] * components[10] + b[13] * components[11];
	(newMat)[10] = b[2] * components[8] + b[6] * components[9]
			+ b[10] * components[10] + b[14] * components[11];
	(newMat)[11] = b[3] * components[8] + b[7] * components[9]
			+ b[11] * components[10] + b[15] * components[11];

	(newMat)[12] = b[0] * components[12] + b[4] * components[13]
			+ b[8] * components[14] + b[12] * components[15];
	(newMat)[13] = b[1] * components[12] + b[5] * components[13]
			+ b[9] * components[14] + b[13] * components[15];
	(newMat)[14] = b[2] * components[12] + b[6] * components[13]
			+ b[10] * components[14] + b[14] * components[15];
	(newMat)[15] = b[3] * components[12] + b[7] * components[13]
			+ b[11] * components[14] + b[15] * components[15];

	memcpy(this->components, newMat.components, 16 * sizeof(double));

	return *this;
}

void Matrix::multiply(Matrix &b)
{
	(*this) *= b;
}

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
	components[0] *= a;
	components[5] *= b;
	components[10] *= c;
}

void Matrix::rotateHK(double hRot, double kRot)
{
    double hRad = hRot * M_PI / 180;
    double kRad = kRot * M_PI / 180;
    
    this->rotate(hRad, kRad, 0);
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

	Matrix alphaThis = *alphaMat * *this;
	Matrix betaAlphaThis = *betaMat * alphaThis;
	Matrix allThis = *gammaMat * betaAlphaThis;

	memcpy(this->components, allThis.components, 16 * sizeof(double));

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

	this->multiply(matrix);
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

Matrix Matrix::inverse3DMatrix()
{
    double determinant = 0;
    determinant = components[0] * components[5] * components[10] +
    components[1] * components[6] * components[8] +
    components[2] * components[4] * components[9] -
    components[0] *components[6] * components[9] -
    components[1] * components[4] * components[10] -
    components[2] * components[5] * components[8];
    
    if (determinant == 0)
        return Matrix();
    
    Matrix inverse = Matrix();
    
    inverse[0] = invertValue(components[5], components[10], components[9], components[6]);
    inverse[1] = invertValue(components[4], components[10], components[8], components[6]);
    inverse[2] = invertValue(components[4], components[9], components[8], components[5]);

    inverse[4] = invertValue(components[1], components[10], components[9], components[2]);
    inverse[5] = invertValue(components[0], components[10], components[8], components[2]);
    inverse[6] = invertValue(components[0], components[9], components[8], components[1]);

    inverse[8] = invertValue(components[1], components[6], components[5], components[2]);
    inverse[9] = invertValue(components[0], components[6], components[4], components[2]);
    inverse[10] = invertValue(components[0], components[5], components[4], components[1]);
    
    inverse[1] = 0 - inverse[1];
    inverse[4] = 0 - inverse[4];
    inverse[6] = 0 - inverse[6];
    inverse[9] = 0 - inverse[9];
    
    for (int i = 0; i < 16; i += 4)
    {
        for (int j = 0; j < 3; j++)
        {
            inverse[i + j] /= determinant;
        }
    }
    
    return inverse;
}

Matrix Matrix::transpose()
{
    Matrix transpose = Matrix();
    
    transpose[0] = components[0];
    transpose[1] = components[4];
    transpose[2] = components[8];
    transpose[3] = components[12];
    transpose[4] = components[1];
    transpose[5] = components[5];
    transpose[6] = components[9];
    transpose[8] = components[2];
    transpose[9] = components[6];
    transpose[10] = components[10];
    transpose[11] = components[14];
    transpose[12] = components[3];
    transpose[13] = components[7];
    transpose[14] = components[11];
    transpose[15] = components[15];
    
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
