//
//  Matrix.h
//  RaddoseViewer
//
//  Created by Helen Ginn on 19/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#ifndef __RaddoseViewer__Matrix__
#define __RaddoseViewer__Matrix__

#include <iostream>
#include "Vector.h"
#include "parameters.h"
#include "csymlib.h"
#include "LoggableObject.h"

class Matrix
{
private:
    class Proxy
    {
        Matrix &a;
        int idx;
    public:
        Proxy(Matrix &a, int idx) : a(a), idx(idx) {}
        double operator= (double x) { a.components[idx] = x; return a.components[idx]; }
        double operator* (double x) { return a.components[idx] * x; }
        double operator* (Proxy x) { return a.components[idx] * x.a.components[idx]; }
        double operator*= (double x) { a.components[idx] *= x; return a.components[idx]; }
        double operator*= (Proxy x) { a.components[idx] *= x.a.components[idx]; return a.components[idx]; }
        
    };
    
    MatrixPtr unitCell;
    MatrixPtr rotation;
    double eulerA;
    double eulerB;
    double eulerC;
    
    static MatrixPtr identityMatrix;
public:
    double components[16];
    
    // for SPEEDY CALCULATIONS when you need an identity pointer. Don't change it or you'll mess everyone else up.
    static MatrixPtr getIdentityPtr()
    {
        return identityMatrix;
    }
    
    double trace();
    bool isIdentity();
    void multiply(double scalar);
    void subtract(MatrixPtr secondMatrix);
    Matrix(void);
    Matrix(double *components);
    MatrixPtr copy(void);
    void copyComponents(MatrixPtr mat2);
    void printDescription(bool detailed = false);
    std::string description(bool detailed = false, bool submatrix = false);
    Matrix inverse2DMatrix();
    MatrixPtr inverse3DMatrix();
    static MatrixPtr matFromCCP4(CSym::ccp4_symop *symop);
    MatrixPtr transpose();
	static void symmetryOperatorsForSpaceGroup(std::vector<MatrixPtr> *matrices, CSym::CCP4SPG *spaceGroup, std::vector<double> cell);

    void rotate(double alpha, double beta, double gamma);
    void rotateRoundUnitVector(vec unitVector, double radians);
    void multiply(Matrix &b);
    void multiplyVector(vec *vector);
    void preMultiply(Matrix &b);
    void scale(double scale);
    void scale(double a, double b, double c);
    void identity(void);
	static MatrixPtr matrixFromUnitCell(std::vector<double> &unitCell);
    void orientationMatrixUnitCell(double *a, double *b, double *c);
    void changeOrientationMatrixDimensions(std::vector<double> cell);
    void setComplexMatrix(MatrixPtr unitCell, MatrixPtr rotation);
    void maxMillers(int (&millers)[3], double maxResolution);
    static MatrixPtr randomOrientation();
    static MatrixPtr randomOrientationMatrix();
;
    void rotate2D(double angle);
    
    double getEwaldSphereNoMatrix(vec index);
    
    void eulerAngles(double *theta, double *phi, double *psi, bool force = false);
    double similarityToRotationMatrix(MatrixPtr mat2, double tolerance, bool force = false);
    void unitCellLengths(double *lengths);
    double *array();
    void recalculateOrientationMatrix();
    std::string summary();
    void setIdentity();
    bool writeToHdf5(std::string address);

    bool isComplex()
    {
        if (unitCell)
            return true;
        
        return false;
    }
    
    MatrixPtr getRotation()
    {
        return rotation;
    }
    
    MatrixPtr getUnitCell()
    {
        return unitCell;
    }
    
    void setComponents(double *newComponents)
    {
        memcpy(components, newComponents, sizeof(double) * 16);
    }
    
    double determinant();
    Matrix operator*=(Matrix &b);
    Matrix operator*(Matrix &b);
    Matrix testMultiply(Matrix &b);
    double &operator[](int index) {return components[index]; };
};

#endif /* defined(__RaddoseViewer__Matrix__) */
