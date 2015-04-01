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
#include <cctbx/miller.h>

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
    
    
public:
    double components[16];
    
    Matrix(void);
    Matrix(double *components);
    MatrixPtr copy(void);
    void printDescription();
    std::string description();
    Matrix inverse2DMatrix();
    Matrix inverse3DMatrix();
    Matrix transpose();
    Matrix operator=(Matrix &b);
    cctbx::miller::index<> multiplyIndex(cctbx::miller::index<> *index);

    void translate(double x, double y, double z);
    void rotateHK(double hRot, double kRot);
    void rotate(double alpha, double beta, double gamma);
    void rotateRoundUnitVector(double *unitVector, double radians);
    void multiply(Matrix &b);
    void multiplyVector(vec *vector);
    void preMultiply(Matrix &b);
    void scale(double scale);
    void scale(double a, double b, double c);
    void identity(void);
    void rotateModelAxes(double alpha, double beta, double gamma);
    void newMultiplyVector(double *vector[]);
    static Matrix matrixFromUnitCell(double a, double b, double c, double alpha, double beta, double gamma);
    void orientationMatrixUnitCell(double *a, double *b, double *c);
    void changeOrientationMatrixDimensions(double newA);
    
    void rotate2D(double angle);
    void translation(double **vector);
    
    double getEwaldSphere(vec *vector);
    double getEwaldSphereNoMatrix(vec index);
    
    double *array();
    void print(void);
    
    Matrix operator*=(Matrix &b);
    Matrix operator*(Matrix &b);
    double &operator[](int index) {return components[index]; };
};

#endif /* defined(__RaddoseViewer__Matrix__) */
