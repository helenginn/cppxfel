//
//  Scaler.cpp
//  cppxfel
//
//  Created by Helen Ginn on 01/04/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "Scaler.h"
#include "parameters.h"
#include <scitbx/lbfgsb.h>

double Scaler::stepForParam(int paramNum)
{
    if (paramNum == PARAM_SCALE_FACTOR)
        return 0.001;
    
    if (paramNum == PARAM_B_FACTOR)
        return 0.001;
    
    if (paramNum == PARAM_WAVELENGTH)
        return 0.00001;
    
    if (paramNum == PARAM_HROT || paramNum == PARAM_KROT)
        return 0.00001;
    
    return false;
}

bool Scaler::isRefiningParameter(int paramNum)
{
    if (paramNum == PARAM_SCALE_FACTOR)
        return true;
    
    if (paramNum == PARAM_B_FACTOR)
        return true;
    
    if (paramNum == PARAM_WAVELENGTH)
        return false;
    
    if (paramNum == PARAM_HROT || paramNum == PARAM_KROT)
        return false;
    
    return false;
}

double Scaler::gNorm(int paramNum)
{
    double count = 0;
    double Gsquares = 0;
    
    for (int i = paramNum; i < parameterCount(); i += paramsPerImage())
    {
        Gsquares += g[i] * g[i];
        count++;
    }
    
    return pow(Gsquares, 0.5);
}

double Scaler::hessianGradientForParam(MtzPtr mtz, int paramNum)
{/*
    double *params = new double[PARAM_NUM];
    mtz->getParams(&params);
    
    double currentParam = params[paramNum];
    double paramPlus = currentParam + stepForParam(paramNum);
    double paramMinus = currentParam - stepForParam(paramNum);
    
    params[paramNum] = paramMinus;
    mtz->setParams(params);
    mtz->refreshPartialities(params);
    double below = gradientForImageParameter(mtz, paramNum);
    
    params[paramNum] = paramPlus;
    mtz->setParams(params);
    mtz->refreshPartialities(params);
    double above = gradientForImageParameter(mtz, paramNum);
    
    double second_derivative = (above - below) / (paramPlus - paramMinus);
    
    params[paramNum] = currentParam;
    mtz->setParams(params);
    mtz->refreshPartialities(params);
    
    delete [] params;
    
    if (second_derivative == 0)
        return 1;
    
    return abs(1 / second_derivative) * hessianGradientForParamType(paramNum);*/
    
    return hessianGradientForParamType(paramNum);
}

double Scaler::hessianGradientForParamType(int paramNum)
{
    // bigger the hessian the greater the movement
    
    if (paramNum == PARAM_SCALE_FACTOR)
        return 1;

    if (paramNum == PARAM_WAVELENGTH)
        return 0.00001;
    
    if (paramNum == PARAM_HROT || paramNum == PARAM_KROT)
        return 1;

    if (paramNum == PARAM_B_FACTOR)
        return 100000;
    
    return 1;
}

Scaler::Scaler(std::vector<MtzPtr> mtzs, MtzManager **grouped)
{
    groupedMtz = grouped;

    for (int i = 0; i < mtzs.size(); i++)
    {
        std::vector<Holder *> refHolders, imageHolders;
        mtzs[i]->findCommonReflections(*grouped, imageHolders, refHolders);
        
        mtzData[mtzs[i]] = refHolders;
    }
}

double Scaler::evaluateForImage(MtzPtr mtz)
{
    double numerator = 0;
    double denominator = 0;
    
    for (int i = 0; i < mtzData[mtz].size(); i++)
    {
        mtzData[mtz][i]->rMergeContribution(&numerator, &denominator);
    }
    
    if (denominator == 0)
        return 0;
    
    double fx = numerator / denominator;
    
    if (fx != fx)
        return 0;
    
    return fx;
}

double Scaler::gradientForImageParameter(MtzPtr mtz, int paramNum)
{
    bool needsRefresh = paramNum < 7;
    
    double *params = new double[PARAM_NUM];
    mtz->getParams(&params);
    
    double currentParam = params[paramNum];
    double paramPlus = currentParam + 0.0001;
    double paramMinus = currentParam - 0.0001;
    
    params[paramNum] = paramMinus;
    mtz->setParams(params);
    if (needsRefresh)
        mtz->refreshPartialities(params);
    double below = evaluateForImage(mtz);

    params[paramNum] = paramPlus;
    mtz->setParams(params);
    if (needsRefresh)
        mtz->refreshPartialities(params);
    double above = evaluateForImage(mtz);
    
    double gradient = (above - below) / (paramPlus - paramMinus);
  
    params[paramNum] = currentParam;
    mtz->setParams(params);
    if (needsRefresh)
        mtz->refreshPartialities(params);

    delete [] params;
    
    return gradient;
}

double Scaler::evaluate()
{
    double numerator = 0;
    double denominator = 0;
    
    for (int i = 0; i < (*groupedMtz)->holderCount(); i++)
    {
        (*groupedMtz)->holder(i)->rMergeContribution(&numerator, &denominator);
    }
    
    double fx = numerator / denominator;
    
    return fx;
}

int Scaler::paramsPerImage()
{
    return PARAM_NUM;
}

int Scaler::parameterCount()
{
    return PARAM_NUM * (int)mtzData.size();
}

void Scaler::loadParametersFromMtzs()
{
    int i = 0;
    int paramNum = paramsPerImage();
    
    for (MtzDataMap::iterator it = mtzData.begin(); it != mtzData.end(); it++)
    {
        double *params = new double[paramNum];
        MtzPtr mtz = it->first;
        mtz->getParams(&params);
        
        for (int j = 0; j < paramNum; j++)
        {
            x[i * paramNum + j] = params[j];
        }
    
        i++;
        delete [] params;
    }
}

void Scaler::loadParametersIntoMtzs()
{
    int i = 0;
    int paramNum = paramsPerImage();
    
    for (MtzDataMap::iterator it = mtzData.begin(); it != mtzData.end(); it++)
    {
        double *params = new double[paramNum];
        MtzPtr mtz = it->first;
        
        for (int j = 0; j < paramNum; j++)
        {
            params[j] = x[i * paramNum + j];
        }
        
        mtz->setParams(params);
        
        i++;
        delete [] params;
    }
}

void Scaler::calculateDiagonals()
{
    int i = 0;
    int paramNum = paramsPerImage();
    
    for (MtzDataMap::iterator it = mtzData.begin(); it != mtzData.end(); it++)
    {
        MtzPtr mtz = it->first;
        
        for (int j = 0; j < paramNum; j++)
        {
            int num = i * paramNum + j;
            
            if (isRefiningParameter(j))
            {
                double diagonal = hessianGradientForParam(mtz, j);
                diag[num] = diagonal;
            }
        }

        i++;
    }
}

void Scaler::calculateGradients()
{
    int i = 0;
    int paramNum = paramsPerImage();
    
    for (MtzDataMap::iterator it = mtzData.begin(); it != mtzData.end(); it++)
    {
        MtzPtr mtz = it->first;
        
        for (int j = 0; j < paramNum; j++)
        {
            double gradient = 0;
            int num = i * paramNum + j;
            
            if (isRefiningParameter(j))
                gradient = gradientForImageParameter(mtz, j);
            
            g[num] = gradient;
            
        }
        
        i++;
    }
}

double Scaler::xStepNorm(double step, int paramNum)
{
    double stepSquares = 0;
    double count = 0;
    
    for (int i = paramNum; i < parameterCount(); i += paramsPerImage())
    {
        double diagonal = diag[i];
        double gradient = g[i];
        double step = diagonal * gradient;
        
        double addition = pow(step, 2);
        stepSquares += addition;
        count++;
    }
    
    return pow(stepSquares, 0.5);
}

void Scaler::minimizeRMerge()
{
    int n = parameterCount();
    
    ostringstream logged;
  /*
    scitbx::af::shared<int> nbd(n);
    scitbx::af::shared<double> lowerLims(n);
    scitbx::af::shared<double> upperLims(n);*/
    
    x = scitbx::af::shared<double>(n);
    scitbx::af::shared<double> bestX(n);
    diag = scitbx::af::shared<double>(n);
    g = scitbx::af::shared<double>(n);
    
    double *gb = &(*g.begin());
    double *xb = &(*x.begin());
    double *diagb = &(*diag.begin());
    
    logged << "Number of variables: " << x.size() << std::endl;
/*
    for (int i = 0; i < n; i++)
        nbd[i] = 0;
    
    double factr = 1.e+7;
    double pgtol = 0;
    int iprint = 0;*/
    
    scitbx::lbfgs::minimizer<double> *minimizer =
    new scitbx::lbfgs::minimizer<double>(n);
    
    scitbx::lbfgs::traditional_convergence_test<double, int> is_converged(n);
    
    loadParametersFromMtzs();
   
    Logger::mainLogger->addStream(&logged);
    
    for (int i = 0; i < n; i++)
    {
        g[i] = 0;
        diag[i] = 1;
    }
    
    calculateDiagonals();
    
    int count = 0;
    double bestFx = FLT_MAX;
    
    try
    {
        for (int num = 0;; num++)
        {
            double fx = evaluate(); // calculate fx
            
            // gradient values (set to 0 initially, then populate arrays)
            
            calculateGradients();
            
           /* for (int i = 0; i < mtzData.size(); i++)
            {
                int num = i * paramsPerImage() + PARAM_SCALE_FACTOR;
            }*/
            
            double step = minimizer->stp();
            
            std::cout << "fx: " << fx << "; step: " << step << std::endl;
            
            std::cout << "Gnorms: ";
            
            for (int i = 0; i < paramsPerImage(); i++)
            {
                std::cout << gNorm(i) << ", ";
            }
            /*
            std::cout << std::endl << "Xnorms: ";
            
            for (int i = 0; i < paramsPerImage(); i++)
            {
                std::cout << xStepNorm(step, i) << ", ";
            }
            
            std::cout << std::endl;*/
            
            if (fx < bestFx)
            {
                bestX = x;
                bestFx = fx;
            }
            
            if (minimizer->requests_diag() || count == 0)
            {
                std::cout << "Supplying diagonals" << std::endl;
                minimizer->run(xb, fx, gb, diagb);
            }
            else
                minimizer->run(xb, fx, gb);
            
            loadParametersIntoMtzs();
            
            count++;
            
            if (minimizer->iter() > 50 || count > 50)
                break;
        }
    } catch (scitbx::lbfgs::error &e)
    {
        std::cout << e.what() << std::endl;
    } catch (void *e)
    {
        std::cout << "Unknown error!" << std::endl;
    }
    
    x = bestX;
    loadParametersIntoMtzs();
    
    double fx = evaluate();
    std::cout << "Final fx: " << fx << std::endl;
    
    delete minimizer;
}