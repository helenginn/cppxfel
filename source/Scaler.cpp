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

bool Scaler::isRefiningParameter(int paramNum)
{
    if (paramNum == PARAM_SCALE_FACTOR)
        return true;
    
    return false;
}

double Scaler::hessianGradientForParamType(int paramNum)
{
    return 1;
}

Scaler::Scaler(std::vector<MtzPtr> mtzs, MtzManager *grouped)
{
    groupedMtz = grouped;

    for (int i = 0; i < mtzs.size(); i++)
    {
        std::vector<Holder *> refHolders, imageHolders;
        mtzs[i]->findCommonReflections(grouped, imageHolders, refHolders);
        
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
    double *params = new double[PARAM_NUM];
    mtz->getParams(&params);
    
    double currentParam = params[paramNum];
    double paramPlus = currentParam + 0.001;
    double paramMinus = currentParam - 0.001;
    
    params[paramNum] = paramMinus;
    mtz->setParams(params);
    mtz->refreshPartialities(params);
    double below = evaluateForImage(mtz);

    params[paramNum] = paramPlus;
    mtz->setParams(params);
    mtz->refreshPartialities(params);
    double above = evaluateForImage(mtz);
  
    double gradient = (above - below) / (paramPlus - paramMinus);
  
    params[paramNum] = currentParam;
    mtz->setParams(params);
    mtz->refreshPartialities(params);

    delete [] params;
    
    return gradient;
}

double Scaler::evaluate()
{
    double numerator = 0;
    double denominator = 0;
    
    for (int i = 0; i < groupedMtz->holderCount(); i++)
    {
        groupedMtz->holder(i)->rMergeContribution(&numerator, &denominator);
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
            
            if (isRefiningParameter(j))
                gradient = gradientForImageParameter(mtz, j);
            
            g[i * paramNum + j] = gradient;
        }
        
        i++;
    }
}

void Scaler::minimizeRMerge()
{
    int n = parameterCount();
    
    ostringstream logged;
    
    scitbx::af::shared<int> nbd(n);
    scitbx::af::shared<double> lowerLims(n);
    scitbx::af::shared<double> upperLims(n);
    
    x = scitbx::af::shared<double>(n);
    
    logged << "Number of variables: " << x.size() << std::endl;
    
    for (int i = 0; i < n; i++)
        nbd[i] = 0;
    
    double factr = 1.e+7;
    double pgtol = 0;
    int iprint = 0;
    
    scitbx::lbfgsb::minimizer<double> *minimizer =
    new scitbx::lbfgsb::minimizer<double>(n, 5, lowerLims, upperLims,
                                          nbd, false, factr, pgtol, iprint);
    
    scitbx::lbfgs::traditional_convergence_test<double, int> is_converged(n);
    
    scitbx::af::ref<double> xb(x.begin(), n);
    
    loadParametersFromMtzs();
   
    Logger::mainLogger->addStream(&logged);
    
    g = scitbx::af::shared<double>(n);
    scitbx::af::ref<double> gb(g.begin(), n);
    
    for (int i = 0; i < n; i++)
        g[i] = 0;
    
    try
    {
        for (int num = 0;; num++)
        {
            float fx = evaluate(); // calculate fx
            
            // gradient values (set to 0 initially, then populate arrays)
            
            calculateGradients();
            
            for (int i = 0; i < mtzData.size(); i++)
            {
                int num = i * parameterCount() + PARAM_SCALE_FACTOR;
                
           //     std::cout << g[num] << std::endl;
            }
            
            std::cout << "fx: " << fx << std::endl;
            
            minimizer->process(xb, fx, gb);
            
            if (minimizer->is_terminated())
                break;
            
            if (minimizer->n_iteration() > 50)
                break;
        }
    } catch (scitbx::lbfgs::error &e)
    {
        std::cout << e.what() << std::endl;
    } catch (void *e)
    {
        std::cout << "Unknown error!" << std::endl;
    }

}