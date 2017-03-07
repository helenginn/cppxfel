//
//  AmbiguityBreaker.cpp
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "MtzGrouper.h"
#include "MtzMerger.h"
#include "AmbiguityBreaker.h"
#include "StatisticsManager.h"
#include "FileParser.h"
#include <vector>

/*
double AmbiguityBreaker::dotProduct(int imageNumI, int imageNumJ)
{
    int arrayStartI = imageNumI * ambiguityCount;
    int arrayStartJ = imageNumJ * ambiguityCount;
    
    double fx = 0;
    
    for (int i = 0; i < ambiguityCount; i++)
    {
        double a = x[arrayStartI + i];
        double b = x[arrayStartJ + i];
        
        double addition = a * b;
        
        fx += addition;
    }
    
    return fx;
}

double AmbiguityBreaker::evaluation()
{
    double fx = 0;
    
    int imageCount = (int)mtzs.size();
    
    for (int i = 0; i < imageCount - 1; i++)
    {
        for (int j = i + 1; j < imageCount; j++)
        {
            double correl = gridCorrelation(i, j);

            if (correl == -1)
                continue;
            
            double dot = dotProduct(i, j);
            
            fx += pow(correl - dot, 2);
        }
    }
    
    return fx;
}

double AmbiguityBreaker::gradientForImage(int imageNum, int axis)
{
    double g = 0;
    
    for (int j = 0; j < mtzs.size(); j++)
    {
        if (j == imageNum)
            continue;
        
        double correl = gridCorrelation(imageNum, j);
        
        if (correl == -1)
            continue;

        int arrayStartJ = j * ambiguityCount;
        
        double dot = dotProduct(imageNum, j);
        
        double jAxis = x[arrayStartJ + axis];
        
        g += (-2) * jAxis * (correl - dot);
    }
    
    return g;
}

double AmbiguityBreaker::gridCorrelation(int imageNumI, int imageNumJ)
{
    return statsManager->gridCorrelation(imageNumI, imageNumJ);
}*/

void AmbiguityBreaker::calculateCorrelations(AmbiguityBreaker *me, int offset)
{
    int maxThreads = 1;//FileParser::getMaxThreads();
    std::ostringstream logged;

    for (int i = offset; i < me->mtzs.size(); i += maxThreads)
    {
        MtzPtr iMtz = me->mtzs[i];
        iMtz->setActiveAmbiguity(0);
        
        for (int j = 0; j < i; j++)
        {
            MtzPtr jMtz = me->mtzs[j];
            
            for (int k = 0; k < me->ambiguityCount; k++)
            {
                int hits = 0;
                double cc = StatisticsManager::cc_pearson(&*iMtz, &*jMtz, true, &hits);
                
                me->correlationMaps[k][iMtz][jMtz] = cc;
                me->correlationMaps[k][jMtz][iMtz] = cc;
                iMtz->incrementActiveAmbiguity();
            }
            
            iMtz->setActiveAmbiguity(0);
        }
    }
    
    if (offset == 0)
    {
        logged << std::endl;
        Logger::log(logged);
    }
    
}

void AmbiguityBreaker::makeCorrelationGrid()
{
    logged << "************************************" << std::endl;
    logged << "****  Calculating correlations  ****" << std::endl;
    logged << "************************************" << std::endl << std::endl;
    logged << "There are " << mtzs.size() << " mtz files." << std::endl;
    sendLog();
    
    
    //int maxThreads = FileParser::getMaxThreads();
    
    for (int j = 0; j < ambiguityCount; j++)
    {
        correlationMaps.push_back(CorrelationMap());
    }
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        for (int j = 0; j < ambiguityCount; j++)
        {
            correlationMaps[j][mtzs[i]] = CrystalCorrelation();
        }
    }
    
    boost::thread_group threads;
    for (int i = 0; i < 1; i++)
    {
        boost::thread *thr = new boost::thread(calculateCorrelations, this, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
}

// Call constructor and then run()

AmbiguityBreaker::AmbiguityBreaker(vector<MtzPtr> newMtzs)
{
    setMtzs(newMtzs);
}

void AmbiguityBreaker::setMtzs(vector<MtzPtr> newMtzs)
{
    mtzs = newMtzs;
    
    if (!mtzs.empty())
    {
        ambiguityCount = mtzs[0]->ambiguityCount();
    }
    else
    {
        logged << "No MTZs loaded; not breaking indexing ambiguity." << std::endl;
        sendLogAndExit();
    }
    
    bool shouldApplyUnrefinedPartiality = FileParser::getKey("APPLY_UNREFINED_PARTIALITY", false);
    
    if (shouldApplyUnrefinedPartiality)
    {
        logged << "Applying unrefined partiality to images." << std::endl;
        sendLog();
        
        for (int i = 0; i < mtzs.size(); i++)
        {
            mtzs[i]->applyUnrefinedPartiality();
        }
    }
    else
    {
        logged << "Merging without applying any starting partiality model." << std::endl;
        sendLog();
    }
    
    bool trustAnyway = FileParser::getKey("TRUST_INDEXING_SOLUTION", false);
    
    if (trustAnyway)
    {
        logged << "TRUST_INDEXING_SOLUTION is set to ON so ambiguity breaking is skipped." << std::endl;
        sendLog();
        
        return;
    }

    
    if (ambiguityCount > 1)
    {
        makeCorrelationGrid();
    }
}


void AmbiguityBreaker::run()
{
    bool trustAnyway = FileParser::getKey("TRUST_INDEXING_SOLUTION", false);
    
    if (!trustAnyway && ambiguityCount > 1)
    {
        breakAmbiguity();
        printResults();
    }
    merge();
}

void AmbiguityBreaker::merge()
{   
    MtzGrouper *grouper = new MtzGrouper();
    grouper->setScalingType(ScalingTypeAverage);
    grouper->setWeighting(WeightTypePartialitySigma);
    grouper->setMtzManagers(mtzs);
    
    grouper->setCutResolution(false);
    
    grouper->setExcludeWorst(false);
    
    MtzPtr mergedMtz, unmergedMtz;
    grouper->merge(&mergedMtz, &unmergedMtz, -1, false);
    merged = mergedMtz;
    mergedMtz->writeToFile("remerged.mtz");
    
    delete grouper;
}

void AmbiguityBreaker::printResults()
{
    int ambiguityNums[8];
    memset(ambiguityNums, 0, sizeof(int) * 8);
    
    for (int j = 0; j < mtzs.size(); j++)
    {
        ambiguityNums[mtzs[j]->getActiveAmbiguity()]++;
    }
    
    for (int i = 0; i < ambiguityCount; i++)
    {
        logged << "Ambiguity " << i << " : " << ambiguityNums[i] << " crystals." << std::endl;
        sendLog();
    }
}

void AmbiguityBreaker::breakAmbiguity()
{
    bool trustAnyway = FileParser::getKey("TRUST_INDEXING_SOLUTION", false);
    
    if (trustAnyway)
    {
        return;
    }
    
    logged << "***********************************" << std::endl;
    logged << "*** Breaking indexing ambiguity ***" << std::endl;
    logged << "***********************************" << std::endl << std::endl;
    sendLog();

    for (int i = 0; i < mtzs.size(); i++)
    {
        int random = rand() % (ambiguityCount);
        
        mtzs[i]->setActiveAmbiguity(random);
    }
    
    bool unchanged = false;
    int cycles = 0;
    
    while (!unchanged)
    {
        cycles++;
        unchanged = true;
        int changedNum = 0;
        double totalMaxAverage = 0;
        
        for (int i = 0; i < mtzs.size(); i++)
        {
            int myAmbiguity = mtzs[i]->getActiveAmbiguity();
            
            double ccSums[8];
            int ccCounts[8];
            
            memset(ccSums, 0, sizeof(double) * 8);
            memset(ccCounts, 0, sizeof(int) * 8);
            
            for (int j = 0; j < mtzs.size(); j++)
            {
                if (i == j) continue;
                
                int theirAmbiguity = mtzs[j]->getActiveAmbiguity();
                double chosenCC = correlationMaps[theirAmbiguity][mtzs[i]][mtzs[j]];
                
                if (chosenCC < 0)
                    continue;
                
                ccSums[theirAmbiguity] += chosenCC;
                ccCounts[theirAmbiguity]++;
            }
            
            double maxAve = 0;
            int maxAveValue = 0;
            
            for (int j = 0; j < ambiguityCount; j++)
            {
                ccSums[j] /= ccCounts[j];
                
                if (ccSums[j] >= maxAve)
                {
                    maxAve = ccSums[j];
                    maxAveValue = j;
                }
            }
            
            totalMaxAverage += maxAve;
            
            if (maxAveValue != myAmbiguity)
            {
                mtzs[i]->setActiveAmbiguity(maxAveValue);
                unchanged = false;
                changedNum++;
            }
        }
        
        totalMaxAverage /= mtzs.size();
        
        logged << "Cycle " << cycles << " - switched " << changedNum << " indexing ambiguities. Best average now " << totalMaxAverage << std::endl;
        sendLog();
    }
    
    /*
    int n = ambiguityCount * (int)mtzs.size();

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
    
    for (int i = 0; i < n; i++)
        x[i] = 1;
    
    for (int i = 0; i < x.size(); i++)
    {
        double random = (double)rand() / RAND_MAX;
        x[i] = random;
    }
    
    sendLog();
    
    scitbx::af::shared<double> g(n);
    scitbx::af::ref<double> gb(g.begin(), n);
    
    for (int i = 0; i < n; i++)
        g[i] = 0;
    
    try
    {
        for (int num = 0;; num++)
        {
            float fx = evaluation(); // calculate fx
            
            // gradient values (set to 0 initially, then populate arrays)
            
            for (int i = 0; i < n; i++)
            {
                g[i] = 0;
            }
            
            for (int i = 0; i < mtzs.size(); i++)
            {
                for (int j = 0; j < ambiguityCount; j++)
                {
                    g[i * ambiguityCount + j] = gradientForImage(i, j);
                }
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
    }*/
}