//
//  MtzMerger.cpp
//  cppxfel
//
//  Created by Helen Ginn on 09/05/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "FileParser.h"
#include "MtzMerger.h"
#include "MtzManager.h"
#include "FileReader.h"
#include <fstream>
#include "misc.h"
#include "Miller.h"
#include "Reflection.h"
#include "csymlib.h"
#include "ccp4_general.h"
#include "ccp4_parser.h"
#include "FreeMillerLibrary.h"

// MARK: Miscellaneous

void MtzMerger::splitAllMtzs(std::vector<MtzPtr> &firstHalfMtzs, std::vector<MtzPtr> &secondHalfMtzs)
{
    int all = (int)allMtzs.size();
    int half = (int)allMtzs.size() / 2;
    
    firstHalfMtzs.reserve(half);
    secondHalfMtzs.reserve(all - half);
    
    firstHalfMtzs.insert(firstHalfMtzs.begin(), allMtzs.begin(), allMtzs.begin() + half);
    secondHalfMtzs.insert(secondHalfMtzs.begin(), allMtzs.begin() + half, allMtzs.end());
    
}

std::string MtzMerger::makeFilename(std::string prefix)
{
    std::string aFilename = prefix + i_to_str(cycle) + ".mtz";
    
    if (cycle < 0)
    {
        aFilename = prefix + ".mtz";
    }
    
    return aFilename;
}

double MtzMerger::maxResolution()
{
    double maxInput = FileParser::getKey("MERGE_TO_RESOLUTION", 1.4);
    
    double maxRes = 0;
    
    if (!lowMemoryMode)
    {
        for (int i = 0; i < someMtzs.size(); i++)
        {
            double thisRes = someMtzs[i]->maxResolution();
            
            if (thisRes > maxRes)
            {
                maxRes = thisRes;
            }
        }
    }
    else if (lowMemoryMode)
    {
        maxRes = maxInput;
    }

    if (maxInput != 0 && maxRes > 1 / maxInput)
    {
        maxRes = maxInput;
    }
    
    return maxRes;
}

// MARK: write anomalous MTZ

void MtzMerger::createAnomalousDiffMtz(MtzPtr negative, MtzPtr positive)
{
    for (int i = 0; i < mergedMtz->reflectionCount(); i++)
    {
        ReflectionPtr refl = mergedMtz->reflection(i);
        MillerPtr meanMiller = refl->miller(0);
        int reflId = (int)refl->getReflId();
        
        ReflectionPtr posRefl;
        ReflectionPtr negRefl;
        
        negative->findReflectionWithId(reflId, &negRefl);
        positive->findReflectionWithId(reflId, &posRefl);
        
        if (!posRefl || !negRefl)
        {
            mergedMtz->removeReflection(i);
            i--;
            continue;
        }
        
        
        double posIntensity = posRefl->miller(0)->intensity();
        double negIntensity = negRefl->miller(0)->intensity();
        double diffIntensity = posIntensity - negIntensity;
        double meanSigma = meanMiller->getCountingSigma();
        
        meanMiller->setRawIntensity(diffIntensity);
        meanMiller->setCountingSigma(meanSigma);
    }
}

void MtzMerger::writeAnomalousMtz(MtzPtr negative, MtzPtr positive, MtzPtr mean, std::string filename)
{
    double doubleCell[6];
    float cell[6], wavelength, fdata[9];
    int num = 0;

    /* variables for symmetry */
    CCP4SPG *mtzspg = mean->getLowGroup();
    float rsm[192][4][4];
    char ltypex[2];
    
    /* variables for MTZ data structure */
    MTZ *mtzout;
    MTZXTAL *xtal;
    MTZSET *set;
    MTZCOL *colout[9];
    
    /*  Removed: General CCP4 initializations e.g. HKLOUT on command line */
    
    mean->getUnitCell(&doubleCell[0], &doubleCell[1], &doubleCell[2], &doubleCell[3], &doubleCell[4],
                             &doubleCell[5]);
    
    for (int i = 0; i < 6; i++)
    {
        cell[i] = (float)doubleCell[i];
    }
    
    wavelength = mean->getWavelength();
    
    std::string fullPath = FileReader::addOutputDirectory(filename);
    
    mtzout = MtzMalloc(0, 0);
    ccp4_lwtitl(mtzout, "Anomalous dataset ", 0);
    mtzout->refs_in_memory = 0;
    mtzout->fileout = MtzOpenForWrite(fullPath.c_str());
    
    // then add symm headers...
    for (int i = 0; i < mtzspg->nsymop; ++i)
        CCP4::rotandtrn_to_mat4(rsm[i], mtzspg->symop[i]);
    strncpy(ltypex, mtzspg->symbol_old, 1);
    ccp4_lwsymm(mtzout, mtzspg->nsymop, mtzspg->nsymop_prim, rsm, ltypex,
                mtzspg->spg_ccp4_num, mtzspg->symbol_old, mtzspg->point_group);
    
    // then add xtals, datasets, cols
    xtal = MtzAddXtal(mtzout, "XFEL crystal", "XFEL project", cell);
    set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
    colout[0] = MtzAddColumn(mtzout, set, "H", "H");
    colout[1] = MtzAddColumn(mtzout, set, "K", "H");
    colout[2] = MtzAddColumn(mtzout, set, "L", "H");
    colout[3] = MtzAddColumn(mtzout, set, "I(+)", "K");
    colout[4] = MtzAddColumn(mtzout, set, "SIGI(+)", "M");
    colout[5] = MtzAddColumn(mtzout, set, "I(-)", "K");
    colout[6] = MtzAddColumn(mtzout, set, "SIGI(-)", "M");
    colout[7] = MtzAddColumn(mtzout, set, "IMEAN", "J");
    colout[8] = MtzAddColumn(mtzout, set, "SIGIMEAN", "Q");
    
    for (int i = 0; i < mean->reflectionCount(); i++)
    {
        ReflectionPtr refl = mean->reflection(i);
        MillerPtr meanMiller = refl->miller(0);
        int reflId = (int)refl->getReflId();
        
        int h = meanMiller->getH();
        int k = meanMiller->getK();
        int l = meanMiller->getL();
        
        if (ccp4spg_is_centric(mtzspg, h, k, l))
        {
            continue;
        }
        
        ReflectionPtr posRefl;
        ReflectionPtr negRefl;
        
        negative->findReflectionWithId(reflId, &negRefl);
        positive->findReflectionWithId(reflId, &posRefl);
        
        if (!posRefl || !negRefl)
        {
            continue;
        }
        
        double meanIntensity = meanMiller->intensity();
        double meanSigma = meanMiller->getCountingSigma();
        
        double posIntensity = posRefl->miller(0)->intensity();
        double posSigma = posRefl->miller(0)->getCountingSigma();

        double negIntensity = negRefl->miller(0)->intensity();
        double negSigma = negRefl->miller(0)->getCountingSigma();

        if (posIntensity != posIntensity || negIntensity != negIntensity || meanIntensity != meanIntensity)
        {
            continue;
        }
        
        if (meanSigma != meanSigma || posSigma != posSigma || negSigma != negSigma)
        {
            continue;
        }
        
        num++;
        
        fdata[0] = h;
        fdata[1] = k;
        fdata[2] = l;
        fdata[3] = negIntensity;
        fdata[4] = negSigma;
        fdata[5] = posIntensity;
        fdata[6] = posSigma;
        fdata[7] = meanIntensity;
        fdata[8] = meanSigma;
        ccp4_lwrefl(mtzout, fdata, colout, 9, num);
    }
    
    // print header information, just for info
    //	ccp4_lhprt(mtzout, 1);
    MtzPut(mtzout, " ");
    MtzFree(mtzout);

}

// MARK: rejecting images.

MtzRejectionReason MtzMerger::isMtzAccepted(MtzPtr mtz)
{
    if (lowMemoryMode)
    {
        mtz->loadReflections(true);
    }
    
    if (!excludeWorst)
    {
        return MtzRejectionNotRejected;
    }
    
    if (mtz->getRefCorrelation() < correlationThreshold)
    {
        return MtzRejectionCorrelation;
    }
    
    double rSplit = mtz->rSplit(0, 0);
    
    if (rSplit > rFactorThreshold)
    {
        return MtzRejectionRFactor;
    }
    
    double partCorrel = mtz->getRefPartCorrel();
    
    if (partCorrel < partCorrelThreshold)
    {
        return MtzRejectionPartCorrel;
    }
    
    if (mtz->accepted() < minReflectionCounts)
    {
        return MtzRejectionMinRefl;
    }
    
    if (mtz->isRejected())
    {
        return MtzRejectionOther;
    }
    
    if (lowMemoryMode)
    {
        mtz->dropReflections();
    }
    
    return MtzRejectionNotRejected;
}

void MtzMerger::pruneMtzs()
{
    someMtzs.clear();
    std::map<MtzRejectionReason, int> rejectNums;
    
    for (int i = 0; i < allMtzs.size(); i++)
    {
        
        MtzRejectionReason rejectReason = isMtzAccepted(allMtzs[i]);
        
        if (rejectReason == MtzRejectionNotRejected)
        {
            someMtzs.push_back(allMtzs[i]);
        }
        else
        {
            logged << "Rejecting " << allMtzs[i]->getFilename() << " " << rejectReason << std::endl;
            sendLog(LogLevelDetailed);
            
            if (rejectNums.count(rejectReason) == 0)
            {
                rejectNums[rejectReason] = 0;
            }
            
            rejectNums[rejectReason]++;
        }
    }
    
    if (!silent)
    {
        logged << "N: --------------------------" << std::endl;
        logged << "N: Crystal rejection summary:" << std::endl;
        logged << "N: --------------------------" << std::endl;
        logged << "N: Correlation too low: " << rejectNums[MtzRejectionCorrelation] << std::endl;
        logged << "N: R factor too high: " << rejectNums[MtzRejectionRFactor] << std::endl;
        logged << "N: Part correl too low: " << rejectNums[MtzRejectionPartCorrel] << std::endl;
        logged << "N: Not enough reflections: " << rejectNums[MtzRejectionMinRefl] << std::endl;
        logged << "N: Some other reason: " << rejectNums[MtzRejectionOther] << std::endl;
        logged << "N: --------------------------" << std::endl;
        logged << "N: Total accepted: " << someMtzs.size() << std::endl;
        logged << "N: --------------------------" << std::endl;
        
        sendLog();
    }
}

// MARK: scaling

void MtzMerger::scaleIndividual(MtzPtr mtz)
{
    double scale = 1;
    
    if (lowMemoryMode)
    {
        mtz->loadReflections(true);
    }
    
    if (scalingType == ScalingTypeAverage)
    {
        scale = 1000 / mtz->averageIntensity();
    }
    else if (scalingType == ScalingTypeReference)
    {
        MtzManager *reference = MtzManager::getReferenceManager();
        
        scale = mtz->gradientAgainstManager(reference, false);
    }
    
    mtz->applyScaleFactor(scale);
}

void MtzMerger::scale()
{
    if (!needToScale)
    {
        return;
    }
    
    for (int i = 0; i < someMtzs.size(); i++)
    {
        scaleIndividual(someMtzs[i]);
    }
}

// MARK: params_cycle_X.csv.

void MtzMerger::writeParameterCSV()
{
    if (silent)
        return;
    
    std::ofstream paramLog;
    std::string paramLogName = "params_cycle_" + i_to_str(cycle) + ".csv";
    
    std::string fullPath = FileReader::addOutputDirectory(paramLogName);
    
    std::ofstream params;
    params.open(fullPath.c_str());
    
    params << MtzManager::parameterHeaders() << std::endl;
    
    if (cycle >= 0)
    {
        for (int i = 0; i < allMtzs.size(); i++)
        {
            params << allMtzs[i]->writeParameterSummary() << std::endl;
        }
    }
    
    params.close();
    logged << "N: --------------------------" << std::endl;
    logged << "N: Written parameter values to " << fullPath << std::endl;
    logged << "N: --------------------------" << std::endl;

    sendLog();
}

// MARK: Grouping millers.

void MtzMerger::incrementRejectedReflections()
{
    reflCountMutex->lock();
    
    rejectedReflections++;
    
    reflCountMutex->unlock();
}

void MtzMerger::makeEmptyReflectionShells(MtzPtr whichMtz)
{
    double maxRes = maxResolution();
    int maxMillers[3];

    CCP4SPG *spg = someMtzs[0]->getLowGroup();
    MatrixPtr anyMat = someMtzs[0]->getMatrix();
    std::vector<double> unitCell = someMtzs[0]->getUnitCell();
    anyMat->maxMillers(maxMillers, maxRes);
    
    for (int h = - maxMillers[0]; h <= maxMillers[0]; h++)
    {
        for (int k = -maxMillers[1]; k <= maxMillers[1]; k++)
        {
            for (int l = -maxMillers[2]; l <= maxMillers[2]; l++)
            {
                if (ccp4spg_is_in_asu(spg, h, k, l) && !ccp4spg_is_sysabs(spg, h, k, l))
                {
                    MillerPtr miller = MillerPtr(new Miller(&*whichMtz, h, k, l, false));
                    miller->setRawIntensity(std::nan(" "));

                    ReflectionPtr newReflection = ReflectionPtr(new Reflection());
                    newReflection->setUnitCellDouble(&unitCell[0]);
                    newReflection->setSpaceGroup(spg->spg_num);
                    newReflection->addMiller(miller);
                    miller->setParent(newReflection);
                    miller->setMatrix(whichMtz->getMatrix());
                    newReflection->calculateResolution(&*whichMtz);
                    
                    if (newReflection->getResolution() > (1 / maxRes))
                    {
                        continue;
                    }
                    
                    whichMtz->addReflection(newReflection);
                }
            }
        }
    }
    
    
}

void MtzMerger::groupMillerThread(int offset)
{
    int wentMissing = 0;
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < someMtzs.size(); i += maxThreads)
    {
        MtzPtr mtz = someMtzs[i];
        mtz->flipToActiveAmbiguity();
        
        if (lowMemoryMode)
        {
            mtz->loadReflections(true);
        }
        
        for (int j = 0; j < mtz->reflectionCount(); j++)
        {
            ReflectionPtr refl = mtz->reflection(j);
            ReflectionPtr partnerRefl;
            int reflId = (int)refl->getReflId();
            
            mergedMtz->findReflectionWithId(reflId, &partnerRefl);
            
            if (partnerRefl)
            {
                for (int k = 0; k < refl->millerCount(); k++)
                {
                    MillerPtr miller = refl->miller(k);
                    
                    miller->setRejected(RejectReasonMerge, false);
                    
                    if (miller->isRejected())
                    {
                        incrementRejectedReflections();
                        continue;
                    }
                    
                    if (!miller->accepted())
                    {
                        continue;
                    }
                    
                    if (freeOnly && !miller->isFree())
                    {
                        continue;
                    }
                    
                    partnerRefl->addLiteMiller(miller);
                }
            }
            else
            {
                wentMissing++;
            }
        }
        
        if (lowMemoryMode)
        {
            mtz->dropReflections();
        }
        
        mtz->resetFlip();
    }
    
    if (wentMissing)
    {
        logged << "N: Warning: " << wentMissing << " millers went missing during merge." << std::endl;
        sendLog(LogLevelDetailed);
    }
}

void MtzMerger::groupMillerThreadWrapper(MtzMerger *object, int offset)
{
    object->groupMillerThread(offset);
}

void MtzMerger::groupMillers()
{
    mergedMtz = MtzPtr(new MtzManager());
    mergedMtz->copySymmetryInformationFromManager(someMtzs[0]);
    mergedMtz->setDefaultMatrix();
    
    makeEmptyReflectionShells(mergedMtz);
    
    boost::thread_group threads;
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(groupMillerThreadWrapper, this, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
}

// MARK: Merging millers.

void MtzMerger::mergeMillersThread(int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < mergedMtz->reflectionCount(); i += maxThreads)
    {
        double intensity = 0;
        double sigma = 0;
        int rejected = 0;
        
        ReflectionPtr refl = mergedMtz->reflection(i);
        
        if (refl->liteMillerCount() == 0)
        {
            continue;
        }
        
        refl->liteMerge(&intensity, &sigma, &rejected, friedel);
     
        // this could be better coded
        for (int r = 0; r < rejected; r++)
        {
            incrementRejectedReflections();
        }
        
        // this should exist. we made it earlier.
        MillerPtr miller = refl->miller(0);
        
        miller->setRawIntensity(intensity);
        miller->setCountingSigma(sigma);
        miller->setSigma(1);
        miller->setPartiality(1);
        
        refl->clearLiteMillers();
    }
}

void MtzMerger::mergeMillersThreadWrapper(MtzMerger *object, int offset)
{
    object->mergeMillersThread(offset);
 
}

void MtzMerger::mergeMillers()
{
    boost::thread_group threads;
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(mergeMillersThreadWrapper, this, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
}

// MARK: remove reflections.

void MtzMerger::removeReflections()
{
    int count = 0;
    
    for (int i = mergedMtz->reflectionCount() - 1; i >= 0 ; i--)
    {
        ReflectionPtr refl = mergedMtz->reflection(i);
        if (!refl->acceptedCount())
        {
            mergedMtz->removeReflection(i);
            count++;
        }
    }
}

// MARK: fixSigmas

void MtzMerger::fixSigmas()
{
    double iOverSigiSum = 0;
    int reflNum = 0;
    
    for (int i = 0; i < mergedMtz->reflectionCount(); i++)
    {
        MillerPtr miller = mergedMtz->reflection(i)->miller(0);
        
        if (miller->getCountingSigma() > 0)
        {
            iOverSigiSum += (miller->intensity() / miller->getCountingSigma());
            reflNum++;
        }
    }
    
    if (reflNum == 0)
    {
        return;
    }
    
    double iOverSigi = iOverSigiSum / (double)reflNum;
    
    for (int i = 0; i < mergedMtz->reflectionCount(); i++)
    {
        MillerPtr miller = mergedMtz->reflection(i)->miller(0);
        
        if (miller->getCountingSigma() == -1)
        {
            double intensity = miller->intensity();
            miller->setCountingSigma(intensity / iOverSigi);
        }
    }
}

int MtzMerger::totalObservations()
{
    int total = 0;
    
    for (int i = 0; i < mergedMtz->reflectionCount(); i++)
    {
        total += mergedMtz->reflection(i)->liteMillerCount();
    }
    
    return total;
}

// MARK: Constructor.

void MtzMerger::setCycle(int num)
{
    cycle = num;

    double corrThreshold = FileParser::getKey("CORRELATION_THRESHOLD", CORRELATION_THRESHOLD);
    double initialCorrelationThreshold = FileParser::getKey("INITIAL_CORRELATION_THRESHOLD", INITIAL_CORRELATION_THRESHOLD);
    int thresholdSwap = FileParser::getKey("THRESHOLD_SWAP", THRESHOLD_SWAP);
    
    if (cycle < thresholdSwap || cycle == -1)
    {
        correlationThreshold = initialCorrelationThreshold;
    }
    else
    {
        correlationThreshold = corrThreshold;
    }
}

MtzMerger::MtzMerger()
{
    rFactorThreshold = FileParser::getKey("R_FACTOR_THRESHOLD", 100.0);
    partCorrelThreshold = FileParser::getKey("PARTIALITY_CORRELATION_THRESHOLD", -2);
    minReflectionCounts = FileParser::getKey("MINIMUM_REFLECTION_CUTOFF", MINIMUM_REFLECTION_CUTOFF);
    lowMemoryMode = FileParser::getKey("LOW_MEMORY_MODE", false);
    excludeWorst = true;
    scalingType = ScalingTypeAverage;
    reflCountMutex = new std::mutex;
    rejectedReflections = 0;
    silent = false;
    friedel = -1;
    freeOnly = false;
    needToScale = true;
}

// MARK: Things to call from other classes.

void MtzMerger::merge()
{
    mergedMtz = MtzPtr(new MtzManager());
    writeParameterCSV();
    pruneMtzs();
    
    if (!needToScale)
    {
        scale();
    }
   

    if (someMtzs.size() == 0)
    {
        logged << "N: Error: No MTZs, cannot merge." << std::endl;
        sendLog();
        return;
    }
    
    groupMillers();
    
    int imageNum = (int)someMtzs.size();
    double rejectsPerImage = (double) rejectedReflections / (double) imageNum;
    
    int observations = totalObservations();
    
    mergeMillers();

    int totalRefls = mergedMtz->reflectionCount();
    
    if (!silent)
    {
        logged << "N: Total observations: " << observations << std::endl;
        logged << "N: Total unique reflections: " << totalRefls << std::endl;
        logged << "N: Multiplicity: " << (double)observations / (double)totalRefls << std::endl;
        logged << "N: Rejects per image: " << rejectsPerImage << std::endl;
    }
    
    removeReflections();
    
    int someRefls = mergedMtz->reflectionCount();
    
    if (!silent)
    {
        logged << "N: Reflections present: " << someRefls << std::endl;
        logged << "N: Completeness: " << 100 * (double)someRefls / (double)totalRefls << "%." << std::endl;
    }
    
    sendLog();
    
    fixSigmas();
    
    mergedMtz->setFilename(filename);
    mergedMtz->writeToFile(filename, !silent);
}

void MtzMerger::mergeFull(bool anomalous)
{
    time_t startcputime;
    time(&startcputime);

    if (freeOnly)
    {
        if (!FreeMillerLibrary::active())
        {
            return;
        }
    }

    scale();
    
    std::vector<MtzPtr> firstHalfMtzs, secondHalfMtzs;
    splitAllMtzs(firstHalfMtzs, secondHalfMtzs);
    
    MtzMerger firstMerge = MtzMerger();
    firstMerge.setAllMtzs(firstHalfMtzs);
    firstMerge.setFilename(makeFilename("half1Merge"));
    firstMerge.setCycle(cycle);
    firstMerge.setFreeOnly(freeOnly);
    firstMerge.setNeedToScale(false);
    firstMerge.setSilent(true);
    firstMerge.setExcludeWorst(excludeWorst);
    anomalous ? firstMerge.mergeAnomalous() : firstMerge.merge();
    MtzPtr idxMerge = firstMerge.getMergedMtz();

    MtzMerger secondMerge = MtzMerger();
    secondMerge.setAllMtzs(secondHalfMtzs);
    secondMerge.setFilename(makeFilename("half2Merge"));
    secondMerge.setExcludeWorst(excludeWorst);
    secondMerge.setFreeOnly(freeOnly);
    secondMerge.setNeedToScale(false);
    secondMerge.setSilent(true);
    secondMerge.setCycle(cycle);
    anomalous ? secondMerge.mergeAnomalous() : secondMerge.merge();
    MtzPtr invMerge = secondMerge.getMergedMtz();

    if (!filename.length())
    {
        filename = makeFilename("allMerge");
    }

    anomalous ? mergeAnomalous() : merge();
    
    double maxRes = 1 / maxResolution();
    
    logged << "N: === R split" << (freeOnly ? " (free)" : "") << " ===" << std::endl;
    sendLog();
    double rSplit = idxMerge->rSplitWithManager(&*invMerge, false, false, 0, maxRes, 20, NULL, true);
    logged << "N: === CC half"  << (freeOnly ? " (free)" : "") << " ===" << std::endl;
    sendLog();
    double correlation = idxMerge->correlationWithManager(&*invMerge, false, false, 0, maxRes, 20, NULL, true);
    
    std::string set = "all";
    
    if (freeOnly)
    {
        set = "free";
    }
    
    if (anomalous)
    {
        set = "anom";
    }
    
    logged << "N: Final stats (" << set << "): " << rSplit << ", " << correlation << std::endl;
    sendLog();
    
    time_t endcputime;
    time(&endcputime);
    
    time_t difference = endcputime - startcputime;
    double seconds = difference;
    
    int finalSeconds = (int) seconds % 60;
    int minutes = seconds / 60;
    
    logged << "N: Clock time " << minutes << " minutes, " << finalSeconds << " seconds to merge (" << set << ")" << std::endl;
    sendLog();
}

void MtzMerger::mergeAnomalous()
{
    MtzMerger firstMerge = MtzMerger();
    firstMerge.setAllMtzs(allMtzs);
    firstMerge.setFriedel(0);
    firstMerge.setFilename(makeFilename("tmp1Merge"));
    firstMerge.setCycle(cycle);
    firstMerge.setNeedToScale(false);
    firstMerge.setSilent(true);
    firstMerge.setExcludeWorst(excludeWorst);
    firstMerge.merge();
    MtzPtr negativeMerge = firstMerge.getMergedMtz();
    
    MtzMerger secondMerge = MtzMerger();
    secondMerge.setAllMtzs(allMtzs);
    secondMerge.setFriedel(1);
    secondMerge.setFilename(makeFilename("tmp2Merge"));
    secondMerge.setExcludeWorst(excludeWorst);
    secondMerge.setNeedToScale(false);
    secondMerge.setSilent(true);
    secondMerge.setCycle(cycle);
    secondMerge.merge();
    MtzPtr positiveMerge = secondMerge.getMergedMtz();
    
    merge();
    writeAnomalousMtz(negativeMerge, positiveMerge, mergedMtz, makeFilename("anomMerge"));
    createAnomalousDiffMtz(negativeMerge, positiveMerge);
}