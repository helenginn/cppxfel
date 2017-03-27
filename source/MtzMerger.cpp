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
#include "StatisticsManager.h"

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
    if (FileParser::hasKey("MERGE_TO_RESOLUTION"))
    {
        /* Return this resolution if it exists - user preference takes precedence */
        return FileParser::getKey("MERGE_TO_RESOLUTION", 1.4);
    }
    
    double maxRes = 1 / 2.0;
    
    if (!lowMemoryMode)
    {
        for (int i = 0; i < allMtzs.size(); i++)
        {
            double thisRes = allMtzs[i]->maxResolution();
            
            if (thisRes > maxRes)
            {
                maxRes = thisRes;
            }
        }
    }
    else if (lowMemoryMode)
    {
        /* Return a default, any default! */
        maxRes = 1.4;
    }
    
    return 1 / maxRes;
}

// MARK: write type of MTZ

void MtzMerger::createUnmergedMtz()
{
    float cell[6], wavelength, fdata[9];
    int num = 0;
    
    /* variables for symmetry */
    CCP4SPG *mtzspg = mergedMtz->getSpaceGroup();
    float rsm[192][4][4];
    char ltypex[2];
    
    /* variables for MTZ data structure */
    MTZ *mtzout;
    MTZXTAL *xtal;
    MTZSET *set;
    MTZCOL *colout[5];
    
    /*  Removed: General CCP4 initializations e.g. HKLOUT on command line */
    
	std::vector<double> aCell = mergedMtz->getUnitCell();
    
    for (int i = 0; i < 6; i++)
    {
        cell[i] = (float)aCell[i];
    }
    
	wavelength = 0;

    std::string fullPath = FileReader::addOutputDirectory("u_" + filename);
    
    mtzout = MtzMalloc(0, 0);
    ccp4_lwtitl(mtzout, "Unmerged dataset ", 0);
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
    colout[3] = MtzAddColumn(mtzout, set, "I", "J");
    colout[4] = MtzAddColumn(mtzout, set, "SIGI", "Q");

    for (int i = 0; i < mergedMtz->reflectionCount(); i++)
    {
        ReflectionPtr refl = mergedMtz->reflection(i);
        MillerPtr miller = refl->miller(0);
        int h = miller->getH();
        int k = miller->getK();
        int l = miller->getL();
        
        for (int j = 0; j < refl->liteMillerCount(); j++)
        {
            LiteMiller lite = refl->liteMiller(j);
            double meanIntensity = lite.intensity;
            double meanSigma = lite.weight;
            
            if (meanSigma != meanSigma || meanIntensity != meanIntensity)
            {
                continue;
            }
            
            num++;
            
            fdata[0] = h;
            fdata[1] = k;
            fdata[2] = l;
            fdata[3] = meanIntensity;
            fdata[4] = meanSigma;
            ccp4_lwrefl(mtzout, fdata, colout, 5, num);
        }
    }
    
    // print header information, just for info
    //	ccp4_lhprt(mtzout, 1);
    MtzPut(mtzout, " ");
    MtzFree(mtzout);
}

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
    float cell[6], wavelength, fdata[9];
    int num = 0;

    /* variables for symmetry */
    CCP4SPG *mtzspg = mean->getSpaceGroup();
    float rsm[192][4][4];
    char ltypex[2];
    
    /* variables for MTZ data structure */
    MTZ *mtzout;
    MTZXTAL *xtal;
    MTZSET *set;
    MTZCOL *colout[9];
    
    /*  Removed: General CCP4 initializations e.g. HKLOUT on command line */
    
	std::vector<double> aCell = mean->getUnitCell();

    for (int i = 0; i < 6; i++)
    {
        cell[i] = (float)aCell[i];
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
    /*
        if (ccp4spg_is_centric(mtzspg, h, k, l))
        {
            continue;
        }*/
        
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
    if (!excludeWorst)
    {
        return MtzRejectionNotRejected;
    }
    
    if (mtz->getRefCorrelation() < correlationThreshold)
    {
        return MtzRejectionCorrelation;
    }
    
    double rSplit = mtz->getLastRSplit();
    
    if (rSplit > rFactorThreshold)
    {
        return MtzRejectionRFactor;
    }
    
    double partCorrel = mtz->getRefPartCorrel();
    
    if (partCorrel < partCorrelThreshold)
    {
        return MtzRejectionPartCorrel;
    }
    
    double scale = mtz->getScale();
    double rejectBelow = FileParser::getKey("REJECT_BELOW_SCALE", 0.0);
    
    if (scale < rejectBelow)
    {
        return MtzRejectionOther;
    }
    
    if (mtz->accepted() < minReflectionCounts)
    {
        return MtzRejectionMinRefl;
    }
    
    if (mtz->isRejected())
    {
        return MtzRejectionOther;
    }
    
    return MtzRejectionNotRejected;
}

bool MtzMerger::mtzIsPruned(MtzPtr mtz)
{
    MtzRejectionReason rejectReason = isMtzAccepted(mtz);
    
    std::lock_guard<std::mutex> lg(*rejectMutex);
    
    logged << "Rejecting " << mtz->getFilename() << " " << rejectReason << std::endl;
    sendLog(LogLevelDetailed);
    
    if (rejectNums.count(rejectReason) == 0)
    {
        rejectNums[rejectReason] = 0;
    }
    
    rejectNums[rejectReason]++;
    
    return (rejectReason != MtzRejectionNotRejected);
}

void MtzMerger::summary()
{
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
        logged << "N: Total accepted: " << rejectNums[MtzRejectionNotRejected] << std::endl;
        logged << "N: --------------------------" << std::endl;
        
        sendLog();
    }
}

// MARK: scaling

void MtzMerger::scaleIndividual(MtzPtr mtz)
{
    double scale = 1;
    
    if (scalingType == ScalingTypeAverage)
    {
        scale = 1000 / mtz->averageIntensity();
    }
    else if (scalingType == ScalingTypeReference)
    {
        MtzManager *reference = MtzManager::getReferenceManager();
        
        mtz->scaleToMtz(reference, true);
    }
    else if (scalingType == ScalingTypeBFactor)
    {
        double bFactor, scale;
        mtz->bFactorAndScale(&scale, &bFactor);
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


    for (int i = 0; i < allMtzs.size(); i++)
    {
        params << allMtzs[i]->writeParameterSummary() << std::endl;
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

    CCP4SPG *spg = allMtzs[0]->getSpaceGroup();
    MatrixPtr anyMat = allMtzs[0]->getMatrix();
    std::vector<double> unitCell = allMtzs[0]->getUnitCell();
    anyMat->maxMillers(maxMillers, maxRes);
    
    for (int h = - maxMillers[0]; h <= maxMillers[0]; h++)
    {
        for (int k = -maxMillers[1]; k <= maxMillers[1]; k++)
        {
            for (int l = -maxMillers[2]; l <= maxMillers[2]; l++)
            {
                vec hkl = new_vector(h, k, l);
                anyMat->multiplyVector(&hkl);
                
                if (length_of_vector(hkl) > (1 / maxRes))
                {
                    continue;
                }
                
                if (ccp4spg_is_in_asu(spg, h, k, l) && !ccp4spg_is_sysabs(spg, h, k, l))
                {
                    MillerPtr miller = MillerPtr(new Miller(&*whichMtz, h, k, l, false));
                    miller->setRawIntensity(std::nan(" "));

                    ReflectionPtr newReflection = ReflectionPtr(new Reflection());
                    newReflection->setUnitCell(unitCell);
                    newReflection->setSpaceGroup(spg->spg_num);
                    newReflection->addMiller(miller);
                    miller->setParent(newReflection);
                    miller->setMatrix(whichMtz->getMatrix());
                    newReflection->calculateResolution(&*whichMtz);
                    
                    whichMtz->addReflection(newReflection);
                }
            }
        }
    }
    
    
}

void MtzMerger::addMtzMillers(MtzPtr mtz)
{
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
                bool accept = true;
                
                MillerPtr miller = refl->miller(k);
                
                miller->setRejected(RejectReasonMerge, false);
                
                if (miller->isRejected())
                {
                    incrementRejectedReflections();
                    accept = false;
                }
                
                if (!miller->accepted())
                {
                    accept = false;
                }
                
                if (freeOnly && !miller->isFree())
                {
                    accept = false;
                }
                
                if (accept)
                {
                    if (miller->isSpecial())
                    {
                        logged << "Miller " << prettyDesc(miller->getHKL()) << " from " << miller->getMtzParent()->getFilename() << std::endl;
                        logged << "Corrected intensity: " << miller->intensity() << std::endl;
                        logged << "Scale: " << miller->getScale() << std::endl;
                        logged << "Weight: " << miller->getWeight() << std::endl;
                        logged << "Partiality: " << miller->getPartiality() << std::endl;

                        sendLog();
                    }
                    partnerRefl->addLiteMiller(miller);
                }
            }
        }
    }
}

void MtzMerger::groupMillerThread(int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < allMtzs.size(); i += maxThreads)
    {
        MtzPtr mtz = allMtzs[i];
        
        if (lowMemoryMode)
        {
            mtz->loadReflections();
        }
        
        if (mtzIsPruned(mtz))
        {
            continue;
        }
        
        mtz->flipToActiveAmbiguity();
        
        if (needToScale)
        {
            scaleIndividual(mtz);
        }
        
        addMtzMillers(mtz);
        
        if (lowMemoryMode)
        {
            mtz->dropReflections();
        }
        
        mtz->resetFlip();
    }
}

void MtzMerger::groupMillerThreadWrapper(MtzMerger *object, int offset)
{
    object->groupMillerThread(offset);
}

void MtzMerger::groupMillers()
{
    mergedMtz = MtzPtr(new MtzManager());
    mergedMtz->copySymmetryInformationFromManager(allMtzs[0]);
    mergedMtz->setDefaultMatrix();
    
    makeEmptyReflectionShells(mergedMtz);
    rejectNums = std::map<MtzRejectionReason, int>();
    
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
    
    bool mergeMedian = FileParser::getKey("MERGE_MEDIAN", false);
    
    for (int i = offset; i < mergedMtz->reflectionCount(); i += maxThreads)
    {
        double intensity = 0;
        double sigma = 0;
        double countingSigma = 0;
        int rejected = 0;
        
        ReflectionPtr refl = mergedMtz->reflection(i);
        
        if (refl->liteMillerCount() == 0)
        {
            continue;
        }
        
        int *rejPtr = preventRejections ? NULL : &rejected;
        
        if (!mergeMedian)
        {
            refl->liteMerge(&intensity, &countingSigma, &sigma, rejPtr, friedel);
        }
        else
        {
            refl->medianMerge(&intensity, &countingSigma, rejPtr, friedel);
        }
        
        float intFloat = (float)intensity;
    
        if (!std::isfinite(intFloat))
        {
            continue;
        }
        
        // this could be better coded
        for (int r = 0; r < rejected; r++)
        {
            incrementRejectedReflections();
        }
        
        // this should exist. we made it earlier.
        MillerPtr miller = refl->miller(0);
        
        miller->setRawIntensity(intensity);
        miller->setCountingSigma(countingSigma);
        miller->setSigma(sigma);
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
    for (int i = mergedMtz->reflectionCount() - 1; i >= 0 ; i--)
    {
        ReflectionPtr refl = mergedMtz->reflection(i);
        if (!refl->anyAccepted())
        {
			mergedMtz->removeReflection(i);
        }
    }
}

// MARK: fixSigmas

void MtzMerger::fixSigmas()
{
    
    double minRes = 0;
    double maxRes = maxResolution();
    
    std::vector<double> bins;
    StatisticsManager::generateResolutionBins(minRes, maxRes, 20, &bins);
    
    for (int bin = 0; bin < bins.size() - 1; bin++)
    {
		double iSum = 0;
		double sigiSum = 0;
        int reflNum = 0;
        std::vector<MillerPtr> millersToCorrect;
        
        for (int i = 0; i < mergedMtz->reflectionCount(); i++)
        {
            if (!mergedMtz->reflection(i)->betweenResolutions(bins[bin], bins[bin + 1]))
            {
                continue;
            }
            
            MillerPtr miller = mergedMtz->reflection(i)->miller(0);
            
            if (miller->getCountingSigma() > 0)
            {
				iSum += miller->intensity();
				sigiSum += miller->getCountingSigma();
                reflNum++;
            }
            else
            {
                millersToCorrect.push_back(miller);
            }
        }
        
        if (reflNum == 0)
        {
            return;
        }
        
        double iOverSigi = iSum / (double)sigiSum;
        
        logged << "<Isigi> for " << bins[bin] << " to " << bins[bin + 1] << " Å is " << iOverSigi << std::endl;
        sendLog();
        
        for (int i = 0; i < millersToCorrect.size(); i++)
        {
            MillerPtr miller = millersToCorrect[i];
            
            if (miller->getCountingSigma() < -0.5)
            {
                double intensity = miller->intensity();
                miller->setCountingSigma(fabs(intensity / iOverSigi));
            }
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
    rejectMutex = new std::mutex;
    rejectedReflections = 0;
    silent = false;
    friedel = -1;
    freeOnly = false;
    needToScale = true;
    preventRejections = false;
}

// MARK: Things to call from other classes.

void MtzMerger::merge()
{
	if (MtzManager::getReferenceManager())
	{
		double refScale = 1000 / MtzManager::getReferenceManager()->averageIntensity();
		MtzManager::getReferenceManager()->applyScaleFactor(refScale);
	}

    writeParameterCSV();
    
    if (allMtzs.size() <= 1)
    {
        logged << "N: Not enough MTZs, cannot merge." << std::endl;
        return;
    }
    
	mergedMtz = MtzPtr(new MtzManager());
	groupMillers();
    summary();
    size_t imageNum = allMtzs.size();
    
    double rejectsPerImage = (double) rejectedReflections / (double) imageNum;
    
    int observations = totalObservations();
    
    if (needToScale)
    {
        createUnmergedMtz();
    }
    
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

void MtzMerger::copyDetails(MtzMerger &second)
{
    setCycle(second.cycle);
    setFreeOnly(second.freeOnly);
    setNeedToScale(second.needToScale);
    setSilent(true);
    setExcludeWorst(second.excludeWorst);
    setScalingType(second.scalingType);
    setPreventRejections(second.preventRejections);
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

    std::vector<MtzPtr> firstHalfMtzs, secondHalfMtzs;
    splitAllMtzs(firstHalfMtzs, secondHalfMtzs);
    
    MtzMerger firstMerge = MtzMerger();
    firstMerge.setAllMtzs(firstHalfMtzs);
    firstMerge.setFilename(makeFilename("half1Merge"));
    firstMerge.copyDetails(*this);
    anomalous ? firstMerge.mergeAnomalous() : firstMerge.merge();
    MtzPtr idxMerge = firstMerge.getMergedMtz();

    MtzMerger secondMerge = MtzMerger();
    secondMerge.setAllMtzs(secondHalfMtzs);
    secondMerge.setFilename(makeFilename("half2Merge"));
    secondMerge.copyDetails(*this);
    anomalous ? secondMerge.mergeAnomalous() : secondMerge.merge();
    MtzPtr invMerge = secondMerge.getMergedMtz();

    if (!filename.length())
    {
        filename = makeFilename("allMerge");
    }

    setNeedToScale(false);

    anomalous ? mergeAnomalous() : merge();

	if (!mergedMtz)
	{
		return;
	}

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
    firstMerge.copyDetails(*this);
    firstMerge.merge();
    MtzPtr negativeMerge = firstMerge.getMergedMtz();
    
    MtzMerger secondMerge = MtzMerger();
    secondMerge.setAllMtzs(allMtzs);
    secondMerge.setFriedel(1);
    secondMerge.setFilename(makeFilename("tmp2Merge"));
    secondMerge.copyDetails(*this);
    secondMerge.merge();
    MtzPtr positiveMerge = secondMerge.getMergedMtz();
    
    merge();
    writeAnomalousMtz(negativeMerge, positiveMerge, mergedMtz, makeFilename("anomMerge"));
    createAnomalousDiffMtz(negativeMerge, positiveMerge);
}