/*
 * Holder.cpp
 *
 *  Created on: 25 Oct 2014
 *      Author: helenginn
 */

#include "Reflection.h"

#include <vector>
#include <cmath>
#include "csymlib.h"
#include "Miller.h"

#include "FileParser.h"
#include "StatisticsManager.h"

bool Reflection::hasSetup = false;
bool Reflection::setupUnitCell = false;
unsigned char Reflection::spgNum = 0;
std::mutex Reflection::setupMutex;
std::vector<MatrixPtr> Reflection::flipMatrices;
MatrixPtr Reflection::customAmbiguity;
MatrixPtr Reflection::unitCellMatrix;
double Reflection::rejectSigma = 0;
bool Reflection::shouldReject = 0;
CSym::CCP4SPG *Reflection::ccp4_space_group;

MatrixPtr Reflection::matrixForAmbiguity(int i)
{
    if (i >= ambiguityCount())
    {
        std::cout << "Ambiguity issue!" << std::endl;
    }
    
    if (flipMatrices.size() > i)
    {
        return flipMatrices[i];
    }
    
    if (i == 0)
    {
        MatrixPtr identity = MatrixPtr(new Matrix());
        
        return identity;
    }
    
    if (i == 1)
    {
        if (ambiguityCount() == 2 && customAmbiguity)
        {
            return customAmbiguity;
        }
        
        if (ambiguityCount() == 2 || ambiguityCount() == 4)
        {
            MatrixPtr khMinusL = MatrixPtr(new Matrix());
            (*khMinusL)[0] = 0;
            (*khMinusL)[4] = 1;
            (*khMinusL)[1] = 1;
            (*khMinusL)[5] = 0;
            (*khMinusL)[10] = -1;
            
            return khMinusL;
        }
        
        if (ambiguityCount() == 3)
        {
            // return -h -k -l
            MatrixPtr minusHminusKL = MatrixPtr(new Matrix());
            (*minusHminusKL)[0] = -1;
            (*minusHminusKL)[5] = -1;
            (*minusHminusKL)[10] = 1;
            
            return minusHminusKL;
        }
    }
    
    if (i == 2)
    {
        if (ambiguityCount() == 3 && customAmbiguity)
        {
            return customAmbiguity;
        }

        if (ambiguityCount() == 3)
        {
            if (spgNum == 149 || spgNum == 151 || spgNum == 153)
            {
                // return k h -l
                MatrixPtr khMinusL = MatrixPtr(new Matrix());
                (*khMinusL)[0] = 0;
                (*khMinusL)[4] = 1;
                (*khMinusL)[1] = 1;
                (*khMinusL)[5] = 0;
                (*khMinusL)[10] = -1;
                
                return khMinusL;
            }
            
            if (spgNum == 152 || spgNum == 152 || spgNum == 154)
            {
                // return -k -h -l
                MatrixPtr minusAllHKL = MatrixPtr(new Matrix());
                (*minusAllHKL)[0] = -1;
                (*minusAllHKL)[5] = -1;
                (*minusAllHKL)[10] = -1;
                
                return minusAllHKL;
            }
        }

        if (ambiguityCount() == 4)
        {
            // return k h -l
            MatrixPtr khMinusL = MatrixPtr(new Matrix());
            (*khMinusL)[0] = 0;
            (*khMinusL)[4] = 1;
            (*khMinusL)[1] = 1;
            (*khMinusL)[5] = 0;
            (*khMinusL)[10] = -1;
            
            return khMinusL;
        }
    }
    
    if (i == 3)
    {
        if (ambiguityCount() == 4 && customAmbiguity)
        {
            return customAmbiguity;
        }

        if (ambiguityCount() == 4)
        {
            // return -k -h -l
            MatrixPtr minusAllHKL = MatrixPtr(new Matrix());
            (*minusAllHKL)[0] = -1;
            (*minusAllHKL)[5] = -1;
            (*minusAllHKL)[10] = -1;
            
            return minusAllHKL;
        }

    }
    
    if (i == 4)
    {
        if (customAmbiguity)
        {
            return customAmbiguity;
        }
    }
    
    return MatrixPtr(new Matrix());
}


int Reflection::ambiguityCount()
{
    int basicNum = 1;
    
    if (spgNum >= 75 && spgNum <= 80)
    {
        basicNum = 2;
    }
    
    if (spgNum >= 195 && spgNum <= 199)
        basicNum = 2;
    
    if (spgNum >= 168 && spgNum <= 173)
        basicNum = 2;
    
    if (spgNum == 146)
        basicNum = 2;
    
    if (spgNum >= 149 && spgNum <= 154)
        basicNum = 3;
    
    if (spgNum >= 143 && spgNum <= 145)
        basicNum = 4;
    
    return basicNum + (customAmbiguity != MatrixPtr());
}


void Reflection::setSpaceGroup(int spaceGroupNum)
{
    if (hasSetup)
        return;
    
    setupMutex.lock();
    
    if (hasSetup)
    {
        setupMutex.unlock();
        return;
    }

    ccp4_space_group = ccp4spg_load_by_ccp4_num(spaceGroupNum);
    spgNum = spaceGroupNum;
    
    if (hasSetup)
        return;
    
    std::vector<double> ambiguityDouble = FileParser::getKey("CUSTOM_AMBIGUITY", std::vector<double>());
    
    if (ambiguityDouble.size() >= 9)
    {
        customAmbiguity = MatrixPtr(new Matrix(&ambiguityDouble[0]));
    }
    
    int totalAmbiguities = ambiguityCount();
    
    for (int i = 0; i < totalAmbiguities; i++)
    {
        flipMatrices.push_back(matrixForAmbiguity(i));
    }
    
    hasSetup = true;
    
    setupMutex.unlock();
}

int Reflection::reflectionIdForCoordinates(int h, int k, int l)
{
    int index = (h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (k + OFFSET) * MULTIPLIER + (l + OFFSET);
    
    return index;
}

int roundToInt(double value)
{
    if (value < 0)
    {
        value -= 0.1;
        return (int)value;
    }
    else
    {
        value += 0.1;
        return (int)value;
    }
}

void Reflection::generateReflectionIds()
{
    if (millerCount() == 0)
    {
        std::cout << "Warning! Miller count is 0" << std::endl;
    }
    
    int h = miller(0)->getH();
    int k = miller(0)->getK();
    int l = miller(0)->getL();
    
    vec miller = new_vector(h, k, l);
    
    for (int i = 0; i < ambiguityCount(); i++)
    {
        MatrixPtr ambiguityMat = matrixForAmbiguity(i);
        ambiguityMat->multiplyVector(&miller);
        int _h = roundToInt(miller.h);
        int _k = roundToInt(miller.k);
        int _l = roundToInt(miller.l);
        int h2, k2, l2;
        
        CSym::ccp4spg_put_in_asu(ccp4_space_group, _h, _k, _l, &h2, &k2, &l2);
        
        int newId = reflectionIdForCoordinates(h2, k2, l2);
        
        reflectionIds.push_back(newId);
    }
}

void Reflection::setUnitCellDouble(double *theUnitCell)
{
    unitCellMatrix = Matrix::matrixFromUnitCell(theUnitCell);
    
    setupUnitCell = true;
}


void Reflection::setUnitCell(float *theUnitCell)
{
    double params[6];
    params[0] = theUnitCell[0];
    params[1] = theUnitCell[1];
    params[2] = theUnitCell[2];
    params[3] = theUnitCell[3];
    params[4] = theUnitCell[4];
    params[5] = theUnitCell[5];
    
    unitCellMatrix = Matrix::matrixFromUnitCell(params);
    
    setupUnitCell = true;
}

Reflection::Reflection(float *unitCell, CSym::CCP4SPG *spg)
{
    if (spg != NULL)
        setSpaceGroup(spg->spg_num);

    if (unitCell != NULL)
    {
        setUnitCell(unitCell);
    }
    
    // TODO Auto-generated constructor stub
    
    resolution = 0;
    activeAmbiguity = 0;
    
    millerMutex = MutexPtr(new std::mutex());
    rejectSigma = FileParser::getKey("OUTLIER_REJECTION_SIGMA", OUTLIER_REJECTION_SIGMA);
    shouldReject = FileParser::getKey("OUTLIER_REJECTION", true);
}

MillerPtr Reflection::acceptedMiller(int num)
{
    int accepted = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (miller(i)->accepted() && accepted == num)
            return miller(i);
        
        if (miller(i)->accepted())
            accepted++;
    }
    
    return MillerPtr();
}


MillerPtr Reflection::miller(int i)
{
    return millers[i];
}

void Reflection::addMiller(MillerPtr miller)
{
    miller->setResolution(resolution);
    millers.push_back(miller);
    
    if (reflectionIds.size() == 0)
    {
        generateReflectionIds();
    }
}

void Reflection::addMillerCarefully(MillerPtr miller)
{
    millerMutex->lock();
    
    addMiller(miller);
    
    millerMutex->unlock();
}

bool Reflection::betweenResolutions(double lowAngstroms, double highAngstroms)
{
    double minD, maxD = 0;
    StatisticsManager::convertResolutions(lowAngstroms,
                                          highAngstroms, &minD, &maxD);
    
    if (resolution > maxD || resolution < minD)
        return false;
    
    return true;
}

int Reflection::millerCount()
{
    return (int)millers.size();
}

void Reflection::removeMiller(int index)
{
    millers.erase(millers.begin() + index);
}

ReflectionPtr Reflection::copy(bool copyMillers)
{
    ReflectionPtr newReflection = ReflectionPtr(new Reflection());
    
    newReflection->spgNum = spgNum;
    newReflection->activeAmbiguity = activeAmbiguity;
    newReflection->reflectionIds = reflectionIds;
    newReflection->resolution = resolution;
    
    for (int i = 0; i < millerCount(); i++)
    {
        MillerPtr newMiller;
        
        if (copyMillers)
        {
            newMiller = miller(i)->copy();
        }
        else
        {
            newMiller = miller(i);
        }
        
        newReflection->addMiller(newMiller);
    }

    return newReflection;
}

double Reflection::meanPartiality(bool withCutoff)
{
    int num = millerCount();
    double total_partiality = 0;
    int count = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if ((miller->accepted() && withCutoff) || !withCutoff)
        {
            total_partiality += miller->getPartiality();
            count++;
        }
    }
    
    total_partiality /= count;
    
    return total_partiality;
}

double Reflection::meanIntensity(bool withCutoff, int start, int end)
{
    if (end == 0)
        end = withCutoff ? acceptedCount() : millerCount();
    
    double total_intensity = 0;
    double total_weights = 0;
    
    for (int i = start; i < end; i++)
    {
        MillerPtr miller = withCutoff ? this->acceptedMiller(i) : this->miller(i);

        double weight = miller->getWeight(withCutoff);
        double intensity = miller->intensity(withCutoff);
        
        if (weight <= 0)
            continue;
        
        total_intensity += intensity * weight;
        total_weights += weight;
    }
    
    total_intensity /= total_weights;
    
    return total_intensity;
}

double Reflection::meanIntensityWithExclusion(std::string *filename, int start, int end)
{
    if (filename == NULL)
        return meanIntensity(true, start, end);
    
    if (end == 0)
        end = acceptedCount();
    double total_intensity = 0;
    double weight = 0;
    
    for (int i = start; i < end; i++)
    {
        MillerPtr miller = this->acceptedMiller(i);
        
        total_intensity += miller->intensity() * miller->getPartiality();
        weight += miller->getPartiality();
    }
    
    total_intensity /= weight;
    
    return total_intensity;
}

double Reflection::mergeSigma()
{
    double mean = mergedIntensity(WeightTypePartialitySigma);
    
    double weights = 0;
    double sumSquares = 0;
    int count = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (!millers[i]->accepted())
            continue;
        
        double weight = millers[i]->getWeight();
        count++;
        
        sumSquares += pow(millers[i]->intensity() - mean, 2) * weight;
        weights += weight;
    }
    
    double stdev = sqrt(sumSquares / weights);
    
    double error = stdev / sqrt(count - 1);
    
    if (count == 1)
        error = -1;
    
    if (count == 0)
        return nan(" ");
    
    
    return error;
}

double Reflection::meanSigma(bool friedel)
{
    int num = (int)millerCount();
    int count = 0;
    
    double total_sigi = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr aMiller = miller(i);
        
        if (aMiller->accepted())
        {
            total_sigi += aMiller->getSigma();
            count++;
        }
    }
    
    total_sigi /= count;
    
    if (total_sigi == 0)
        return nan(" ");
    
    return total_sigi;
}

double Reflection::meanWeight(bool withCutoff)
{
    int num = (int)millerCount();
    int count = 0;
    
    double total_weight = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if ((miller->accepted() && withCutoff) || !withCutoff)
        {
            total_weight += miller->getWeight(withCutoff);
            count++;
        }
    }
    
    total_weight /= count;
    
    return total_weight;
}

double Reflection::meanSigma()
{
    int num = (int)millers.size();
    int count = 0;
    
    double total_sigi = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if (miller->accepted())
        {
            total_sigi += miller->getSigma();
            count++;
        }
    }
    
    total_sigi /= count;
    
    return total_sigi;
}

void Reflection::calculateResolution(MtzManager *mtz)
{
    int h = millers[0]->getH();
    int k = millers[0]->getK();
    int l = millers[0]->getL();
    
    vec hkl = new_vector(h, k, l);
    unitCellMatrix->multiplyVector(&hkl);
    resolution = length_of_vector(hkl);
    
    for (int i = 0; i < millerCount(); i++)
    {
        miller(i)->setResolution(resolution);
    }
}


void Reflection::incrementAmbiguity()
{
    int count = ambiguityCount();
    int newActive = activeAmbiguity + 1;
    activeAmbiguity = (newActive % count);
}

void Reflection::setFlipAsActiveAmbiguity()
{
    setFlip(activeAmbiguity);
}

void Reflection::resetFlip()
{
    for (int j = 0; j < millerCount(); j++)
    {
        miller(j)->setFlipMatrix(0);
    }
}

void Reflection::setFlip(int i)
{
    for (int j = 0; j < millerCount(); j++)
    {
        miller(j)->setFlipMatrix(i);
    }
}

MatrixPtr Reflection::getFlipMatrix(int i)
{
    if (flipMatrices.size() == 0)
        return Matrix::getIdentityPtr();
    
    return flipMatrices[i];
}

void Reflection::reflectionDescription()
{
    int acceptedCount = 0;
    
    std::ostringstream logged;
    
    for (int i = 0; i < millerCount(); i++)
    {
        MillerPtr miller = this->miller(i);
        logged << miller->getH() << "\t" << miller->getK() << "\t" << miller->getL() << "\t"
        << miller->getRawIntensity() << "\t" << miller->getPartiality()
        << "\t" << miller->getSigma() << "\t" << miller->getMtzParent()->getFilename()
        << std::endl;
        if (miller->accepted())
            acceptedCount++;
    }
    logged << std::endl;
    
    Logger::mainLogger->addStream(&logged);
}

void Reflection::clearMillers()
{
    millers.clear();
    vector<MillerPtr>().swap(millers);
    
}

Reflection::~Reflection()
{
    clearMillers();
}

int Reflection::indexForReflection(int h, int k, int l, CSym::CCP4SPG *spgroup,
                               bool inverted)
{
    int _h, _k, _l;
    
    ccp4spg_put_in_asu(spgroup, h, k, l, &_h, &_k, &_l);
    
    int multiplier = MULTIPLIER;
    int offset = OFFSET;
    
    int index = (_h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (_k + OFFSET) * MULTIPLIER + (_l + OFFSET);
    if (inverted)
        index = (_h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
        + (_l + OFFSET) * MULTIPLIER + (_k + OFFSET);
    
    if (spgroup->spg_num == 197)
    {
        if (inverted == 0)
            index = (_h + offset) * pow((double) multiplier, (int) 2)
            + (_k + offset) * multiplier + (_l + offset);
        else
            index = (_h + offset) * pow((double) multiplier, (int) 2)
            + (_l + offset) * multiplier + (_k + offset);
    }
    else if (spgroup->spg_num == 146)
    {
        if (inverted == 0)
            index = (_h + offset) * pow((double) multiplier, (int) 2)
            + (_k + offset) * multiplier + (_l + offset);
        else
            index = (_k + offset) * pow((double) multiplier, (int) 2)
            + (_h + offset) * multiplier + ((-_l) + offset);
    }
    
    return index;
}

double Reflection::mergedIntensity(WeightType weighting)
{
    double sum_intensities = 0;
    double sum_weights = 0;
    
    for (int i = 0; i < millers.size(); i++)
    {
        if (!miller(i)->accepted())
            continue;
        
        double weight = miller(i)->getWeight(weighting);
        
        sum_intensities += millers[i]->intensity() * weight;
        sum_weights += weight;
    }
    
    double mean_intensity = sum_intensities / sum_weights;
    
    return mean_intensity;
}

double Reflection::standardDeviation(WeightType weighting)
{
    double mean_intensity = mergedIntensity(weighting);
    
    // we do not weight standard deviation
    double squares = 0;
    int num = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (!miller(i)->accepted())
            continue;
        
        squares += pow(miller(i)->intensity() - mean_intensity, 2);
        num++;
    }
    
    double stdev = sqrt(squares / (double) num);
    
    if (num == 1)
        return 100;
    
    return stdev;
}

double Reflection::rMergeContribution(double *numerator, double *denominator)
{
    if (acceptedCount() < 2)
        return 0;

    double mean_intensity = mergedIntensity(WeightTypePartialitySigma);
    
    double littleNum = 0;
    double littleDenom = 0;
    int count = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (!miller(i)->accepted())
            continue;
        
        double intensity = miller(i)->intensity();
        double weight = miller(i)->getWeight();
        
        if (intensity < 0 || mean_intensity < 0)
            continue;
        
        count++;
        littleNum += weight * fabs(intensity - mean_intensity);
        littleDenom += weight * fabs(intensity);
    }
    
    *numerator += littleNum;
    *denominator += littleDenom;
    
    return littleNum / littleDenom;
}

void Reflection::medianMerge(double *intensity, double *sigma, int *rejected, signed char friedel)
{
    std::vector<double> intensities, weights;
    
    for (int i = 0; i < liteMillers.size(); i++)
    {
        if (friedel != -1)
        {
            if (liteMillers[i].friedel != friedel)
            {
                continue;
            }
        }
        
        intensities.push_back(liteMillers[i].intensity);
    }
    
    *intensity = median(&intensities);
    *sigma = standard_deviation(&intensities);
    
    if (rejected != NULL)
    {
        *rejected = 0;        
    }
    
    /*
    std::ostringstream logged;
    logged << "Int " << *intensity << " from " << intensities.size() << " reflections for (" << miller(0)->getH() << ", " << miller(0)->getK() << ", " << miller(0)->getL() << ")" << std::endl;
    Logger::log(logged);*/
}


void Reflection::liteMerge(double *intensity, double *sigma, int *rejected, signed char friedel)
{
    std::vector<double> intensities, weights;
    
    if (rejected != NULL)
    {
        *rejected = 0;
    }
    
    for (int i = 0; i < liteMillers.size(); i++)
    {
        if (friedel != -1)
        {
            if (liteMillers[i].friedel != friedel)
            {
                continue;
            }
        }
        
        intensities.push_back(liteMillers[i].intensity);
        weights.push_back(liteMillers[i].weight);
    }
    
    double mean = weighted_mean(&intensities, &weights);
    double stdev = standard_deviation(&intensities, NULL, mean);
    
    bool shouldRejectLocal = (rejected != NULL) * shouldReject;
    
    if (shouldRejectLocal && liteMillers.size() >= MIN_MILLER_COUNT)
    {
        intensities.clear();
        weights.clear();
        
        int minIntensity = mean - stdev * rejectSigma;
        int maxIntensity = mean + stdev * rejectSigma;
        
        for (int i = 0; i < liteMillers.size(); i++)
        {
            if (friedel != -1)
            {
                if (liteMillers[i].friedel != friedel)
                {
                    continue;
                }
            }

            double testIntensity = liteMillers[i].intensity;
            
            if (testIntensity > maxIntensity || testIntensity < minIntensity)
            {
                (*rejected)++;
                continue;
            }
            
            intensities.push_back(testIntensity);
            weights.push_back(liteMillers[i].weight);
        }

        mean = weighted_mean(&intensities, &weights);
        stdev = standard_deviation(&intensities, NULL, mean);
    }
    
    stdev /= sqrt(intensities.size());
    
    if (intensities.size() == 1)
    {
        stdev = -1;
    }
    
    *intensity = mean;
    *sigma = stdev;
}

void Reflection::clearLiteMillers()
{
    liteMillers.clear();
    std::vector<LiteMiller>().swap(liteMillers);
}

void Reflection::merge(WeightType weighting, double *intensity, double *sigma,
                   bool calculateRejections)
{
    std::ostringstream logged;
    
    logged << "Rejection info:" << std::endl;
    
    int isSpecial = (getReflId() == reflectionIdForCoordinates(5, 10, 7));
    isSpecial = 0;
    
    if (calculateRejections)
    {
        for (int i = 0; i < millerCount(); i++)
        {
            if (isSpecial)
            {
                std::ostringstream aLog;
                aLog << miller(i)->getH() << " " << miller(i)->getK() << " " << miller(i)->getL() << " " << miller(i)->getMtzParent()->getFilename() << " " << miller(i)->intensity() << " " << miller(i)->getRejectionFlags() << " " << miller(i)->getPartiality() << std::endl;
                Logger::log(aLog);
            }
            
     //       miller(i)->setRejected(false);
            
            miller(i)->setRejected(RejectReasonMerge, false);
        }
    }
    
    double mean_intensity = mergedIntensity(weighting);
    double stdev = standardDeviation(weighting);
    
    if (acceptedCount() < MIN_MILLER_COUNT || !REJECTING_MILLERS
        || !calculateRejections)
    {
        *intensity = mean_intensity;
        *sigma = stdev;
        return;
    }
    
    shouldReject = FileParser::getKey("OUTLIER_REJECTION", true);
    
    if (shouldReject)
    {
        double error = stdev * rejectSigma;
        
        logged << "Std error: " << error << std::endl;
        
        int rejectedCount = 0;
        
        double minIntensity = mean_intensity - error;
        double maxIntensity = mean_intensity + error;
        
        for (int i = 0; i < millerCount(); i++)
        {
            if (!miller(i)->accepted())
                continue;
            
            if (miller(i)->intensity() > minIntensity
                && miller(i)->intensity() < maxIntensity)
                continue;
            
            miller(i)->setRejected(RejectReasonMerge, true);
            rejectedCount++;
        }
        
        logged << "Rejected " << rejectedCount << " reflections between " << minIntensity << " and " << maxIntensity << std::endl;
    }
    
    mean_intensity = meanIntensity();
    double newSigma = stdev;
    
    *intensity = mean_intensity;
    *sigma = newSigma;
    
    Logger::mainLogger->addStream(&logged, isSpecial ? LogLevelNormal : LogLevelDebug);
}

int Reflection::rejectCount()
{
    int count = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (miller(i)->isRejected())
        {
            count++;
        }
    }
    
    return count;
}

int Reflection::acceptedCount()
{
    int count = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (miller(i)->accepted())
        {
            count++;
        }
    }
    
    return count;
}

bool Reflection::anyAccepted()
{
    for (int i = 0; i < millers.size(); i++)
    {
        if (miller(i)->accepted())
            return true;
    }
    
    return false;
}

double Reflection::observedPartiality(MtzManager *reference, Miller *miller)
{
    ReflectionPtr refReflection;
    double reflId = getReflId();
    reference->findReflectionWithId(reflId, &refReflection);
    
    if (refReflection)
        return miller->observedPartiality(refReflection->meanIntensity());
    
    return nan(" ");
}

void Reflection::printDescription()
{
    std::cout << "Mean intensity " << this->meanIntensity() << std::endl;
    std::cout << "Miller corrected intensities: ";
    
    for (int i = 0; i < millerCount(); i++)
    {
        std::cout << miller(i)->intensity() << ", ";
    }
    
    std::cout << std::endl;
}

void Reflection::detailedDescription()
{
    double meanIntensity = this->meanIntensity();
    
    std::cout << "Mean intensity " << meanIntensity << std::endl;
    std::cout << "Resolution " << 1 / this->getResolution() << " Å" << std::endl;
    
    for (int i = 0; i < millerCount(); i++)
    {
        double rawIntensity = miller(i)->getRawIntensity();
        double fraction = rawIntensity / meanIntensity;
        bool friedel = false; int isym = 0;
        miller(i)->positiveFriedel(&friedel, &isym);
        
        std::cout << miller(i)->getPartiality() << "\t" << fraction << "\t" << rawIntensity << "\t" << isym << std::endl;
    }
    
    std::cout << std::endl;
}

int Reflection::checkSpotOverlaps(std::vector<SpotPtr> *spots, bool actuallyDelete)
{
    int count = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (miller(i)->isOverlappedWithSpots(spots, actuallyDelete))
        {
            count++;
        }
    }
    
    return count;
}

int Reflection::checkOverlaps()
{
    int count = 0;
    
    for (int i = 0; i < millerCount(); i++)
    {
        if (miller(i)->isOverlapped())
        {
            count++;
            removeMiller(i);
        }
    }
    
    return count;
}

void Reflection::addLiteMiller(MillerPtr miller)
{
    double intensity = miller->intensity();
    double weight = miller->getWeight();
    /*
    double reflId = reflectionIdForCoordinates(5, 10, 7);
    
    if (reflId == getReflId())
    {
   //     std::ostringstream logged;
   //     logged << "5 22 7: " << intensity << ", " << weight << ", " << miller->getScale() << std::endl;
   //     Logger::log(logged);
    }
    */
    LiteMiller liteMiller;
    liteMiller.intensity = intensity;
    miller->positiveFriedel(&(liteMiller.friedel));
    liteMiller.weight = weight;
    
    millerMutex->lock();
    
    liteMillers.push_back(liteMiller);
    
    millerMutex->unlock();
}

