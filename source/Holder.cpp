/*
 * Holder.cpp
 *
 *  Created on: 25 Oct 2014
 *      Author: helenginn
 */

#include "Holder.h"

#include <vector>
#include <cmath>
#include "headers/csymlib.h"

#include "FileParser.h"
#include "StatisticsManager.h"

Holder::Holder()
{
    // TODO Auto-generated constructor stub
    
    refIntensity = 0;
    refSigma = 0;
    refl_id = 0;
    inv_refl_id = 0;
    resolution = 0;
}

MillerPtr Holder::miller(int i)
{
    return millers[i];
}

void Holder::addMiller(MillerPtr miller)
{
    miller->setResolution(resolution);
    millers.push_back(miller);
}

bool Holder::betweenResolutions(double lowAngstroms, double highAngstroms)
{
    double minD, maxD = 0;
    StatisticsManager::convertResolutions(lowAngstroms,
                                          highAngstroms, &minD, &maxD);
    
    if (resolution > maxD || resolution < minD)
        return false;
    
    return true;
}

int Holder::millerCount()
{
    return (int)millers.size();
}

void Holder::removeMiller(int index)
{
    millers.erase(millers.begin() + index);
}

Holder *Holder::copy()
{
    return copy(false);
}

Holder *Holder::copy(bool copyMillers)
{
    Holder *newHolder = new Holder();
    
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
        
        newHolder->addMiller(newMiller);
    }
    
    newHolder->inverse = inverse;
    newHolder->refIntensity = refIntensity;
    newHolder->refSigma = refSigma;
    newHolder->refl_id = refl_id;
    newHolder->inv_refl_id = inv_refl_id;
    newHolder->resolution = resolution;
    
    return newHolder;
}

double Holder::meanPartiality()
{
    int num = (int)millers.size();
    double total_partiality = 0;
    int count = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if (miller->accepted())
        {
            total_partiality += miller->getPartiality();
            count++;
        }
    }
    
    total_partiality /= count;
    
    return total_partiality;
}

double Holder::meanIntensity()
{
    int num = (int)millers.size();
    double total_intensity = 0;
    double weight = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = this->miller(i);
        
        if (miller->accepted())
        {
            total_intensity += miller->intensity() * miller->getPartiality();
            weight += miller->getPartiality();
        }
    }
    
    total_intensity /= weight;
    
    return total_intensity;
}

double Holder::meanIntensityWithExclusion(string *filename)
{
    if (filename == NULL)
        return meanIntensity();
    
    int num = (int)millers.size();
    double total_intensity = 0;
    double weight = 0;
    int acceptedCount = 0;
    bool noFilenames = true;
    
    for (int i = 0; i < num; i++)
    {
        if (miller(i)->accepted())
            acceptedCount++;
        
        if (miller(i)->getFilename().substr(0, 3) == "img")
            noFilenames = false;
        
        if (acceptedCount >= 2)
            break;
    }
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = this->miller(i);
        
        if (miller->accepted())
        {
            if (acceptedCount > 2 && miller->getFilename() == *filename)
                continue;
            
            total_intensity += miller->intensity() * miller->getPartiality();
            weight += miller->getPartiality();
        }
    }
    
    total_intensity /= weight;
    
    return total_intensity;
}

double Holder::meanSigma(bool friedel)
{
    int num = (int)millers.size();
    int count = 0;
    
    double total_sigi = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if (miller->accepted())
        {
            if (miller->positiveFriedel() != friedel)
                continue;
            
            total_sigi += miller->getSigma();
            count++;
        }
    }
    
    total_sigi /= count;
    
    return total_sigi;
    
    return total_sigi;
}

double Holder::meanWeight()
{
    int num = (int)millers.size();
    int count = 0;
    
    double total_weight = 0;
    
    for (int i = 0; i < num; i++)
    {
        MillerPtr miller = millers[i];
        
        if (miller->accepted())
        {
            total_weight += miller->getWeight();
            count++;
        }
    }
    
    total_weight /= count;
    
    return total_weight;
}

double Holder::meanSigma()
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

void Holder::calculateResolution(MtzManager *mtz)
{
    double a, b, c, alpha, beta, gamma;
    
    mtz->getUnitCell(&a, &b, &c, &alpha, &beta, &gamma);
    
    Matrix mat = Matrix::matrixFromUnitCell(a, b, c, alpha, beta, gamma);
    
    int h = millers[0]->h;
    int k = millers[0]->k;
    int l = millers[0]->l;
    
    vec hkl = new_vector(h, k, l);
    
    mat.multiplyVector(&hkl);
    
    resolution = length_of_vector(hkl);
    
    for (int i = 0; i < millerCount(); i++)
    {
        miller(i)->setResolution(resolution);
    }
}

void Holder::flipMiller(MillerPtr miller, int spg_num)
{
    if (spg_num == 197)
    {
        int tmp2 = miller->h;
        miller->h = miller->k;
        miller->k = tmp2;
        miller->l = - miller->l;
    }
    else if (spg_num == 146)
    {
        int tmp2 = miller->k;
        miller->k = miller->h;
        miller->h = tmp2;
        miller->l = 0 - miller->l;
    }
    else if (spg_num == 4 || spg_num == 21)
    {
        int tmp2 = miller->h;
        miller->h = miller->l;
        miller->l = tmp2;
    }
    else
    {
//        std::cout << "Warning! Unable to flip miller" << std::endl;
    }
}

void Holder::flip()
{
    flip(MtzManager::getReferenceManager());
    
    inverse = !inverse;
}

void Holder::flip(MtzManager *mtz)
{
    CCP4SPG *spaceGroup = mtz->getLowGroup();
    
    int tmp = getInvReflId();
    setInvReflId(refl_id);
    setReflId(tmp);
    
    for (int j = 0; j < millerCount(); j++)
    {
        flipMiller(miller(j), spaceGroup->spg_num);
    }
}

void Holder::holderDescription()
{
    int acceptedCount = 0;
    
    std::ostringstream logged;
    
    for (int i = 0; i < millerCount(); i++)
    {
        MillerPtr miller = this->miller(i);
        logged << miller->h << "\t" << miller->k << "\t" << miller->l << "\t"
        << miller->getRawIntensity() << "\t" << miller->getPartiality()
        << "\t" << miller->getSigma() << "\t" << miller->getFilename()
        << std::endl;
        if (miller->accepted())
            acceptedCount++;
    }
    logged << std::endl;
    
    Logger::mainLogger->addStream(&logged);
}

void Holder::clearMillers()
{
    millers.clear();
    vector<MillerPtr>().swap(millers);
    
}

Holder::~Holder()
{
    millers.clear();
    vector<MillerPtr>().swap(millers);
}

int Holder::indexForReflection(int h, int k, int l, CSym::CCP4SPG *spgroup,
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

double Holder::mergedIntensity(WeightType weighting)
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

double Holder::standardDeviation(WeightType weighting)
{
    /*	if (acceptedCount() < MIN_MILLER_COUNT)
     {
     return meanSigma() / meanPartiality();
     }
     */
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

void Holder::merge(WeightType weighting, double *intensity, double *sigma,
                   bool calculateRejections)
{
    ostringstream logged;
    
    logged << "Rejection info:" << std::endl;
    
    if (calculateRejections)
    {
        for (int i = 0; i < millerCount(); i++)
        {
            miller(i)->setRejected(false);
            //	miller(i)->setRejected("partiality", false);
        }
    }
    
    double mean_intensity = mergedIntensity(weighting);
    double stdev = standardDeviation(weighting);
    
    if (acceptedCount() < MIN_MILLER_COUNT || !REJECTING_MILLERS
        || !calculateRejections)
    {
        *intensity = mean_intensity;
        *sigma = stdev; // meanSigma() / meanPartiality();
        return;
    }
    
    double rejectSigma = FileParser::getKey("OUTLIER_REJECTION_SIGMA",
                                           OUTLIER_REJECTION_SIGMA);
    
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
        
        miller(i)->setRejected(true);
        rejectedCount++;
    }
    
    logged << "Rejected " << rejectedCount << " reflections between " << minIntensity << " and " << maxIntensity << std::endl;
    
    mean_intensity = meanIntensity();
    double newSigma = stdev;
    
    *intensity = mean_intensity;
    *sigma = newSigma;
    
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
}

int Holder::rejectCount()
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

int Holder::acceptedCount()
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

bool Holder::anyAccepted()
{
    for (int i = 0; i < millers.size(); i++)
    {
        if (miller(i)->accepted())
            return true;
    }
    
    return false;
}

double Holder::observedPartiality(MtzManager *reference, Miller *miller)
{
    Holder *refHolder;
    double reflId = refl_id;
    reference->findHolderWithId(reflId, &refHolder);
    
    if (refHolder != NULL)
        return miller->observedPartiality(refHolder->meanIntensity());
    
    return nan(" ");
}

void Holder::printDescription()
{
    std::cout << "Mean intensity " << this->meanIntensity() << std::endl;
    std::cout << "Miller corrected intensities: ";
    
    for (int i = 0; i < millerCount(); i++)
    {
        std::cout << miller(i)->intensity() << ", ";
    }
    
    std::cout << std::endl;
}

int Holder::checkOverlaps()
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
