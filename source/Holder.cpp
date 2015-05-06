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
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller.h>

#include "FileParser.h"
#include "StatisticsManager.h"

using cctbx::sgtbx::space_group_type;
using cctbx::miller::sym_equiv_indices;
using cctbx::miller::asym_index;
using cctbx::sgtbx::reciprocal_space::asu;

asu Holder::asymmetricUnit = asu();
bool Holder::hasSetup = false;
bool Holder::setupUnitCell = false;
space_group_type Holder::spgType = cctbx::sgtbx::space_group_type();
cctbx::sgtbx::space_group Holder::spaceGroup = cctbx::sgtbx::space_group();
cctbx::uctbx::unit_cell Holder::unitCell;

MatrixPtr Holder::matrixForAmbiguity(int i)
{
    if (i >= ambiguityCount())
        std::cout << "Ambiguity issue!" << std::endl;
    
    if (i == 0)
    {
        MatrixPtr identity = MatrixPtr(new Matrix());
        
        return identity;
    }
    
    if (i == 1)
    {
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
    
    return MatrixPtr(new Matrix());
}

int Holder::ambiguityCount()
{
    if (spgNum >= 195 && spgNum <= 199)
        return 2;
    
    if (spgNum >= 168 && spgNum <= 173)
        return 2;
    
    if (spgNum == 146)
        return 2;
    
    if (spgNum >= 149 && spgNum <= 154)
        return 3;
    
    if (spgNum >= 143 && spgNum <= 145)
        return 4;
    
    return 1;
}

void Holder::setSpaceGroup(CSym::CCP4SPG *ccp4spg, cctbx::sgtbx::space_group_type newSpgType, asu newAsymmetricUnit)
{
    if (ccp4spg == NULL)
        return;
    
    char *hallSymbol = ccp4spg_symbol_Hall(ccp4spg);
    spgNum = ccp4spg->spg_num;
    
    if (hasSetup)
        return;
    
    spaceGroup = space_group(hallSymbol);
    spgType = newSpgType;
    asymmetricUnit = newAsymmetricUnit;
    
    hasSetup = true;
}

void Holder::setSpaceGroup(int spaceGroupNum)
{
    if (hasSetup)
        return;
    
    spgNum = spaceGroupNum;
    space_group_symbols spaceGroupSymbol = space_group_symbols(spaceGroupNum);
    std::string hallSymbol = spaceGroupSymbol.hall();
    
    spaceGroup = space_group(hallSymbol);
    spgType = cctbx::sgtbx::space_group_type(spaceGroup);
    asymmetricUnit = asu(spgType);
    
    hasSetup = true;
}

int Holder::reflectionIdForMiller(cctbx::miller::index<> cctbxMiller)
{
    int h = cctbxMiller[0];
    int k = cctbxMiller[1];
    int l = cctbxMiller[2];
    
    int index = (h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (k + OFFSET) * MULTIPLIER + (l + OFFSET);
    
    return index;
}

void Holder::generateReflectionIds()
{
    if (millerCount() == 0)
    {
        std::cout << "Warning! Miller count is 0" << std::endl;
    }
    
    int h = miller(0)->getH();
    int k = miller(0)->getK();
    int l = miller(0)->getL();
    
    cctbx::miller::index<> cctbxMiller = cctbx::miller::index<>(h, k, l);
        for (int i = 0; i < ambiguityCount(); i++)
    {
        MatrixPtr ambiguityMat = matrixForAmbiguity(i);
        cctbx::miller::index<> cctbxTwinnedMiller = ambiguityMat->multiplyIndex(&cctbxMiller);
        
        asym_index asymmetricMiller = asym_index(spaceGroup, asymmetricUnit, cctbxTwinnedMiller);
        
        
    //    sym_equiv_indices equivMaker = sym_equiv_indices(spaceGroup, cctbxTwinnedMiller);
    //    cctbx::miller::index<> asymmetricMiller = equivMaker(0).h();
       
        int newId = reflectionIdForMiller(asymmetricMiller.h());
        
        reflectionIds.push_back(newId);
    }
}

void Holder::setUnitCell(float *theUnitCell)
{
    scitbx::af::double6 params;
    params[0] = theUnitCell[0];
    params[1] = theUnitCell[1];
    params[2] = theUnitCell[2];
    params[3] = theUnitCell[3];
    params[4] = theUnitCell[4];
    params[5] = theUnitCell[5];
    
    unitCell = cctbx::uctbx::unit_cell(params);
    
    setupUnitCell = true;
}

Holder::Holder(float *unitCell, CSym::CCP4SPG *spg)
{
    spgNum = 0;
    
    if (spg != NULL)
        setSpaceGroup(spg->spg_num);

    if (unitCell != NULL)
    {
        setUnitCell(unitCell);
    }
    
    // TODO Auto-generated constructor stub
    
    refIntensity = 0;
    refSigma = 0;
    refl_id = 0;
    inv_refl_id = 0;
    resolution = 0;
    activeAmbiguity = 0;
}

MillerPtr Holder::miller(int i)
{
    return millers[i];
}

void Holder::addMiller(MillerPtr miller)
{
    miller->setResolution(resolution);
    millers.push_back(miller);
    
    if (reflectionIds.size() == 0)
    {
        generateReflectionIds();
    }
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
    
    newHolder->spgNum = spgNum;
    newHolder->activeAmbiguity = activeAmbiguity;
    newHolder->reflectionIds = reflectionIds;
    newHolder->refIntensity = refIntensity;
    newHolder->refSigma = refSigma;
    newHolder->refl_id = refl_id;
    newHolder->inv_refl_id = inv_refl_id;
    newHolder->resolution = resolution;
    
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

double Holder::mergeSigma()
{
    double mean = mergedIntensity(WeightTypePartialitySigma);
    
    double weights = 0;
    double sumSquares = 0;
    double count = 0;
    
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
    
    double error = stdev / sqrt(count);
    
    return error;
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
            bool theFriedel = false;
            bool shouldUse = miller->positiveFriedel(&theFriedel);
            
            if (!shouldUse || theFriedel != friedel)
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
  /*  double a, b, c, alpha, beta, gamma;
    
    mtz->getUnitCell(&a, &b, &c, &alpha, &beta, &gamma);
    
    Matrix mat = Matrix::matrixFromUnitCell(a, b, c, alpha, beta, gamma);*/
    
    int h = millers[0]->getH();
    int k = millers[0]->getK();
    int l = millers[0]->getL();
    
    cctbx::miller::index<> anyMiller = cctbx::miller::index<>(h, k, l);
    
    resolution = unitCell.two_stol(anyMiller);
    
 /*   asym_index asymmetricMiller = asym_index(spaceGroup, asymmetricUnit, anyMiller);
    cctbx::miller::index<> asymmetricIndex = asymmetricMiller.h();
    
    vec hkl = new_vector(asymmetricIndex.as_tiny()[0], asymmetricIndex.as_tiny()[1], asymmetricIndex.as_tiny()[2]);
    
    mat.multiplyVector(&hkl);
    
    resolution = length_of_vector(hkl);*/
    
    for (int i = 0; i < millerCount(); i++)
    {
        miller(i)->setResolution(resolution);
    }
}


void Holder::incrementAmbiguity()
{
    int count = ambiguityCount();
    int newActive = activeAmbiguity + 1;
    activeAmbiguity = (newActive % count);
}

void Holder::setFlipAsActiveAmbiguity()
{
    setFlip(activeAmbiguity);
}

void Holder::resetFlip()
{
    for (int j = 0; j < millerCount(); j++)
    {
        miller(j)->setFlipMatrix(MatrixPtr(new Matrix()));
    }
}

void Holder::setFlip(int i)
{
    for (int j = 0; j < millerCount(); j++)
    {
        MatrixPtr ambiguityMat = matrixForAmbiguity(i);
        
        miller(j)->setFlipMatrix(ambiguityMat);
    }
}


void Holder::holderDescription()
{
    int acceptedCount = 0;
    
    std::ostringstream logged;
    
    for (int i = 0; i < millerCount(); i++)
    {
        MillerPtr miller = this->miller(i);
        logged << miller->getH() << "\t" << miller->getK() << "\t" << miller->getL() << "\t"
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
    clearMillers();
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

double Holder::rMergeContribution(double *numerator, double *denominator)
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
        
        count++;
        littleNum += abs(miller(i)->intensity() - mean_intensity);
        littleDenom += abs(miller(i)->intensity());
    }
    
    *numerator += littleNum;
    *denominator += littleDenom;
    
    return littleNum / littleDenom;
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

void Holder::detailedDescription()
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

void Holder::setAdditionalWeight(double weight)
{
    for (int i = 0; i < millerCount(); i++)
    {
        miller(i)->setAdditionalWeight(weight);
    }
}
