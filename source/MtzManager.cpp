#include "MtzManager.h"

#include <string>
#include <iostream>
#include "headers/cmtzlib.h"
#include "headers/csymlib.h"
#include "headers/ccp4_spg.h"
#include "headers/ccp4_general.h"
#include "headers/ccp4_parser.h"
#include "Vector.h"
#include "FileReader.h"
#include <cerrno>
#include <fstream>
#include <iostream>
#include <boost/variant.hpp>
#include "StatisticsManager.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "definitions.h"
#include "FileParser.h"

using namespace std;
using namespace CMtz;
using namespace CSym;

MtzManager *MtzManager::referenceManager;
//ScoreType MtzManager::scoreType = ScoreTypeCorrelation;
double MtzManager::highRes = 3.5;
double MtzManager::lowRes = 0;

string MtzManager::describeScoreType()
{
    switch (scoreType)
    {
        case ScoreTypeCorrelation:
            return string("correl");
        case ScoreTypePartialityCorrelation:
            return string("part");
        case ScoreTypePartialityLeastSquares:
            return string("part");
        case ScoreTypeMinimizeRSplit:
            return string("rfactor");
        case ScoreTypeMinimizeRSplitLog:
            return string("logR");
        case ScoreTypeCorrelationLog:
            return string("logCC");
        default:
            return string("unknown");
    }
}

void MtzManager::makeSuperGaussianLookupTable(double exponent)
{
    if (optimisingExponent)
        return;
    
    if (exponent == lastExponent)
        return;
    
    superGaussianScale = pow((M_PI / 2), (2 / exponent - 1));
    
    superGaussianTable.clear();
    
    const double min = 0;
    const double max = MAX_SUPER_GAUSSIAN;
    const double step = SUPER_GAUSSIAN_STEP;
    int count = 0;
    
    for (double x = min; x < max; x += step)
    {
        double value = super_gaussian(x, 0, superGaussianScale, exponent);
        
        //      std::cout << x << "\t" << value << std::endl;
        
        superGaussianTable.push_back(value);
        
        count++;
    }
    
    lastExponent = exponent;
}

std::string MtzManager::filenameRoot()
{
    std::vector<std::string> components = FileReader::split(filename, '.');
    return components[0];
}

double MtzManager::extreme_index(MTZ *mtz, int max)
{
    double extreme = 0;
    MTZCOL *col_h = MtzColLookup(mtz, "H");
    MTZCOL *col_k = MtzColLookup(mtz, "K");
    MTZCOL *col_l = MtzColLookup(mtz, "L");
    
    if (max == 0)
    {
        extreme = col_h->min;
        extreme = (extreme > col_k->min) ? col_k->min : extreme;
        extreme = (extreme > col_l->min) ? col_l->min : extreme;
    }
    else
    {
        extreme = col_h->max;
        extreme = (extreme < col_k->max) ? col_k->max : extreme;
        extreme = (extreme < col_l->max) ? col_l->max : extreme;
    }
    
    return extreme;
}

void MtzManager::findMultiplier(MTZ *mtz, int *multiplier, int *offset)
{
    int hkl_min = extreme_index(mtz, 0);
    int hkl_max = extreme_index(mtz, 1);
    
    *offset = (hkl_min < 0) ? -hkl_min : 0;
    
    *multiplier = hkl_max - ((hkl_min < 0) ? hkl_min : 0);
    
    *multiplier = MULTIPLIER;
    *offset = OFFSET;
}

void MtzManager::hkls_for_reflection(MTZ *mtz, float *adata, int *h, int *k,
                                     int *l, int *multiplier, int *offset)
{
    if (multiplier == 0)
        findMultiplier(mtz, multiplier, offset);
    
    MTZCOL *col_h = MtzColLookup(mtz, "H");
    MTZCOL *col_k = MtzColLookup(mtz, "K");
    MTZCOL *col_l = MtzColLookup(mtz, "L");
    
    *h = adata[col_h->source - 1];
    *k = adata[col_k->source - 1];
    *l = adata[col_l->source - 1];
}

int MtzManager::index_for_reflection(int h, int k, int l, bool inverted)
{
    int _h = h;
    int _k = k;
    int _l = l;
    
    if (inverted)
    {
        if (low_group->spg_num == 197)
        {
            _h = k;
            _k = h;
            _l = -l;
        }
        
        if (low_group->spg_num == 146)
        {
            _h = k;
            _k = h;
            _l = -l;
        }
        
        if (low_group->spg_num == 4 || low_group->spg_num == 21)
        {
            _h = l;
            _k = k;
            _l = h;
        }
    }
    
    ccp4spg_put_in_asu(low_group, _h, _k, _l, &h, &k, &l);
    
    int index = (h + OFFSET) * pow((double) MULTIPLIER, (int) 2)
    + (k + OFFSET) * MULTIPLIER + (l + OFFSET);
    
    return index;
}

bool MtzManager::holder_comparison(Holder *i, Holder *j)
{
    return (i->getReflId() < j->getReflId());
}

void MtzManager::insertionSortHolders(void)
{
    std::sort(holders.begin(), holders.end(), holder_comparison);
}

void MtzManager::sortLastHolder(void)
{
    int i, j;
    Holder *tmp;
    
    i = (int)holders.size() - 1;
    j = i;
    tmp = holders[i];
    while (j > 0 && tmp->getReflId() < holders[j - 1]->getReflId())
    {
        holders[j] = holders[j - 1];
        j--;
    }
    holders[j] = tmp;
    
    return;
    
    
}

void MtzManager::loadParametersMap()
{
    optimisingWavelength = FileParser::getKey("OPTIMISING_WAVELENGTH",
                                              !OPTIMISED_WAVELENGTH);
    optimisingBandwidth = FileParser::getKey("OPTIMISING_BANDWIDTH",
                                             !OPTIMISED_BANDWIDTH);
    optimisingMosaicity = FileParser::getKey("OPTIMISING_MOSAICITY",
                                             !OPTIMISED_MOSAICITY);
    optimisingOrientation = FileParser::getKey("OPTIMISING_ORIENTATION",
                                               !OPTIMISED_ROT);
    optimisingRlpSize = FileParser::getKey("OPTIMISING_RLP_SIZE",
                                           !OPTIMISED_SPOT_SIZE);
    optimisingExponent = FileParser::getKey("OPTIMISING_EXPONENT",
                                            !OPTIMISED_EXPONENT);
    
    wavelength = FileParser::getKey("INITIAL_WAVELENGTH", 0.0);
    usingFixedWavelength = (wavelength != 0);
    bandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
    mosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);
    spotSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
    
    exponent = FileParser::getKey("INITIAL_EXPONENT", INITIAL_EXPONENT);
    
    allowTrust = FileParser::getKey("ALLOW_TRUST", true);
    bool alwaysTrust = FileParser::getKey("TRUST_INDEXING_SOLUTION", false);
    
    if (alwaysTrust)
        trust = TrustLevelGood;
    
    stepSizeWavelength = FileParser::getKey("STEP_SIZE_WAVELENGTH",
                                            MEAN_STEP);
    stepSizeBandwidth = FileParser::getKey("STEP_SIZE_BANDWIDTH",
                                           BANDWIDTH_STEP);
    stepSizeMosaicity = FileParser::getKey("STEP_SIZE_MOSAICITY",
                                           MOS_STEP);
    stepSizeRlpSize = FileParser::getKey("STEP_SIZE_RLP_SIZE", SPOT_STEP);
    stepSizeOrientation = FileParser::getKey("STEP_SIZE_ORIENTATION",
                                             ROT_STEP);
    stepSizeExponent = FileParser::getKey("STEP_SIZE_EXPONENT",
                                          EXPONENT_STEP);
    
    toleranceWavelength = FileParser::getKey("TOLERANCE_WAVELENGTH",
                                             MEAN_TOLERANCE);
    toleranceBandwidth = FileParser::getKey("TOLERANCE_BANDWIDTH",
                                            BANDWIDTH_TOLERANCE);
    toleranceMosaicity = FileParser::getKey("TOLERANCE_MOSAICITY",
                                            MOS_TOLERANCE);
    toleranceRlpSize = FileParser::getKey("TOLERANCE_RLP_SIZE",
                                          SPOT_SIZE_TOLERANCE);
    toleranceOrientation = FileParser::getKey("TOLERANCE_ORIENTATION",
                                              ROT_TOLERANCE);
    toleranceExponent = FileParser::getKey("TOLERANCE_EXPONENT",
                                           EXPO_TOLERANCE);
    
    int defaultScoreInt = FileParser::getKey("DEFAULT_TARGET_FUNCTION",
                                             (int) DEFAULT_SCORE_TYPE);
    defaultScoreType = (ScoreType) defaultScoreInt;
    
    usePartialityFunction = FileParser::getKey("USE_PARTIALITY_FUNCTION",
                                               FURTHER_OPTIMISATION);
    
    maxResolutionAll = FileParser::getKey("MAX_RESOLUTION_ALL",
                                          MAX_OPTIMISATION_RESOLUTION);
    maxResolutionRlpSize = FileParser::getKey("MAX_RESOLUTION_RLP_SIZE",
                                              MAX_SPOT_SIZE_OPT_RESOLUTION);
}

MtzManager::MtzManager(void)
{
    lastReference = NULL;
    filename = "";
    inverted = 0;
    holders.resize(0);
    low_group = NULL;
    bandwidth = INITIAL_BANDWIDTH;
    hRot = 0;
    kRot = 0;
    mosaicity = INITIAL_MOSAICITY;
    spotSize = INITIAL_SPOT_SIZE;
    wavelength = 0;
    refCorrelation = 0;
    inverse = false;
    flipped = false;
    exponent = INITIAL_EXPONENT;
    finalised = false;
    scoreType = ScoreTypeCorrelation;
    trust = TrustLevelBad;
    maxResolutionAll = MAX_OPTIMISATION_RESOLUTION;
    maxResolutionRlpSize = MAX_SPOT_SIZE_OPT_RESOLUTION;
    freePass = false;
    defaultScoreType = DEFAULT_SCORE_TYPE;
    usePartialityFunction = FURTHER_OPTIMISATION;
    rejected = false;
    scale = 1;
    lastExponent = 0;
    previousReference = NULL;
    previousInverse = false;
    bFactor = 0;
    
    optimisingWavelength = !OPTIMISED_WAVELENGTH;
    optimisingBandwidth = !OPTIMISED_BANDWIDTH;
    optimisingMosaicity = !OPTIMISED_MOSAICITY;
    optimisingOrientation = !OPTIMISED_ROT;
    optimisingRlpSize = !OPTIMISED_SPOT_SIZE;
    optimisingExponent = !OPTIMISED_EXPONENT;
    
    stepSizeWavelength = MEAN_STEP;
    stepSizeBandwidth = BANDWIDTH_STEP;
    stepSizeMosaicity = MOS_STEP;
    stepSizeRlpSize = SPOT_STEP;
    stepSizeOrientation = ROT_STEP;
    stepSizeExponent = EXPONENT_STEP;
    
    toleranceWavelength = MEAN_TOLERANCE;
    toleranceBandwidth = BANDWIDTH_TOLERANCE;
    toleranceMosaicity = MOS_TOLERANCE;
    toleranceRlpSize = SPOT_SIZE_TOLERANCE;
    toleranceOrientation = ROT_TOLERANCE;
    toleranceExponent = EXPO_TOLERANCE;
    
    usingFixedWavelength = USE_FIXED_WAVELENGTH;
    
    matrix = MatrixPtr();
}

MtzManager *MtzManager::copy()
{
    MtzManager *newManager = new MtzManager();
    newManager->filename = filename;
    newManager->inverted = inverted;
    newManager->inverse = isInverse();
    double lowNum = low_group->spg_num;
    newManager->setSpaceGroup(lowNum);
    newManager->bandwidth = bandwidth;
    newManager->hRot = hRot;
    newManager->kRot = kRot;
    newManager->mosaicity = mosaicity;
    newManager->spotSize = spotSize;
    newManager->wavelength = wavelength;
    newManager->refCorrelation = refCorrelation;
    newManager->flipped = flipped;
    newManager->exponent = exponent;
    newManager->matrix = matrix;
    
    for (int i = 0; i < 3; i++)
    {
        newManager->cellDim[i] = cellDim[i];
        newManager->cellAngles[i] = cellAngles[i];
        
    }
    
    for (int i = 0; i < holders.size(); i++)
    {
        Holder *newHolder = holders[i]->copy();
        newManager->holders.push_back(newHolder);
    }
    
    return newManager;
}


void MtzManager::setUnitCell(double a, double b, double c, double alpha, double beta,
                             double gamma)
{
    cellDim[0] = a;
    cellDim[1] = b;
    cellDim[2] = c;
    cellAngles[0] = alpha;
    cellAngles[1] = beta;
    cellAngles[2] = gamma;
}

void MtzManager::setUnitCell(vector<double> unitCell)
{
    cellDim[0] = unitCell[0];
    cellDim[1] = unitCell[1];
    cellDim[2] = unitCell[2];
    cellAngles[0] = unitCell[3];
    cellAngles[1] = unitCell[4];
    cellAngles[2] = unitCell[5];
}

void MtzManager::clearHolders()
{
    for (int i = 0; i < holders.size(); i++)
    {
        delete holders[i];
    }
    
    holders.clear();
    vector<Holder *>().swap(holders);
}

void MtzManager::removeHolder(int i)
{
    delete holders[i];
    
    holders.erase(holders.begin() + i);
}

void MtzManager::addHolder(Holder *holder)
{
    Holder *newHolder = NULL;
    int insertionPoint = findHolderWithId(holder->getReflId(), &newHolder, true);
    
    holders.insert(holders.begin() + insertionPoint, holder);
    
 //   holders.push_back(holder);
 //   this->sortLastHolder();
}

void MtzManager::addHolders(vector<Holder *>holders)
{
    for (int i = 0; i < holders.size(); i++)
    {
        addHolder(holders[i]);
    }
}


int MtzManager::symmetryRelatedReflectionCount()
{
    int count = 0;
    
    for (int i = 0; i < holderCount(); i++)
    {
        if (holder(i)->millerCount() > 1)
            count++;
    }
    
    return count;
}

int MtzManager::refinedParameterCount()
{
    int count = 0;
    
    count += optimisingBandwidth;
    count += optimisingExponent;
    count += optimisingMosaicity;
    count += optimisingOrientation * 2;
    count += optimisingRlpSize;
    count += optimisingWavelength;
    
    return count;
}

void MtzManager::setMatrix(double *components)
{
    matrix = MatrixPtr(new Matrix(components));
}

void MtzManager::setMatrix(MatrixPtr newMat)
{
    matrix = newMat;
}

void MtzManager::setDefaultMatrix()
{
    double a = cellDim[0];
    double b = cellDim[1];
    double c = cellDim[2];
    double alpha = cellAngles[0];
    double beta = cellAngles[1];
    double gamma = cellAngles[2];
    
    Matrix mat = Matrix::matrixFromUnitCell(a, b, c, alpha, beta, gamma);
    matrix = mat.copy();
}

void MtzManager::getUnitCell(double *a, double *b, double *c, double *alpha,
                             double *beta, double *gamma)
{
    *a = cellDim[0];
    *b = cellDim[1];
    *c = cellDim[2];
    *alpha = cellAngles[0];
    *beta = cellAngles[1];
    *gamma = cellAngles[2];
}

void MtzManager::copySymmetryInformationFromManager(MtzPtr toCopy)
{
    double a, b, c, alpha, beta, gamma;
    toCopy->getUnitCell(&a, &b, &c, &alpha, &beta, &gamma);
    this->setUnitCell(a, b, c, alpha, beta, gamma);
    int spgnum = toCopy->getLowGroup()->spg_ccp4_num;
    CCP4SPG *spg = ccp4spg_load_by_ccp4_num(spgnum);
    this->setLowGroup(spg);
}

void MtzManager::loadReflections(int partials)
{
    PartialityModel model =
    partials ? PartialityModelScaled : PartialityModelSimple;
    
    loadReflections(model);
}

void MtzManager::loadReflections(PartialityModel model)
{
    if (filename.length() == 0)
    {
        std::cerr
        << "Cannot load reflections as no filename has been specified."
        << endl;
        return;
    }
    
    MTZ *mtz = MtzGet(filename.c_str(), 0);
    
    int spgnum = MtzSpacegroupNumber(mtz);
    
    low_group = ccp4spg_load_by_ccp4_num(spgnum);
    
    if (mtz == NULL)
        return;
    
    float *refldata = (float *) CCP4::ccp4_utils_malloc(
                                                        (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    memset(refldata, '\0',
           (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
    
    float *adata = (float *) CCP4::ccp4_utils_malloc(
                                                     (mtz->ncol_read) * sizeof(float));
    
    MtzRrefl(mtz->filein, mtz->ncol_read * mtz->nref_filein, refldata);
    
    MTZCOL *col_f = MtzColLookup(mtz, "I");
    
    if (col_f == NULL)
    {
        col_f = MtzColLookup(mtz, "IMEAN");
    }
    
    if (col_f == NULL)
    {
        std::cout << "Warning: intensity column not labelled I or IMEAN" << std::endl;
        exit(1);
    }
    
    MTZCOL *col_sigf = MtzColLookup(mtz, "SIGI");
    
    if (col_sigf == NULL)
    {
        col_sigf = MtzColLookup(mtz, "SIGIMEAN");
    }
    
    if (col_sigf == NULL)
    {
        std::cout << "Warning: sigma column not labelled SIGI or SIGIMEAN" << std::endl;
        exit(1);
    }
    
    MTZCOL *col_wave = MtzColLookup(mtz, "WAVE");
    MTZCOL *col_partials = MtzColLookup(mtz, "PART");
    MTZCOL *col_shiftx = MtzColLookup(mtz, "SHIFTX");
    MTZCOL *col_shifty = MtzColLookup(mtz, "SHIFTY");
    
    int multiplier = MULTIPLIER;
    
    MTZXTAL **xtals = MtzXtals(mtz);
    float cell[6];
    double coefhkl[6];
    ccp4_lrcell(xtals[0], cell);
    MtzHklcoeffs(cell, coefhkl);
    
    setUnitCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    
    int num = 0;
    
    for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
    {
        memcpy(adata, &refldata[i], mtz->ncol_read * sizeof(float));
        
        int h;
        int k;
        int l;
        int offset;
        
        hkls_for_reflection(mtz, adata, &h, &k, &l, &multiplier, &offset);
        
        int reflection_holder_index = index_for_reflection(h, k, l, false);
        int inv_holder_index = index_for_reflection(h, k, l, true);
        
        float intensity = adata[col_f->source - 1];
        float sigma = adata[col_sigf->source - 1];
        float wavelength = 1;
        float partiality = 1;
        float shiftX = 0;
        float shiftY = 0;
        
        if (col_wave != NULL && col_partials != NULL)
        {
            wavelength = adata[col_wave->source - 1];
            partiality = adata[col_partials->source - 1];
        }
        
        if (col_shiftx != NULL && col_shifty != NULL)
        {
            shiftX = adata[col_shiftx->source - 1];
            shiftY = adata[col_shifty->source - 1];
        }
        
        MillerPtr miller = MillerPtr(new Miller(this, h, k, l));
        miller->setData(intensity, sigma, partiality, wavelength);
        miller->setCountingSigma(sigma);
        miller->setFilename(filename);
        miller->setPartialityModel(model);
        miller->setShift(make_pair(shiftX, shiftY));
        miller->matrix = this->matrix;
        
        Holder *prevHolder;
        
        this->findHolderWithId(reflection_holder_index, &prevHolder);
        
        if (prevHolder != NULL)
        {
            /** Exclude unobserved reflections by checking for nan */
            if (adata[col_f->source - 1] == adata[col_f->source - 1])
            {
                prevHolder->addMiller(miller); // TODO
                miller->setParent(prevHolder);
            }
            
            // reflection is a repeat so set flag.
        }
        
        if (prevHolder == NULL)
        {
            Holder *newHolder = new Holder();
            holders.push_back(newHolder);
            
            holders[num]->addMiller(miller);
            miller->setParent(holders[num]);
            holders[num]->calculateResolution(this);
            
            holders[num]->setReflId(reflection_holder_index);
            holders[num]->setInvReflId(inv_holder_index);
            
            this->sortLastHolder();
            
            num++;
        }
    }
    
    free(refldata);
    free(adata);
    
    //	insertionSortHolders();
    
    std::ostringstream log;
    
    log << "Loaded " << mtz->nref_filein << " reflections (" << accepted()
    << " accepted)." << endl;
    
    Logger::mainLogger->addStream(&log);
    
    MtzFree(mtz);
}

void MtzManager::getRefHolders(vector<Holder *> *refPointer,
                               vector<Holder *> *matchPointer)
{
    if (lastReference != referenceManager)
    {
        refHolders.clear();
        matchHolders.clear();
        
        for (int i = 0; i < holderCount(); i++)
        {
            int reflId = holder(i)->getReflId();
            if (isInverse())
                reflId = holder(i)->getInvReflId();
            
            Holder *refHolder = NULL;
            referenceManager->findHolderWithId(reflId, &refHolder);
            
            if (refHolder != NULL)
            {
                matchHolders.push_back(holder(i));
                refHolders.push_back(refHolder);
            }
        }
        
        lastReference = referenceManager;
    }
    
    refPointer->reserve(refHolders.size());
    refPointer->insert(refPointer->begin(), refHolders.begin(),
                       refHolders.end());
    
    matchPointer->reserve(matchHolders.size());
    matchPointer->insert(matchPointer->begin(), matchHolders.begin(),
                         matchHolders.end());
    
}

void MtzManager::setReference(MtzManager *reference)
{
    MtzManager::referenceManager = reference;
}

void MtzManager::setFilename(string name)
{
    filename = name;
}

string MtzManager::getFilename(void)
{
    return filename;
}

void MtzManager::setSpaceGroup(int spgnum)
{
    low_group = ccp4spg_load_by_ccp4_num(spgnum);
}

int MtzManager::findHolderWithId(int refl_id, Holder **holder, bool insertionPoint)
{
    if (holderCount() == 0)
    {
        *holder = NULL;
        return 0;
    }
    
    int lower = 0;
    int higher = holderCount() - 1;
    int new_bound = (higher + lower) / 2;
    
    if ((refl_id < this->holder(lower)->getReflId())
        || (refl_id > this->holder(higher)->getReflId()))
    {
        if (insertionPoint)
        {
            if (refl_id < this->holder(lower)->getReflId())
                return 0;
            
            if (refl_id > this->holder(higher)->getReflId())
                return holderCount();
        }
        
        *holder = NULL;
        return -1;
    }
    
    while (this->holder(new_bound)->getReflId() != refl_id)
    {
        if (new_bound == higher || new_bound == lower)
        {
            if (this->holder(higher)->getReflId() == refl_id)
            {
                (*holder) = this->holder(higher);
                return -1;
            }
            
            if (this->holder(lower)->getReflId() == refl_id)
            {
                (*holder) = this->holder(lower);
                return -1;
            }
            
            if (insertionPoint)
            {
                int lowest = 0;
                
                int start = lower - 2 >= 0 ? lower - 2 : 0;
                
                for (int i = start; i < higher + 2 && i < holderCount(); i++)
                {
                    if (this->holder(i)->getReflId() < refl_id)
                        lowest = i;
                }
                
                return lowest + 1;
            }
            
            *holder = NULL;
            return -1;
        }
        
        if (this->holder(new_bound)->getReflId() > refl_id)
        {
            higher = new_bound;
        }
        else if (this->holder(new_bound)->getReflId() < refl_id)
        {
            lower = new_bound;
        }
        
        new_bound = (higher + lower) / 2;
    }
    
    (*holder) = this->holder(new_bound);
    
    return -1;
}

void MtzManager::findCommonReflections(MtzManager *other,
                                       vector<Holder *> &holderVector1, vector<Holder *> &holderVector2,
                                       int *num, int inverted, bool excludeSymmetrical)
{
    if (other == previousReference && previousInverse == inverted)
    {
        holderVector1.reserve(matchHolders.size());
        holderVector2.reserve(refHolders.size());
        
        holderVector1.insert(holderVector1.begin(), matchHolders.begin(), matchHolders.end());
        holderVector2.insert(holderVector2.begin(), refHolders.begin(), refHolders.end());
        
        if (num != NULL)
            *num = (int)holderVector1.size();
        
        return;
    }

    matchHolders.clear();
    refHolders.clear();
    
    previousReference = other;
    previousInverse = inverted;
    
    for (int i = 0; i < holderCount(); i++)
    {
        if (excludeSymmetrical
            && holder(i)->getReflId() == holder(i)->getInvReflId())
            continue;
        
        int refl_id = holder(i)->getReflId();
        
        if (inverted)
            refl_id = holder(i)->getInvReflId();
        
        Holder *otherHolder = NULL;
        
        other->findHolderWithId(refl_id, &otherHolder);
        
        if (otherHolder != NULL)
        {
            holderVector1.push_back(holder(i));
            matchHolders.push_back(holder(i));
            holderVector2.push_back(otherHolder);
            refHolders.push_back(otherHolder);
        }
    }
    
    if (num != NULL)
    {
        *num = (int)holderVector1.size();
    }
}

void MtzManager::applyScaleFactorsForBins()
{
    vector<double> bins;
    StatisticsManager::generateResolutionBins(50, 1.4, 5, &bins);
    
    for (int shell = 0; shell < bins.size() - 1; shell++)
    {
        double low = bins[shell];
        double high = bins[shell + 1];
        
        vector<Holder *> refHolders, imgHolders;
        
        this->findCommonReflections(referenceManager, imgHolders, refHolders,
                                    NULL, isInverse());
        
        double weights = 0;
        double refMean = 0;
        double imgMean = 0;
        int count = 0;
        
        for (int i = 0; i < imgHolders.size(); i++)
        {
            if (!imgHolders[i]->anyAccepted())
                continue;
            
            if (imgHolders[i]->betweenResolutions(low, high))
            {
                weights += imgHolders[i]->meanPartiality();
                refMean += refHolders[i]->meanIntensity()
                * imgHolders[i]->meanPartiality();
                imgMean += imgHolders[i]->meanIntensity()
                * imgHolders[i]->meanPartiality();
                count++;
            }
        }
        
        refMean /= weights;
        imgMean /= weights;
        
        double ratio = refMean / imgMean;
        
        /*	std::cout << 1 / pow(low, 2) << "\t" << ratio << "\t" << grad << "\t"
         << count << "\t" << std::endl;
         */
        applyScaleFactor(ratio, low, high);
    }
}

void MtzManager::bFactorAndScale(double *scale, double *bFactor, double exponent,
                                 vector<pair<double, double> > *dataPoints)
{
    clearScaleFactor();
    MtzManager *reference = MtzManager::getReferenceManager();
    
    double grad = this->gradientAgainstManager(*reference);
    this->applyScaleFactor(grad);
    
    vector<Holder *> refHolders, imageHolders;
    this->findCommonReflections(reference, imageHolders, refHolders, NULL,
                                isInverse());
    
    vector<double> smoothBins;
    StatisticsManager::generateResolutionBins(50, 1.8, 24, &smoothBins);
    
    if (dataPoints == NULL)
    {
        vector<pair<double, double> > data = vector<pair<double, double> >();
        dataPoints = &data;
    }
    
    vector<boost::tuple<double, double, double> > pointsToFit;
    
    for (int shell = 0; shell < smoothBins.size() - 3; shell++)
    {
        double low = smoothBins[shell];
        double high = smoothBins[shell + 3];
        
        vector<Holder *> refHolders, imgHolders;
        
        this->findCommonReflections(referenceManager, imgHolders, refHolders,
                                    NULL, isInverse());
        
        double weights = 0;
        double refMean = 0;
        double imgMean = 0;
        int count = 0;
        
        for (int i = 0; i < imgHolders.size(); i++)
        {
            if (!imgHolders[i]->anyAccepted())
                continue;
            
            if (imgHolders[i]->betweenResolutions(low, high))
            {
                weights += imgHolders[i]->meanPartiality();
                
                refMean += refHolders[i]->meanIntensity()
                * imgHolders[i]->meanPartiality();
                
                imgMean += imgHolders[i]->meanIntensity()
                * imgHolders[i]->meanPartiality();
                count++;
            }
        }
        
        refMean /= weights;
        imgMean /= weights;
        
        double ratio = refMean / imgMean;
        ratio = 1 / ratio;
        
        if (ratio != ratio)
            continue;
        
        double resolution = StatisticsManager::midPointBetweenResolutions(
                                                                          smoothBins[shell], smoothBins[shell + 3]);
        
        double res_squared = pow(1 / resolution, 2);
        
        double four_d_squared = 4 * pow(resolution, 2);
        
        double right_exp = pow(1 / four_d_squared, exponent);
        
        double four_d_to_exp = pow(2, exponent) * right_exp;
        
        double intensityRatio = (refMean / imgMean);
        double logIntensityRatio = log(intensityRatio);
        
        if (logIntensityRatio != logIntensityRatio)
            continue;
        
        double weight = 1;
        weight /= res_squared;
        
        boost::tuple<double, double, double> point = boost::make_tuple(four_d_to_exp,
                                                              logIntensityRatio, weight);
        
        pair<double, double> showPoint = std::make_pair(res_squared,
                                                        1 / intensityRatio);
        
        dataPoints->push_back(showPoint);
        pointsToFit.push_back(point);
    }
    
    double gradient = 0;
    double intercept = 0;
    
    regression_line(pointsToFit, intercept, gradient);
    
    double k = exp(intercept);
    
    double b = pow(abs(gradient), 1 / (double) exponent);
    
    if (gradient < 0)
        b = -b;
    
    if (b != b)
        b = 0;
    
    *scale = k;
    *bFactor = b;
}

double MtzManager::minimizeRFactor(MtzManager *otherManager)
{
    return minimizeGradient(otherManager, false);
}

double MtzManager::minimizeGradient(MtzManager *otherManager, bool leastSquares)
{
    double resolution = 0.001;
    double step = 0.1;
    double gradient = gradientAgainstManager(*otherManager);
    
    vector<Holder *> holders1;
    vector<Holder *> holders2;
    vector<double> ints1;
    vector<double> ints2;
    vector<double> weights;
    int num = 0;
    
    MtzManager::findCommonReflections(otherManager, holders1, holders2, &num,
                                      isInverse());
    
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < holders1[i]->millerCount(); j++)
        {
            if (!holders1[i]->miller(j)->accepted())
                continue;
            
            double mean1 = holders1[i]->miller(j)->intensity();
            double mean2 = holders2[i]->meanIntensity();
            double weight = holders2[i]->meanPartiality();
            
            if (mean1 != mean1 || mean2 != mean2 || weight != weight)
                continue;
            
            //		std::cout << mean1 << "\t" << mean2 << std::endl;
            
            ints1.push_back(mean1);
            ints2.push_back(mean2);
            weights.push_back(weight);
        }
    }
    
    double best = r_factor_between_vectors(&ints1, &ints2, &weights, gradient);
    ;
    
    while (step > resolution)
    {
        double gradUp = gradient + step;
        double gradDown = gradient - step;
        
        double upScore = r_factor_between_vectors(&ints1, &ints2, &weights,
                                                  gradUp);
        double downScore = r_factor_between_vectors(&ints1, &ints2, &weights,
                                                    gradDown);
        
        //	std::cout << upScore << "\t" << best << "\t" << downScore << std::endl;
        
        if (upScore < best)
        {
            gradient = gradUp;
            best = upScore;
        }
        else if (downScore < best)
        {
            gradient = gradDown;
            best = downScore;
        }
        else
        {
            step /= 2;
        }
    }
    
    return gradient;
}

double MtzManager::gradientAgainstManager(MtzManager &otherManager,
                                          bool leastSquares, double lowRes, double highRes)
{
    vector<Holder *> holders1;
    vector<Holder *> holders2;
    
    vector<double> ints1, ints2;
    int num = 0;
    
    double minD = 0;
    double maxD = 0;
    StatisticsManager::convertResolutions(lowRes, highRes, &minD, &maxD);
    
    MtzManager::findCommonReflections(&otherManager, holders1, holders2, &num,
                                      isInverse());
    
    if (num <= 1)
        return 1;
    
    double x_squared = 0;
    double x_y = 0;
    
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < holders1[i]->millerCount(); j++)
        {
            if (!holders1[i]->miller(j)->accepted())
                continue;
            
            if (holders1[i]->getResolution() < minD
                || holders1[i]->getResolution() > maxD)
                continue;
            
            double mean1 = holders1[i]->miller(j)->intensity();
            double mean2 = holders2[i]->meanIntensity();
            
            double part1 = holders1[i]->miller(j)->getPartiality();
            
            if (mean1 != mean1 || mean2 != mean2)
                continue;
            
            x_squared += mean1 * mean1 * part1;
            x_y += mean1 * mean2 * part1;
            
            if (holders1[i]->miller(j)->getPartiality() < 0.5)
                continue;
            
            ints1.push_back(mean1);
            ints2.push_back(mean2);
            
        }
    }
    
    double grad = x_y / x_squared;
    
    if (grad < 0)
        grad = -1;
    
    if (leastSquares)
        grad = minimize_gradient_between_vectors(&ints1, &ints2);
    
    return grad;
}

void MtzManager::clearScaleFactor()
{
    for (int i = 0; i < holders.size(); i++)
    {
        for (int j = 0; j < holders[i]->millerCount(); j++)
        {
            holders[i]->miller(j)->setScale(1);
            holders[i]->miller(j)->setBFactor(0);
        }
    }
}

void MtzManager::makeScalesPermanent()
{
    for (int i = 0; i < holders.size(); i++)
    {
        for (int j = 0; j < holder(i)->millerCount(); j++)
        {
            holder(i)->miller(j)->makeScalesPermanent();
        }
    }
}

void MtzManager::applyBFactor(double bFactor)
{
 //   if (bFactor < 0)
 //       bFactor = 0 - bFactor;
    
    for (int i = 0; i < holders.size(); i++)
    {
        for (int j = 0; j < holders[i]->millerCount(); j++)
        {
            holder(i)->miller(j)->setBFactor(bFactor);
        }
    }
}

void MtzManager::applyScaleFactor(double scaleFactor,
                                  double lowRes, double highRes)
{
    if (scaleFactor == scaleFactor)
        scale *= scaleFactor;
    
    for (int i = 0; i < holders.size(); i++)
    {
        for (int j = 0; j < holders[i]->millerCount(); j++)
        {
            if (holder(i)->betweenResolutions(lowRes, highRes))
                holder(i)->miller(j)->applyScaleFactor(scaleFactor);
        }
    }
}

double MtzManager::averageIntensity(void)
{
    double total_intensity = 0;
    double total = 0;
    
    for (int i = 0; i < holders.size(); i++)
    {
        double intensity = holders[i]->meanIntensity();
        if (intensity == intensity)
        {
            double weight = holders[i]->meanPartiality();
            total_intensity += intensity * weight;
            total += weight;
        }
    }
    
    total_intensity /= total;
    return total_intensity;
}

void MtzManager::applyPolarisation(void)
{
    for (int i = 0; i < holders.size(); i++)
    {
        for (int j = 0; j < holders[i]->millerCount(); j++)
        {
            holders[i]->miller(j)->applyPolarisation(1.459);
        }
    }
}

void MtzManager::writeToFile(string newFilename, bool announce, bool shifts)
{
    int columns = 7;
    
    if (shifts) columns += 2;
    
    float cell[6], wavelength;
    float *fdata = new float[columns];
    
    bool flipped = false;
    
    if (isInverse())
    {
        //	flip();
        //	flipped = true;
    }
    
    /* variables for symmetry */
    CCP4SPG *mtzspg = low_group;
    float rsm[192][4][4];
    char ltypex[2];
    
    /* variables for MTZ data structure */
    MTZ *mtzout;
    MTZXTAL *xtal;
    MTZSET *set;
    MTZCOL *colout[9];
    
    
    /*  Removed: General CCP4 initializations e.g. HKLOUT on command line */
    
    cell[0] = cellDim[0];
    cell[1] = cellDim[1];
    cell[2] = cellDim[2];
    cell[3] = cellAngles[0];
    cell[4] = cellAngles[1];
    cell[5] = cellAngles[2];
    wavelength = this->getWavelength();
    
    mtzout = MtzMalloc(0, 0);
    ccp4_lwtitl(mtzout, "Written from Helen's XFEL tasks ", 0);
    mtzout->refs_in_memory = 0;
    mtzout->fileout = MtzOpenForWrite(newFilename.c_str());
    
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
    colout[5] = MtzAddColumn(mtzout, set, "PART", "R");
    colout[6] = MtzAddColumn(mtzout, set, "WAVE", "R");
    
    if (shifts)
    {
        colout[7] = MtzAddColumn(mtzout, set, "SHIFTX", "R");
        colout[8] = MtzAddColumn(mtzout, set, "SHIFTY", "R");
    }
    
    int num = 0;
    
    for (int i = 0; i < holderCount(); i++)
    {
        for (int j = 0; j < holders[i]->millerCount(); j++)
        {
            if (!holder(i)->miller(j))
                std::cout << "!miller(j) in mtz manager" << std::endl;
            
            if (holders[i]->miller(j)->isRejected())
                continue;
            
            double intensity = holders[i]->miller(j)->getRawIntensity();
            double sigma = holders[i]->miller(j)->getSigma();
            double partiality = holders[i]->miller(j)->getPartiality();
            double bFactor = holders[i]->miller(j)->getBFactorScale();
            
            if (intensity != intensity)
            {
                continue;
            }
            
            num++;
            
            int h = holders[i]->miller(j)->h;
            int k = holders[i]->miller(j)->k;
            int l = holders[i]->miller(j)->l;
            int _h, _k, _l;
            ccp4spg_put_in_asu(low_group, h, k, l, &_h, &_k, &_l);
            
            // set fdata
            fdata[0] = h;
            fdata[1] = k;
            fdata[2] = l;
            fdata[3] = intensity / bFactor;
            fdata[4] = sigma;
            fdata[5] = partiality;
            fdata[6] = holders[i]->miller(j)->getWavelength();
            
            if (shifts)
            {
                fdata[7] = holders[i]->miller(j)->getShift().first;
                fdata[8] = holders[i]->miller(j)->getShift().second;
            }
            
            ccp4_lwrefl(mtzout, fdata, colout, columns, num);
        }
    }
    
    if (flipped)
        flip();
    
    MtzPut(mtzout, " ");
    MtzFree(mtzout);
    
    LogLevel shouldAnnounce = announce ? LogLevelNormal : LogLevelDebug;
    
    ostringstream logged;
    logged << "Written to file " << newFilename << std::endl;
    Logger::mainLogger->addStream(&logged, shouldAnnounce);
    
    delete [] fdata;
}

void MtzManager::flip(void)
{
    for (int i = 0; i < holders.size(); i++)
    {
        holders[i]->flip(this);
    }
    
    flipped = !flipped;
    
    this->insertionSortHolders();
}

MtzManager::~MtzManager(void)
{
    for (int i = 0; i < holders.size(); i++)
    {
        delete holders[i];
    }
    
    holders.clear();
    vector<Holder *>().swap(holders);
    
    if (low_group != NULL)
        ccp4spg_free(&low_group);
    
}

void MtzManager::description(void)
{
    logged << "Filename: " << filename << std::endl;
    logged << "Number of holders: " << holderCount() << std::endl;
    logged << "Number of accepted Millers: " << accepted() << std::endl;
    logged << "Average intensity: " << this->averageIntensity() << std::endl;
    
    sendLog();
}

void MtzManager::setSigmaToUnity()
{
    for (int i = 0; i < holders.size(); i++)
    {
        for (int j = 0; j < holders[i]->millerCount(); j++)
        {
            holders[i]->miller(j)->setSigma(1);
            
            // TAKE NEXT LINE OUT
        }
    }
}

void MtzManager::setPartialityToUnity()
{
    for (int i = 0; i < holders.size(); i++)
    {
        for (int j = 0; j < holders[i]->millerCount(); j++)
        {
            holders[i]->miller(j)->setPartiality(1);
            
            // TAKE NEXT LINE OUT
        }
    }
}

int MtzManager::accepted(void)
{
    int acceptedCount = 0;
    
    for (int i = 0; i < holderCount(); i++)
    {
        for (int j = 0; j < holder(i)->millerCount(); j++)
        {
            if (holder(i)->miller(j)->accepted())
                acceptedCount++;
        }
    }
    
    return acceptedCount;
}

void MtzManager::invert()
{
    inverse = !inverse;
    
    for (int i = 0; i < holderCount(); i++)
    {
        holder(i)->setInverse(inverse);
    }
}

void MtzManager::writeToDat()
{
    std::string name = filename;
    int lastindex = (int)name.find_last_of(".");
    std::string rootName = name.substr(0, lastindex);
    std::string datName = rootName + ".dat";
    
    std::ofstream datOutput;
    datOutput.open(datName);
    
    datOutput << this->matrix->description() << std::endl;
    
    for (int i = 0; i < holderCount(); i++)
    {
        Holder *aHolder = holder(i);
        for (int j = 0; j < aHolder->millerCount(); j++)
        {
            MillerPtr miller = aHolder->miller(j);
            double shiftX = miller->getShift().first;
            double shiftY = miller->getShift().second;
            double lastX = miller->getLastX();
            double lastY = miller->getLastY();
            
            double combinedX = shiftX + lastX;
            double combinedY = shiftY + lastY;
            
            datOutput << miller->h << "\t" << miller->k << "\t" << miller->l;
            
            datOutput << "\t" << miller->getWavelength();
            datOutput << "\t" << miller->getRawIntensity();
            datOutput << "\t" << miller->getCountingSigma();
            datOutput << "\t" <<  combinedX << "\t" << combinedY << "\t"
            << miller->getPartiality() << "\t" << miller->getBFactorScale() << "\t" << miller->getResolution() << std::endl;
        }
    }
    
    
    datOutput.close();
}

void MtzManager::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

int MtzManager::rejectOverlaps()
{
    int count = 0;
    
    for (int i = 0; i < holderCount(); i++)
    {
        count += holder(i)->checkOverlaps();
    }
    
    return count;
}

