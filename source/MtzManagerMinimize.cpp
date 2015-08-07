#include "StatisticsManager.h"
#include <cmath>
#include "Vector.h"
#include <algorithm>
#include "Scaler.h"

#include "MtzRefiner.h"
#include "FileParser.h"
#include "MtzManager.h"
#include "parameters.h"
#include "GraphDrawer.h"
#include "NelderMead.h"


#define DEFAULT_RESOLUTION 3.4
typedef boost::tuple<double, double, double, double, double> ResultTuple;

double MtzManager::weightedBestWavelength(double lowRes, double highRes)
{
    vector<double> wavelengths;
    vector<double> percentages;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        Reflection *imageReflection = reflection(i);
        
        if (imageReflection->getResolution() < lowRes)
            continue;
        
        if (imageReflection->getResolution() > highRes)
            continue;
        
        Reflection *refReflection = NULL;
        int reflid = imageReflection->getReflId();
        
        referenceManager->findReflectionWithId(reflid, &refReflection);
        
        if (refReflection != NULL)
        {
            double wavelength = imageReflection->miller(0)->getWavelength();
            
            double percentage = imageReflection->miller(0)->getRawIntensity()
            / imageReflection->miller(0)->getCountingSigma();
            
            if (percentage != percentage)
                continue;
            
            wavelengths.push_back(wavelength);
            percentages.push_back(percentage);
        }
    }
    
    return weighted_mean(&wavelengths, &percentages);
}

double MtzManager::unitCellScore(void *object)
{
    static_cast<MtzManager *>(object)->refreshPartialities(static_cast<MtzManager *>(object)->params);
    
    return exclusionScoreWrapper(object, 0, 0);
}

double MtzManager::bFactorScoreWrapper(void *object)
{
    vector<double> bins;
    StatisticsManager::generateResolutionBins(0, static_cast<MtzManager *>(object)->maxResolutionAll, 1, &bins);
    
    static_cast<MtzManager *>(object)->applyBFactor(static_cast<MtzManager *>(object)->bFactor);
    
    double score = 0;
    
    for (int i = 0; i < bins.size() - 1; i++)
    {
        score += exclusionScoreWrapper(object, bins[i], bins[i + 1]);
        
    }
    
    return score;
}

double MtzManager::scoreNelderMead(void *object)
{
    static_cast<MtzManager *>(object)->refreshCurrentPartialities();
    return static_cast<MtzManager *>(object)->exclusionScoreWrapper(object, 0, 0);
}

double MtzManager::exclusionScoreWrapper(void *object, double lowRes,
                                         double highRes)
{
    ScoreType scoreType = static_cast<MtzManager *>(object)->scoreType;
    
    if (scoreType == ScoreTypeCorrelation
        || scoreType == ScoreTypeCorrelationLog)
    {
        return static_cast<MtzManager *>(object)->exclusionScore(lowRes,
                                                                 highRes, scoreType);
    }
    else if (scoreType == ScoreTypeMinimizeRSplit)
    {
        double rSplit = static_cast<MtzManager *>(object)->rSplit(lowRes, highRes);
        double penalty = static_cast<MtzManager *>(object)->belowPartialityPenalty(0, static_cast<MtzManager *>(object)->penaltyResolution);
        double penaltyWeight = static_cast<MtzManager *>(object)->penaltyWeight;
        
        return rSplit - penaltyWeight * penalty;
    }
    else if (scoreType == ScoreTypeMinimizeRSplitLog)
    {
        return static_cast<MtzManager *>(object)->rSplit(lowRes, highRes, true);
    }
    else if (scoreType == ScoreTypeSymmetry)
    {
        return static_cast<MtzManager *>(object)->rFactorWithManager(RFactorTypeMeas);
    }
    else if (scoreType == ScoreTypePartialityCorrelation
             || scoreType == ScoreTypePartialityLeastSquares)
    {
        return 1
        - static_cast<MtzManager *>(object)->leastSquaresPartiality(
                                                                    scoreType);
    }
    else if (scoreType == ScoreTypeStandardDeviation)
    {
        return static_cast<MtzManager *>(object)->wavelengthStandardDeviation();
    }
    else if (scoreType == ScoreTypeMinimizeRMeas)
    {
        Scaler *scaler = Scaler::getScaler();
        double rSplit = static_cast<MtzManager *>(object)->rSplit(lowRes, highRes);
        if (scaler == NULL)
        {
            return rSplit;
        }
        else
        {
            return scaler->evaluateForImage(static_cast<MtzManager *>(object)) + rSplit / 2;
        }
    }
    
    else
        return static_cast<MtzManager *>(object)->exclusionScore(lowRes,
                                                                 highRes, ScoreTypeCorrelation);
}

double MtzManager::belowPartialityPenalty(double low, double high)
{
    double partialPenalty = 0;
    double weights = 0;
    
    for (int i = 0; i < reflections.size(); i++)
    {
        Reflection *imageReflection = reflections[i];
        
        if (imageReflection->anyAccepted())
            continue;
        
        if (!imageReflection->betweenResolutions(low, high))
            continue;
        
        Reflection *refReflection = NULL;
        int reflid = imageReflection->getReflId();
        
        referenceManager->findReflectionWithId(reflid, &refReflection);
        
        if (refReflection != NULL)
        {
            double rawIntensity = imageReflection->miller(0)->getRawIntensity();
            double meanIntensity = refReflection->meanIntensity();
            
            if (rawIntensity != rawIntensity || meanIntensity != meanIntensity)
                continue;
            
            double ratio = rawIntensity / meanIntensity;
            ratio -= imageReflection->miller(0)->getPartialCutoff();
            
            partialPenalty += ratio;
            weights++;
        }
    }
    
    partialPenalty /= weights;
    
 //   std::cout << "Partial penalty: " << partialPenalty << std::endl;
    
    return partialPenalty;
}

double MtzManager::rSplit(double low, double high, bool square)
{
    double scale = this->gradientAgainstManager(referenceManager);
    
    applyScaleFactor(scale);
    
    double sum_numerator = 0;
    double sum_denominator = 0;
    int count = 0;
    double weights = 0;
    
    double lowCut = low == 0 ? 0 : 1 / low;
    double highCut = high == 0 ? FLT_MAX : 1 / high;
    
    vector<Reflection *> reflections1;
    vector<Reflection *> reflections2;
    
    this->findCommonReflections(referenceManager, reflections1, reflections2, NULL);
    
    for (int i = 0; i < reflections1.size(); i++)
    {
        Reflection *reflection = reflections1[i];
        Reflection *reflection2 = reflections2[i];
        
        if (reflection->getResolution() < lowCut
            || reflection->getResolution() > highCut)
            continue;
        
        double int1 = reflection->meanIntensity();
        if (reflection2->meanIntensity() < REFERENCE_WEAK_REFLECTION)
            continue;

        double int2 = reflection2->meanIntensityWithExclusion(&filename);
        
        double weight = reflection->meanWeight() * reflection->getResolution();
        
        if (int1 != int1 || int2 != int2)
            continue;
        
        if (int1 + int2 < 0)
            continue;
        
        count++;
        weights += weight;
        
        if (!square)
        {
            sum_numerator += fabs(int1 - int2) * weight;
            sum_denominator += (int1 + int2) * weight / 2;
        }
        else
        {
            double averageIntensity = pow((int1 + int2) / 2, 2);
            double squareAddition = pow(int1 - int2, 2);
            
            sum_numerator += squareAddition;
            sum_denominator += averageIntensity;
        }
        
    }
    
    if (square)
    {
        //   sum_numerator = sqrt(sum_numerator);
    }
    
    double r_split = sum_numerator / (sum_denominator * sqrt(2));
    
    return r_split;
    
}

double MtzManager::exclusionScore(double lowRes, double highRes,
                                  ScoreType scoreType)
{
    double score = 0;
    
    bool shouldLog = (scoreType == ScoreTypeCorrelationLog);
    
    double correlation = this->correlationWithManager(referenceManager, false,
                                                      true, lowRes, highRes, 1, NULL, shouldLog);
    
    score += 1 - correlation;
    
    return score;
}

double MtzManager::leastSquaresPartiality(ScoreType typeOfScore)
{
    double score = leastSquaresPartiality(0, 2.5, typeOfScore);
    
    return score;
}

bool percentage(Partial x, Partial y)
{
    return (x.percentage < y.percentage);
}

double MtzManager::leastSquaresPartiality(double low, double high,
                                          ScoreType typeOfScore)
{
    vector<Partial> partials;
    vector<double> partialities;
    vector<double> percentages;
    
    
    for (int i = 0; i < reflections.size(); i++)
    {
        Reflection *imageReflection = reflections[i];
        Reflection *refReflection = NULL;
        int reflid = imageReflection->getReflId();
        
        referenceManager->findReflectionWithId(reflid, &refReflection);
        
        if (refReflection != NULL)
        {
            if (refReflection->meanIntensity() < REFERENCE_WEAK_REFLECTION)
                continue;
            
            if (!refReflection->betweenResolutions(low, high))
                continue;
            
            Partial partial;
            partial.partiality = imageReflection->miller(0)->getPartiality();
            partial.percentage = imageReflection->miller(0)->getRawIntensity()
            / refReflection->meanIntensityWithExclusion(&filename);
            partial.resolution = imageReflection->getResolution();
            
            partials.push_back(partial);
        }
    }
    
    std::sort(partials.begin(), partials.end(), percentage);
    int position = (int)partials.size() - 5;
    if (position <= 0) return 0;
    double maxPercentage = 2.5;
    double minPartiality = FileParser::getKey("PARTIALITY_CUTOFF", PARTIAL_CUTOFF);
    
    for (int i = 0; i < partials.size(); i++)
    {
        if (partials[i].percentage > maxPercentage || partials[i].percentage < 0)
            continue;
        
        if (partials[i].partiality != partials[i].partiality
            || partials[i].partiality > 1 || partials[i].partiality < minPartiality)
            continue;
        
        if (partials[i].percentage != partials[i].percentage)
            continue;
        
        partialities.push_back(partials[i].partiality);
        percentages.push_back(partials[i].percentage);
        
     //   std::cout << partials[i].partiality << "\t" << partials[i].percentage << std::endl;
    }
    
    double correl = 0;
    
    if (typeOfScore == ScoreTypePartialityLeastSquares)
    {
        correl = 1
        - least_squares_between_vectors(&partialities, &percentages, 0);
    }
    else if (typeOfScore == ScoreTypePartialityGradient)
    {
        correl = gradient_between_vectors(&partialities, &percentages);
    }
    else
    {
        correl = correlation_through_origin(&partialities, &percentages);
    }
    
    return correl;
}

double MtzManager::minimizeTwoParameters(double *meanStep1, double *meanStep2,
                                         double **params, int paramNum1, int paramNum2,
                                         double (*score)(void *object, double lowRes, double highRes), void *object,
                                         double lowRes, double highRes, double low)
{
    double param_trials1[9];
    double param_trials2[9];
    double param_scores[9];
    
    int j = 0;
    double param_min_score = low;
    int param_min_num = 4;
    
    double bestParam1 = (*params)[paramNum1];
    double bestParam2 = (*params)[paramNum2];
    
    for (double i = bestParam1 - *meanStep1; j < 3; i += *meanStep1)
    {
        int l = 0;
        
        for (double k = bestParam2 - *meanStep2; l < 3; k += *meanStep2)
        {
            (*params)[paramNum1] = i;
            (*params)[paramNum2] = k;
            this->refreshPartialities((*params));
            param_scores[j * 3 + l] = (*score)(object, lowRes, highRes);
            param_trials1[j * 3 + l] = i;
            param_trials2[j * 3 + l] = k;
            l++;
        }
        j++;
    }
    
    for (int i = 0; i < 9; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }
    
    (*params)[paramNum1] = param_trials1[param_min_num];
    (*params)[paramNum2] = param_trials2[param_min_num];
    
    if (param_min_num == 4)
    {
        *meanStep1 /= 2;
        *meanStep2 /= 2;
    }
    
    //	std::cout << param_trials1[param_min_num] << " " << param_trials2[param_min_num] << std::endl;
    
    return param_min_score;
}

double MtzManager::minimizeParameter(double *meanStep, double **params,
                                     int paramNum, double (*score)(void *object, double lowRes, double highRes),
                                     void *object, double lowRes, double highRes)
{
    double param_trials[3];
    double param_scores[3];
    
    int j = 0;
    double param_min_score = FLT_MAX;
    int param_min_num = 1;
    
    double bestParam = (*params)[paramNum];
    
    for (double i = bestParam - *meanStep; j < 3; i += *meanStep)
    {
        (*params)[paramNum] = i;
        this->refreshPartialities((*params));
        param_scores[j] = (*score)(object, lowRes, highRes);
        param_trials[j] = i;
        j++;
    }
    
    for (int i = 0; i < 3; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }
    
    (*params)[paramNum] = param_trials[param_min_num];
    
    if (param_min_num == 1)
        *meanStep /= 2;
    
    this->refreshPartialities((*params));

    return param_min_score;
}

double MtzManager::minimize()
{
    return minimize(exclusionScoreWrapper, this);
}

double MtzManager::minimize(double (*score)(void *object, double lowRes, double highRes), void *object)
{
    double wavelength = 0;
    
    bool reinitialiseWavelength = FileParser::getKey("REINITIALISE_WAVELENGTH", false);
    double minResolution = FileParser::getKey("MIN_REFINED_RESOLUTION", 0.0);
    
    if (isUsingFixedWavelength())
    {
        wavelength = this->getWavelength();
        
        if (wavelength == 0)
        {
            std::cout << "Warning: using fixed wavelength but wavelength is 0"
            << std::endl;
        }
    }
    else if (isUsingFixedWavelength() && reinitialiseWavelength)
    {
        wavelength = FileParser::getKey("INITIAL_WAVELENGTH", 1.75);
    }
    else
    {
        wavelength = bestWavelength(0, 0, false);
        
        if (this->isRejected())
            return 0;
        
        if (finalised)
            wavelength = this->getWavelength();
    }
    
    
    bandwidth = this->getBandwidth();
    
    bool optimisedMean = !optimisingWavelength || (scoreType == ScoreTypeStandardDeviation);
    bool optimisedBandwidth = !optimisingBandwidth || (scoreType == ScoreTypeStandardDeviation);
    bool optimisedExponent = !optimisingExponent || (scoreType == ScoreTypeStandardDeviation);
    bool optimisedSpotSize = !optimisingRlpSize || (scoreType == ScoreTypeStandardDeviation);
    bool optimisedMos = !optimisingMosaicity || (scoreType == ScoreTypeStandardDeviation);
    bool optimisedHRot = !optimisingOrientation;
    bool optimisedKRot = !optimisingOrientation;
    
    double meanStep = stepSizeWavelength;
    double bandStep = stepSizeBandwidth;
    double expoStep = stepSizeExponent;
    double spotStep = stepSizeRlpSize;
    double mosStep = stepSizeMosaicity;
    double hStep = stepSizeOrientation;
    double kStep = stepSizeOrientation;
    double aStep = FileParser::getKey("STEP_SIZE_UNIT_CELL_A", 0.5);
    double bStep = FileParser::getKey("STEP_SIZE_UNIT_CELL_B", 0.5);
    double cStep = FileParser::getKey("STEP_SIZE_UNIT_CELL_C", 0.5);
  
    params = new double[PARAM_NUM];
    
    params[PARAM_HROT] = this->getHRot();
    params[PARAM_KROT] = this->getKRot();
    params[PARAM_MOS] = this->getMosaicity();
    params[PARAM_SPOT_SIZE] = this->getSpotSize();
    params[PARAM_WAVELENGTH] = wavelength;
    params[PARAM_BANDWIDTH] = bandwidth;
    params[PARAM_B_FACTOR] = bFactor;
    params[PARAM_SCALE_FACTOR] = scale;
    params[PARAM_EXPONENT] = this->getExponent();
    params[PARAM_UNIT_CELL_A] = this->cellDim[0];
    params[PARAM_UNIT_CELL_B] = this->cellDim[1];
    params[PARAM_UNIT_CELL_C] = this->cellDim[2];
    
    bool optimisingUnitCellA = FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
    bool optimisingUnitCellB = FileParser::getKey("OPTIMISING_UNIT_CELL_B", false);
    bool optimisingUnitCellC = FileParser::getKey("OPTIMISING_UNIT_CELL_C", false);
    
    bool optimisedUnitCellA = !FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
    bool optimisedUnitCellB = !FileParser::getKey("OPTIMISING_UNIT_CELL_B", false);
    bool optimisedUnitCellC = !FileParser::getKey("OPTIMISING_UNIT_CELL_C", false);
    
    int count = 0;
    
    refreshPartialities(params);
    
    int miniMethod = FileParser::getKey("MINIMIZATION_METHOD", 0);
    MinimizationMethod method = (MinimizationMethod)miniMethod;
    
    if (method == MinimizationMethodStepSearch)
    {
        while (!(optimisedMean && optimisedBandwidth && optimisedSpotSize
                 && optimisedMos && optimisedExponent && optimisedHRot
                 && optimisedKRot) && count < 50)
        {
            count++;
            
            if (scoreType != ScoreTypeStandardDeviation)
            {
                if (!optimisedMean)
                    minimizeParameter(&meanStep, &params, PARAM_WAVELENGTH, score,
                                      object, minResolution, maxResolutionAll);
                
                if (!optimisedBandwidth)
                    minimizeParameter(&bandStep, &params, PARAM_BANDWIDTH, score,
                                      object, minResolution, maxResolutionAll);
            }
            
            if (!optimisedHRot && !optimisedKRot)
                minimizeTwoParameters(&hStep, &kStep, &params, PARAM_HROT,
                                      PARAM_KROT, score, object, minResolution, maxResolutionAll, FLT_MAX);
            
            if (scoreType != ScoreTypeStandardDeviation)
            {
                if (!optimisedSpotSize)
                    minimizeParameter(&spotStep, &params, PARAM_SPOT_SIZE, score,
                                      object, minResolution, maxResolutionAll);
                
                if (!optimisedExponent)
                    minimizeParameter(&expoStep, &params, PARAM_EXPONENT, score, object,
                                      minResolution, maxResolutionAll);
                
                if (!optimisedMos)
                    minimizeParameter(&mosStep, &params, PARAM_MOS, score, object,
                                      minResolution, maxResolutionAll);
            }
            
            if (hStep < toleranceOrientation)
                optimisedHRot = true;
            
            if (kStep < toleranceOrientation)
                optimisedKRot = true;
            
            if (meanStep < toleranceWavelength)
                optimisedMean = true;
            
            if (bandStep < toleranceBandwidth)
                optimisedBandwidth = true;
            
            if (mosStep < toleranceMosaicity)
                optimisedMos = true;
            
            if (spotStep < toleranceRlpSize)
                optimisedSpotSize = true;
            
            if (expoStep < toleranceExponent)
                optimisedExponent = true;
            
            if (scoreType == ScoreTypeStandardDeviation)
            {
                std::cout << "stdev" << "\t" << params[PARAM_HROT] << "\t" << params[PARAM_KROT] << "\t" << (*score)(object, 0, 0) << std::endl;
            }
            
        }
        
        this->setParams(params);
        this->refreshPartialities(params);
        
        count = 0;
        
        while (!(optimisedUnitCellA && optimisedUnitCellB
                 && optimisedUnitCellC) && count < 50)
        {
            if (!optimisedUnitCellA)
            {
                minimizeParameter(&aStep, &params, PARAM_UNIT_CELL_A, score, object, 0, maxResolutionAll);
            }
            
            if (!optimisedUnitCellB)
            {
                minimizeParameter(&bStep, &params, PARAM_UNIT_CELL_B, score, object, 0, maxResolutionAll);
            }
            
            if (!optimisedUnitCellC)
            {
                minimizeParameter(&cStep, &params, PARAM_UNIT_CELL_C, score, object, 0, maxResolutionAll);
            }
            
            //    std::cout << (*score)(object, 0, maxResolutionAll) << ", " << correlation() << std::endl;
            
            if (aStep < 0.005)
                optimisedUnitCellA = true;
            
            if (bStep < 0.005)
                optimisedUnitCellB = true;
            
            if (cStep < 0.005)
                optimisedUnitCellC = true;
            
            count++;
        }
        
        this->refreshPartialities(params);
        
        bool refineB = FileParser::getKey("REFINE_B_FACTOR", false);
        
        if (refineB && bFactor == 0)
        {
            double bStep = 10;
            double optimisedB = false;
            int count = 0;
            
            while (!optimisedB && count < 30)
            {
                minimizeParam(bStep, bFactor, bFactorScoreWrapper, this);
                
                //    std::cout << bFactor << "\t" << score << std::endl;
                
                if (bStep < 0.01)
                    optimisedB = true;
                
                count++;
            }
        }
        
        this->setParams(params);
        this->refreshPartialities(params);
        this->applyBFactor(bFactor);
    }
    else if (method == MinimizationMethodNelderMead)
    {
        std::vector<double *>paramPtrs;
        paramPtrs.resize(PARAM_NUM);
        
        if (this->wavelength == 0)
            this->wavelength = wavelength;
        
        paramPtrs[PARAM_HROT] = optimisingOrientation ? &this->hRot : NULL;
        paramPtrs[PARAM_KROT] = optimisingOrientation ? &this->kRot : NULL;
        paramPtrs[PARAM_MOS] = optimisingMosaicity ? &this->mosaicity : NULL;
        paramPtrs[PARAM_SPOT_SIZE] = optimisingRlpSize ? &this->spotSize : NULL;
        paramPtrs[PARAM_WAVELENGTH] = optimisingWavelength ? &this->wavelength : NULL;
        paramPtrs[PARAM_BANDWIDTH] = optimisingBandwidth ? &this->bandwidth : NULL;
        paramPtrs[PARAM_B_FACTOR] = NULL;
        paramPtrs[PARAM_SCALE_FACTOR] = NULL;
        paramPtrs[PARAM_EXPONENT] = optimisingExponent ? &this->exponent : NULL;
        paramPtrs[PARAM_UNIT_CELL_A] = optimisingUnitCellA ? &this->cellDim[0] : NULL;
        paramPtrs[PARAM_UNIT_CELL_B] = optimisingUnitCellB ? &this->cellDim[1] : NULL;
        paramPtrs[PARAM_UNIT_CELL_C] = optimisingUnitCellC ? &this->cellDim[2] : NULL;
        
        std::vector<double> ranges;
        ranges.resize(PARAM_NUM);
        
        ranges[PARAM_HROT] = stepSizeOrientation;
        ranges[PARAM_KROT] = stepSizeOrientation;
        ranges[PARAM_MOS] = stepSizeMosaicity;
        ranges[PARAM_SPOT_SIZE] = stepSizeRlpSize;
        ranges[PARAM_WAVELENGTH] = stepSizeWavelength;
        ranges[PARAM_BANDWIDTH] = stepSizeBandwidth;
        ranges[PARAM_B_FACTOR] = 0;
        ranges[PARAM_SCALE_FACTOR] = 0;
        ranges[PARAM_EXPONENT] = stepSizeExponent;
        ranges[PARAM_UNIT_CELL_A] = FileParser::getKey("STEP_SIZE_UNIT_CELL_A", 0.5);
        ranges[PARAM_UNIT_CELL_B] = FileParser::getKey("STEP_SIZE_UNIT_CELL_B", 0.5);
        ranges[PARAM_UNIT_CELL_C] = FileParser::getKey("STEP_SIZE_UNIT_CELL_C", 0.5);
        
        NelderMead refiner(paramPtrs, ranges, this, &this->scoreNelderMead);
        refiner.process();
        double lastScore = rSplit(0, 0);
        
        getParams(&params);
        
    }

    double correl = correlation(true);
    
    this->wavelength = params[PARAM_WAVELENGTH];
    this->bandwidth = params[PARAM_BANDWIDTH];
    this->mosaicity = params[PARAM_MOS];
    this->spotSize = params[PARAM_SPOT_SIZE];
    this->hRot = params[PARAM_HROT];
    this->kRot = params[PARAM_KROT];
    this->exponent = params[PARAM_EXPONENT];
    this->cellDim[0] = params[PARAM_UNIT_CELL_A];
    this->cellDim[1] = params[PARAM_UNIT_CELL_B];
    this->cellDim[2] = params[PARAM_UNIT_CELL_C];
    this->refCorrelation = correl;
    
    delete[] params;
    
    if (scoreType == ScoreTypeSymmetry)
    {
        return rFactorWithManager(RFactorTypeMeas);
    }
    
    logged << "Returning correl: " << correl << std::endl;
    
    return correl;
}


void MtzManager::gridSearch(bool silent)
{
    scoreType = defaultScoreType;
    
    if (defaultScoreType == ScoreTypePartialityCorrelation)
    {
        if (!(refPartCorrel > 0.65 && trust == TrustLevelGood))
            scoreType = ScoreTypeMinimizeRSplit;
    }
    
    chooseAppropriateTarget();
    std::string scoreDescription = this->describeScoreType();
    penaltyWeight = FileParser::getKey("PENALTY_WEIGHT", 0.5);
    penaltyResolution = FileParser::getKey("PENALTY_RESOLUTION", 2.5);
    
    double scale = this->gradientAgainstManager(this->getReferenceManager());
    applyScaleFactor(scale);
    
    this->reallowPartialityOutliers();
    
    std::map<int, std::pair<vector<double>, double> > ambiguityResults;
    
    double *firstParams = new double[PARAM_NUM];
    getParams(&firstParams);
    
    if (trust == TrustLevelGood)
        minimize();
    
    if (trust != TrustLevelGood)
    {
        for (int i = 0; i < ambiguityCount(); i++)
        {
            vector<double> bestParams;
            bestParams.resize(PARAM_NUM);
            setParams(firstParams);
            
            int ambiguity = getActiveAmbiguity();
            double correl = minimize();
            double *params = &(*(bestParams.begin()));
            getParams(&params);
            
            std::pair<vector<double>, double> result = std::make_pair(bestParams, correl);
            
            ambiguityResults[ambiguity] = result;
            incrementActiveAmbiguity();
        }
        
        double bestScore = -1;
        int bestAmbiguity = 0;
        
        logged << "Ambiguity results: ";
        
        for (int i = 0; i < ambiguityCount(); i++)
        {
            std::pair<vector<double>, double> result = ambiguityResults[i];
            
            logged << result.second << " ";
            
            if (result.second > bestScore)
            {
                bestScore = result.second;
                bestAmbiguity = i;
            }
        }
        
        logged << std::endl;
        sendLog(LogLevelDetailed);
        
        setActiveAmbiguity(bestAmbiguity);
        setParams(&(*(ambiguityResults[bestAmbiguity].first).begin()));
        refreshPartialities(&(*(ambiguityResults[bestAmbiguity].first).begin()));
        
        double threshold = 0.9;
        
        if (bestScore > threshold)
        {
            trust = TrustLevelGood;
        }
    }
    
    double newCorrel = (scoreType == ScoreTypeSymmetry) ? 0 : correlation(true);
    this->setFinalised(true);
    
    if (newCorrel == 1 || newCorrel == -1)
        setFinalised(false);
    
    bool partialityRejection = FileParser::getKey("PARTIALITY_REJECTION", false);
    bool correlationRejection = FileParser::getKey("CORRELATION_REJECTION", true);
    
    if (partialityRejection && scoreType != ScoreTypeSymmetry)
        this->excludePartialityOutliers();
    
    if (correlationRejection && scoreType != ScoreTypeSymmetry)
        this->excludeFromLogCorrelation();
    double hits = accepted();
    
    double newerCorrel = (scoreType == ScoreTypeSymmetry) ? rFactorWithManager(RFactorTypeMeas) : correlation(true);
    double partCorrel =  (scoreType == ScoreTypeSymmetry) ? 0 : leastSquaresPartiality(ScoreTypePartialityCorrelation);
    double rSplitValue = (scoreType == ScoreTypeSymmetry) ? 0 : rSplit(0, maxResolutionAll);
    
    this->setRefPartCorrel(partCorrel);
    
    this->setRefCorrelation(newerCorrel);
    
    if (scoreType == ScoreTypeSymmetry)
    {
        this->setRefCorrelation(1 - newerCorrel);
        
        if (newerCorrel != newerCorrel)
            this->setRefCorrelation(0);
    }
    
    this->sendLog(LogLevelDetailed);
    
    Scaler *scaler = Scaler::getScaler();
    
    double rMerge = 0;
    
    if (scaler)
        rMerge = scaler->evaluateForImage(this);
    
    logged << filename << "\t" << scoreDescription << "\t" << "\t"
    << newerCorrel << "\t" << rSplitValue << "\t"
    << partCorrel << "\t" << rMerge << "\t" << hits << std::endl;
    
    this->sendLog(silent ? LogLevelDebug : LogLevelNormal);
    
    delete[] firstParams;
    
    writeToFile(std::string("ref-") + filename);
}

void MtzManager::excludeFromLogCorrelation()
{
    MtzManager &image1 = *this;
    MtzManager &image2 = *(MtzManager::referenceManager);
    
    double lowCut = 0;
    double highCut = 1 / maxResolutionAll;
    
    vector<double> refIntensities;
    vector<double> imageIntensities;
    vector<double> weights;
    vector<Reflection *> imgReflections;
    
    for (int i = 0; i < image1.reflectionCount(); i++)
    {
        Reflection *reflection = image1.reflection(i);
        Reflection *reflection2 = NULL;
        
        if (reflection->getResolution() < lowCut
            || reflection->getResolution() > highCut)
            continue;
        
        int refl = reflection->getReflId();
        
        image2.findReflectionWithId(refl, &reflection2);
        
        if (reflection2 == NULL)
            continue;
        
        double int1 = (reflection->meanIntensity());
        
        double int2 = (reflection2->meanIntensityWithExclusion(&filename));
        
        double weight = reflection->meanPartiality() / reflection2->meanSigma();
        
        if (int1 != int1 || int2 != int2 || weight != weight)
            continue;
        
        if (!std::isfinite(int1) || !std::isfinite(int2))
            continue;
        
        imgReflections.push_back(reflection);
        imageIntensities.push_back(int1);
        refIntensities.push_back(int2);
        weights.push_back(weight);
    }
    
    if (imageIntensities.size() < 100)
        return;
    
    std::map<double, int> correlationResults;
    
    double baseCorrelation = correlation_between_vectors(&refIntensities,
                                                         &imageIntensities, &weights);
    
    if (baseCorrelation > 0.99)
        return;
    
    double correlationToGain = 1 - baseCorrelation;
    
    for (int i = 0; i < refIntensities.size(); i++)
    {
        double newCorrelation = correlation_between_vectors(&refIntensities,
                                                            &imageIntensities, &weights, i);
        double newCorrelationToGain = 1 - newCorrelation;
        
        double ratio = (correlationToGain - newCorrelationToGain)
        / correlationToGain;
        
        correlationResults[ratio] = i;
    }
    
    int count = 0;
    
    for (std::map<double, int>::iterator it = correlationResults.end();
         it != correlationResults.begin(); --it)
    {
        if (count >= 3)
            break;
        
        if (it->first > 0.06)
        {
            imgReflections[correlationResults[it->first]]->miller(0)->setRejected(
                                                                              "correl", true);
            count++;
        }
    }
}

double MtzManager::partialityRatio(Reflection *imgReflection, Reflection *refReflection)
{
    double rawIntensity = imgReflection->miller(0)->getRawIntensity();
    double percentage = rawIntensity /= refReflection->meanIntensity();
    
    double partiality = imgReflection->meanPartiality();
    
    double ratio = percentage / partiality;
    
    if (ratio != ratio || !std::isfinite(ratio))
        return 10;
    
    return ratio;
}

void MtzManager::reallowPartialityOutliers()
{
    for (int i = 0; i < reflectionCount(); i++)
    {
        reflection(i)->miller(0)->setRejected("partiality", false);
    }
}

void MtzManager::excludePartialityOutliers()
{
    vector<Reflection *> refReflections, imgReflections;
    
    applyScaleFactor(this->gradientAgainstManager(referenceManager));
    
    this->findCommonReflections(referenceManager, imgReflections, refReflections,
                                NULL);
    
    vector<double> ratios;
    
    for (int i = 0; i < refReflections.size(); i++)
    {
        Reflection *imgReflection = imgReflections[i];
        Reflection *refReflection = refReflections[i];
        
        if (!imgReflection->anyAccepted())
            continue;
        
        double ratio = partialityRatio(imgReflection, refReflection);
        ratios.push_back(ratio);
    }
    
    double stdev = standard_deviation(&ratios);
    
    double upperBound = 1 + stdev * 1.4;
    double lowerBound = 1 - stdev * 0.7;
    
    if (lowerBound < 0.5)
        lowerBound = 0.5;
    if (upperBound > 1.7)
        upperBound = 1.7;
    
    int rejectedCount = 0;
    
    for (int i = 0; i < refReflections.size(); i++)
    {
        Reflection *imgReflection = imgReflections[i];
        Reflection *refReflection = refReflections[i];
        
        double ratio = partialityRatio(imgReflection, refReflection);
        
        if (ratio > upperBound || ratio < lowerBound)
        {
            if (!imgReflection->miller(0)->accepted())
                continue;
            
            if (imgReflection->betweenResolutions(1.7, 0))
                continue;
            
            imgReflection->miller(0)->setRejected("partiality", true);
            rejectedCount++;
        }
    }
    
    logged << "Rejected " << rejectedCount << " outside partiality range "
    << lowerBound * 100 << "% to " << upperBound * 100 << "%." << std::endl;
}

bool compareResult(ResultTuple one, ResultTuple two)
{
    return boost::get<4>(two) > boost::get<4>(one);
}

void MtzManager::findSteps(int param1, int param2, int param3)
{
    params = new double[PARAM_NUM];
    
    if (isnan(wavelength))
        return;
    
    params[PARAM_HROT] = this->getHRot();
    params[PARAM_KROT] = this->getKRot();
    params[PARAM_MOS] = this->getMosaicity();
    params[PARAM_SPOT_SIZE] = this->getSpotSize();
    params[PARAM_WAVELENGTH] = wavelength;
    params[PARAM_BANDWIDTH] = bandwidth;
    params[PARAM_EXPONENT] = this->getExponent();
    params[PARAM_B_FACTOR] = this->bFactor;
    params[PARAM_SCALE_FACTOR] = this->getScale();
    params[PARAM_UNIT_CELL_A] = this->cellDim[0];
    params[PARAM_UNIT_CELL_B] = this->cellDim[1];
    params[PARAM_UNIT_CELL_C] = this->cellDim[2];
    
    
    refreshPartialities(params);
    
    double *ranges = new double[3];
    double divisions = 9;
    
    for (int i = 0; i < 3; i++)
    {
        int chosen = (i == 0) ? param1 : ((i == 1) ? param2 : param3);
        
        if (chosen == PARAM_HROT || chosen == PARAM_KROT)
                ranges[i] = 0.08;
        else if (chosen == PARAM_WAVELENGTH)
            ranges[i] = 0.03;
        else if (chosen == PARAM_BANDWIDTH)
            ranges[i] = 0.0005;
        else if (chosen == PARAM_SPOT_SIZE)
            ranges[i] = 0.0006;
        else if (chosen == PARAM_MOS)
            ranges[i] = 0.05;
        else if (chosen == PARAM_EXPONENT)
            ranges[i] = 0.8;
        else if (chosen == PARAM_UNIT_CELL_A || chosen == PARAM_UNIT_CELL_B || chosen == PARAM_UNIT_CELL_C)
            ranges[i] = 0.5;
        else if (chosen == PARAM_SCALE_FACTOR || chosen == PARAM_B_FACTOR)
            ranges[i] = 2.0;
        else if (chosen == -1)
            ranges[i] = 0.0;
    }
    
    double iMinParam = 0;
    double iMaxParam = 0;
    double iStep = 1;
    
    double jMinParam = 0;
    double jMaxParam = 0;
    double jStep = 1;
    
    double kMinParam = 0;
    double kMaxParam = 0;
    double kStep = 1;
    
    if (param1 >= 0)
    {
        iMinParam = params[param1] - ranges[0];
        iMaxParam = params[param1] + ranges[0];
        iStep = (iMaxParam - iMinParam) / divisions;
    }
    
    if (param2 >= 0)
    {
        jMinParam = params[param2] - ranges[1];
        jMaxParam = params[param2] + ranges[1];
        jStep = (jMaxParam - jMinParam) / divisions;
    }
    
    if (param3 >= 0)
    {
        kMinParam = params[param3] - ranges[2];
        kMaxParam = params[param3] + ranges[2];
        kStep = (kMaxParam - kMinParam) / divisions;
    }
    
    vector<ResultTuple> results;
    double bestResult = rSplit(0, maxResolutionAll);
    
    for (double iParam = iMinParam; iParam <= iMaxParam; iParam += iStep)
    {
        for (double jParam = jMinParam; jParam <= jMaxParam; jParam += jStep)
        {
            for (double kParam = kMinParam; kParam <= kMaxParam; kParam += kStep)
            {
                if (ranges[0] > 0)
                    params[param1] = iParam;
                if (ranges[1] > 0)
                    params[param2] = jParam;
                if (ranges[2] > 0)
                    params[param3] = kParam;
                this->refreshPartialities(params);
                double score = rSplit(0, maxResolutionAll);
                double rSplitValue = rSplit(0, maxResolutionAll);
                
                ResultTuple result = boost::make_tuple<>(iParam, jParam, kParam, rSplitValue, score);
                results.push_back(result);
                
                if (bestResult > score)
                {
                    logged << "newmin\t" << iParam << "\t" << jParam << "\t" << kParam << "\t" << rSplitValue << std::endl;
                    sendLog(LogLevelDetailed);
                    bestResult = score;
                }
                //      }
            }
        }
    }
    
    if (results.size() == 0)
        return;
    
    std::sort(results.begin(), results.end(), compareResult);
    
    double paramA = (boost::get<0>(results[0]));
    double paramB = (boost::get<1>(results[0]));
    double paramC = (boost::get<2>(results[0]));
    
    if (ranges[0] > 0)
        params[param1] = paramA;
    if (ranges[1] > 0)
        params[param2] = paramB;
    if (ranges[2] > 0)
        params[param3] = paramC;

    setParams(params);
    double rSplit = boost::get<4>(results[0]);
    
    this->refreshPartialities(params);
    
    this->setRefCorrelation(correlation());
    this->setFinalised(true);
    
    delete [] ranges;
    
    this->writeToFile("ref-" + getFilename());
}

void MtzManager::chooseAppropriateTarget()
{
    if (defaultScoreType != ScoreTypeMinimizeRMeas)
        return;
    
    if (MtzRefiner::getCycleNum() == 0)
    {
        scoreType = ScoreTypeMinimizeRSplit;
        return;
    }
    
    if (failedCount == 0)
    {
        scoreType = ScoreTypeMinimizeRMeas;
    }
    
    if (failedCount > 1)
    {
        scoreType = ScoreTypeMinimizeRSplit;
        resetDefaultParameters();
        wavelength = bestWavelength();
        applyUnrefinedPartiality();
    }
}

void MtzManager::resetDefaultParameters()
{
    if (!setInitialValues)
    {
        wavelength = FileParser::getKey("INITIAL_WAVELENGTH", 0.0);
        bandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
        mosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);
        spotSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
        exponent = FileParser::getKey("INITIAL_EXPONENT", INITIAL_EXPONENT);
    }
    
    usingFixedWavelength = (wavelength != 0);
    hRot = 0;
    kRot = 0;
    allowTrust = FileParser::getKey("ALLOW_TRUST", true);
    bool alwaysTrust = FileParser::getKey("TRUST_INDEXING_SOLUTION", false);
    
    if (alwaysTrust)
        trust = TrustLevelGood;
}
