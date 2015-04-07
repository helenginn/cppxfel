#include "StatisticsManager.h"
#include <cmath>
#include "Vector.h"
#include <algorithm>

#include "FileParser.h"
#include "MtzManager.h"
#include "parameters.h"
#include "GraphDrawer.h"

#define DEFAULT_RESOLUTION 3.4
typedef boost::tuple<double, double, double, double> ResultTuple;

double MtzManager::weightedBestWavelength(double lowRes, double highRes)
{
    vector<double> wavelengths;
    vector<double> percentages;
    
    for (int i = 0; i < holderCount(); i++)
    {
        Holder *imageHolder = holder(i);
        
        if (imageHolder->getResolution() < lowRes)
            continue;
        
        if (imageHolder->getResolution() > highRes)
            continue;
        
        Holder *refHolder = NULL;
        int reflid = imageHolder->getReflId();
        
        referenceManager->findHolderWithId(reflid, &refHolder);
        
        if (refHolder != NULL)
        {
            double wavelength = imageHolder->miller(0)->getWavelength();
            
            double percentage = imageHolder->miller(0)->getRawIntensity()
            / imageHolder->miller(0)->getCountingSigma();
            
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
        return static_cast<MtzManager *>(object)->rSplit(lowRes, highRes);
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
    else
        return static_cast<MtzManager *>(object)->exclusionScore(lowRes,
                                                                 highRes, ScoreTypeCorrelation);
}

double MtzManager::rSplit(double low, double high, bool square)
{
    double scale = this->gradientAgainstManager(*referenceManager);
    
    applyScaleFactor(scale);
    
    double sum_numerator = 0;
    double sum_denominator = 0;
    int count = 0;
    double weights = 0;
    
    double lowCut = low == 0 ? 0 : 1 / low;
    double highCut = high == 0 ? FLT_MAX : 1 / high;
    
    vector<Holder *> holders1;
    vector<Holder *> holders2;
    
    this->findCommonReflections(referenceManager, holders1, holders2, NULL);
    
    for (int i = 0; i < holders1.size(); i++)
    {
        Holder *holder = holders1[i];
        Holder *holder2 = holders2[i];
        
        if (holder->getResolution() < lowCut
            || holder->getResolution() > highCut)
            continue;
        
        double int1 = holder->meanIntensity();
        double int2 = holder2->meanIntensityWithExclusion(&filename);
        double weight = holder->meanWeight() * holder->getResolution();
        
        if (int1 != int1 || int2 != int2)
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
    double score = leastSquaresPartiality(0, 1.5, typeOfScore);
    
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
    
    double lowCut = low == 0 ? 0 : 1 / low;
    double highCut = high == 0 ? FLT_MAX : 1 / high;
    
    for (int i = 0; i < holders.size(); i++)
    {
        Holder *imageHolder = holders[i];
        Holder *refHolder = NULL;
        int reflid = imageHolder->getReflId();
        
        referenceManager->findHolderWithId(reflid, &refHolder);
        
        if (refHolder != NULL)
        {
            Partial partial;
            partial.partiality = imageHolder->miller(0)->getPartiality();
            partial.percentage = imageHolder->miller(0)->getRawIntensity()
            / refHolder->meanIntensityWithExclusion(&filename);
            partial.resolution = imageHolder->getResolution();
            
            partials.push_back(partial);
        }
    }
    
    std::sort(partials.begin(), partials.end(), percentage);
    int position = (int)partials.size() - 5;
    if (position <= 0) return 0;
    double maxPercentage = partials[position].percentage;
    
    for (int i = 0; i < partials.size(); i++)
    {
        if (partials[i].percentage > maxPercentage)
            continue;
        
        if (partials[i].resolution > highCut || partials[i].resolution < lowCut)
            continue;
        
        if (partials[i].partiality != partials[i].partiality
            || partials[i].partiality > 1 || partials[i].partiality < 0.02)
            continue;
        
        if (partials[i].percentage != partials[i].percentage)
            continue;
        
        partialities.push_back(partials[i].partiality);
        percentages.push_back(partials[i].percentage);
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
    
    return param_min_score;
}

double MtzManager::minimize(bool minimizeSpotSize, bool suppress)
{
    return minimize(minimizeSpotSize, suppress, exclusionScoreWrapper, this);
}

double MtzManager::minimize(bool minimizeSpotSize, bool suppress,
                            double (*score)(void *object, double lowRes, double highRes), void *object)
{
    double wavelength = 0;
    
    if (isUsingFixedWavelength())
    {
        wavelength = this->getWavelength();
        
        if (wavelength == 0)
        {
            std::cout << "Warning: using fixed wavelength but wavelength is 0"
            << std::endl;
        }
    }
    else
    {
        wavelength = bestWavelength();
        
        if (this->isRejected())
            return 0;
        
        if (finalised)
            wavelength = this->getWavelength();
    }
    
    bandwidth = this->getBandwidth();
    
    bool optimisedMean = !optimisingWavelength;
    bool optimisedBandwidth = !optimisingBandwidth;
    bool optimisedExponent = !optimisingExponent;
    bool optimisedSpotSize = !optimisingRlpSize;
    bool optimisedMos = !optimisingMosaicity;
    bool optimisedHRot = !optimisingOrientation;
    bool optimisedKRot = !optimisingOrientation;
    
    double meanStep = stepSizeWavelength;
    double bandStep = stepSizeBandwidth;
    double expoStep = stepSizeExponent;
    double spotStep = stepSizeRlpSize;
    double mosStep = stepSizeMosaicity;
    double hStep = stepSizeOrientation;
    double kStep = stepSizeOrientation;
    
    bool optimisedUnitCell = false;
    double unitCellStepSize = 0.2;
    double unitCellTolerance = 0.001;
    
    params = new double[PARAM_NUM];
    
    params[PARAM_HROT] = this->getHRot();
    params[PARAM_KROT] = this->getKRot();
    params[PARAM_MOS] = this->getMosaicity();
    params[PARAM_SPOT_SIZE] = this->getSpotSize();
    params[PARAM_WAVELENGTH] = wavelength;
    params[PARAM_BANDWIDTH] = bandwidth;
    params[PARAM_EXPONENT] = this->getExponent();
    
    int count = 0;
    
    while (!(optimisedMean && optimisedBandwidth && optimisedSpotSize
             && optimisedMos && optimisedExponent && optimisedHRot
             && optimisedKRot) && count < 50)
    {
        count++;
        
        if (!optimisedMean)
            minimizeParameter(&meanStep, &params, PARAM_WAVELENGTH, score,
                              object, 0, maxResolutionAll);
        
        if (!optimisedBandwidth)
            minimizeParameter(&bandStep, &params, PARAM_BANDWIDTH, score,
                              object, 0, maxResolutionAll);
        
        if (!optimisedHRot && !optimisedKRot)
            minimizeTwoParameters(&hStep, &kStep, &params, PARAM_HROT,
                                  PARAM_KROT, score, object, maxResolutionRlpSize, maxResolutionAll, FLT_MAX);
        
        if (!optimisedSpotSize)
            minimizeParameter(&spotStep, &params, PARAM_SPOT_SIZE, score,
                              object, 0, maxResolutionRlpSize);
        
        if (!optimisedExponent)
            minimizeParameter(&expoStep, &params, PARAM_EXPONENT, score, object,
                              0, maxResolutionAll);
        
        if (!optimisedMos)
            minimizeParameter(&mosStep, &params, PARAM_MOS, score, object,
                              0, maxResolutionAll);
        
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
    }
    
    bool refineB = FileParser::getKey("REFINE_B_FACTOR", false);
    
    if (refineB && bFactor == 0)
    {
        double bStep = 10;
        double optimisedB = false;
        int count = 0;
        
        while (!optimisedB && count < 30)
        {
            double score = minimizeParam(bStep, bFactor, bFactorScoreWrapper, this);
            
            //    std::cout << bFactor << "\t" << score << std::endl;
            
            if (bStep < 0.01)
                optimisedB = true;
            
            count++;
        }
    }
    
    this->refreshPartialities(params);
    this->applyBFactor(bFactor);
    
    double newScore = (*score)(object, 0, 0);
    double hits = accepted();
    string scoreDescription = this->describeScoreType();
    double correl = correlation(true);
    
    this->wavelength = params[PARAM_WAVELENGTH];
    this->bandwidth = params[PARAM_BANDWIDTH];
    this->mosaicity = params[PARAM_MOS];
    this->spotSize = params[PARAM_SPOT_SIZE];
    this->hRot = params[PARAM_HROT];
    this->kRot = params[PARAM_KROT];
    this->exponent = params[PARAM_EXPONENT];
    this->refCorrelation = correl;
    
    delete[] params;
    
    if (scoreType == ScoreTypeSymmetry)
    {
        return rFactorWithManager(RFactorTypeMeas);
    }
    
    return correl;
}

void MtzManager::gridSearchWrapper(MtzManager *image, bool minimizeSpotSize)
{
    image->gridSearch(minimizeSpotSize);
}

void MtzManager::gridSearch(
                            double (*score)(void *object, double lowRes, double highRes), void *object)
{
    minimize(true, 0, score, object);
}

void MtzManager::gridSearch(bool minimizeSpotSize)
{
    scoreType = defaultScoreType;
    string scoreDescription = this->describeScoreType();
    
    double scale = this->gradientAgainstManager(*referenceManager);
    applyScaleFactor(scale);
    
    this->reallowPartialityOutliers();
    
    std::map<int, std::pair<std::vector<double>, double> > ambiguityResults;
    
    double *firstParams = new double[PARAM_NUM];
    getParams(&firstParams);
    
    if (trust == TrustLevelGood)
        minimize(minimizeSpotSize, true);
    
    if (trust != TrustLevelGood)
    {
        for (int i = 0; i < ambiguityCount(); i++)
        {
            vector<double> bestParams;
            bestParams.resize(PARAM_NUM + 5);
            setParams(firstParams);
            
            int ambiguity = getActiveAmbiguity();
            double correl = minimize();
            double *params = &(*(bestParams.begin()));
            getParams(&params);
            
            int hits = 0;
            double realCorrel = StatisticsManager::cc_pearson(this, getReferenceManager(), true, &hits, NULL, 0, 0, false);
            
            std::pair<std::vector<double>, double> result = std::make_pair(bestParams, correl);
            
            ambiguityResults[ambiguity] = result;
            incrementActiveAmbiguity();
        }
        
        double bestScore = -1;
        int bestAmbiguity = 0;
        
        logged << "Ambiguity results: ";
        
        for (int i = 0; i < ambiguityCount(); i++)
        {
            std::pair<std::vector<double>, double> result = ambiguityResults[i];
            
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
        
        double threshold = 0.9;
        
        if (bestScore > threshold)
        {
            trust = TrustLevelGood;
        }
    }
    /*
     if (trust != TrustLevelGood)
     {
     if (isRLog)
     scoreType = ScoreTypeCorrelation;
     double correl = minimize(minimizeSpotSize, true);
     getParams(&bestParams);
     setParams(firstParams);
     
     setInverse(!inverse);
     double invCorrel = minimize(minimizeSpotSize, true);
     
     if (correl > invCorrel)
     {
     setInverse(!inverse);
     setParams(bestParams);
     refreshPartialities(bestParams);
     finalised = true;
     logged << "Not flipping image" << std::endl;
     }
     else
     {
     logged << "Flipping image" << std::endl;
     }
     
     if (scoreType != ScoreTypeSymmetry)
     {
     double startCorrel = correlation(true);
     double threshold = 0.9;
     
     if (startCorrel > threshold && allowTrust)
     {
     logged << "Correlation > " << threshold << ", trust level = good" << std::endl;
     trust = TrustLevelGood;
     }
     }
     }
     
     getParams(&bestParams);
     */
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
    
    this->setRefCorrelation(newerCorrel);
    
    if (scoreType == ScoreTypeSymmetry)
    {
        this->setRefCorrelation(1 - newerCorrel);
        
        if (newerCorrel != newerCorrel)
            this->setRefCorrelation(0);
    }
    
    this->sendLog(LogLevelDetailed);
    
    logged << filename << "\t" << scoreDescription << "\t" << "\t"
    << newerCorrel << "\t" << rSplitValue << "\t"
    << partCorrel << "\t" << abs(bFactor) << "\t" << hits << std::endl;
    
    this->sendLog(LogLevelNormal);
    
    delete[] firstParams;
    
    writeToFile(string("ref-") + filename);
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
    vector<Holder *> imgHolders;
    
    for (int i = 0; i < image1.holderCount(); i++)
    {
        Holder *holder = image1.holder(i);
        Holder *holder2 = NULL;
        
        if (holder->getResolution() < lowCut
            || holder->getResolution() > highCut)
            continue;
        
        int refl = holder->getReflId();
        
        image2.findHolderWithId(refl, &holder2);
        
        if (holder2 == NULL)
            continue;
        
        double int1 = (holder->meanIntensity());
        
        double int2 = (holder2->meanIntensityWithExclusion(&filename));
        
        double weight = holder->meanPartiality() / holder2->meanSigma();
        
        if (int1 != int1 || int2 != int2 || weight != weight)
            continue;
        
        if (!isfinite(int1) || !isfinite(int2))
            continue;
        
        imgHolders.push_back(holder);
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
            imgHolders[correlationResults[it->first]]->miller(0)->setRejected(
                                                                              "correl", true);
            count++;
        }
    }
}

double MtzManager::partialityRatio(Holder *imgHolder, Holder *refHolder)
{
    double rawIntensity = imgHolder->miller(0)->getRawIntensity();
    double percentage = rawIntensity /= refHolder->meanIntensity();
    
    double partiality = imgHolder->meanPartiality();
    
    double ratio = percentage / partiality;
    
    if (ratio != ratio || !isfinite(ratio))
        return 10;
    
    return ratio;
}

void MtzManager::reallowPartialityOutliers()
{
    for (int i = 0; i < holderCount(); i++)
    {
        holder(i)->miller(0)->setRejected("partiality", false);
    }
}

void MtzManager::excludePartialityOutliers()
{
    vector<Holder *> refHolders, imgHolders;
    
    applyScaleFactor(this->gradientAgainstManager(*referenceManager));
    
    this->findCommonReflections(referenceManager, imgHolders, refHolders,
                                NULL);
    
    vector<double> ratios;
    
    for (int i = 0; i < refHolders.size(); i++)
    {
        Holder *imgHolder = imgHolders[i];
        Holder *refHolder = refHolders[i];
        
        if (!imgHolder->anyAccepted())
            continue;
        
        double ratio = partialityRatio(imgHolder, refHolder);
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
    
    for (int i = 0; i < refHolders.size(); i++)
    {
        Holder *imgHolder = imgHolders[i];
        Holder *refHolder = refHolders[i];
        
        double ratio = partialityRatio(imgHolder, refHolder);
        
        if (ratio > upperBound || ratio < lowerBound)
        {
            if (!imgHolder->miller(0)->accepted())
                continue;
            
            if (imgHolder->betweenResolutions(1.7, 0))
                continue;
            
            imgHolder->miller(0)->setRejected("partiality", true);
            rejectedCount++;
        }
    }
    
    logged << "Rejected " << rejectedCount << " outside partiality range "
    << lowerBound * 100 << "% to " << upperBound * 100 << "%." << std::endl;
}

bool compareResult(ResultTuple one, ResultTuple two)
{
    return boost::get<3>(two) > boost::get<3>(one);
}

void MtzManager::findSteps()
{
    params = new double[PARAM_NUM];
    
    wavelength = bestWavelength();
    
    params[PARAM_HROT] = this->getHRot();
    params[PARAM_KROT] = this->getKRot();
    params[PARAM_MOS] = this->getMosaicity();
    params[PARAM_SPOT_SIZE] = this->getSpotSize();
    params[PARAM_WAVELENGTH] = wavelength;
    params[PARAM_BANDWIDTH] = bandwidth;
    params[PARAM_EXPONENT] = this->getExponent();
    
    refreshPartialities(params);
    
    double iMinParam = -0.04;
    double iMaxParam = 0.04;
    double iStep = (iMaxParam - iMinParam) / 10;
    
    double jMinParam = -0.04;
    double jMaxParam = 0.04;
    double jStep = (jMaxParam - jMinParam) / 10;
    
    double kMinParam = wavelength - 0.005;
    double kMaxParam = wavelength + 0.005;
    double kStep = (kMaxParam - kMinParam) / 20;
    
    std::vector<ResultTuple> results;
    
    for (double iParam = iMinParam; iParam < iMaxParam; iParam += iStep)
    {
        for (double jParam = jMinParam; jParam < jMaxParam; jParam += jStep)
        {
            for (double kParam = kMinParam; kParam < kMaxParam; kParam += kStep)
            {
                params[PARAM_HROT] = iParam;
                params[PARAM_KROT] = jParam;
                params[PARAM_WAVELENGTH] = kParam;
                this->refreshPartialities(params);
                double score = rSplit(0, 0);
                
                ResultTuple result = boost::make_tuple<>(iParam, jParam, kParam, score);
                results.push_back(result);
            }
        }
    }
    
    std::sort(results.begin(), results.end(), compareResult);
    
    double newHRot = (boost::get<0>(results[0]));
    double newKRot = (boost::get<1>(results[0]));
    double newWavelength = (boost::get<2>(results[0]));
    
    params[PARAM_HROT] = newHRot;
    params[PARAM_KROT] = newKRot;
    params[PARAM_WAVELENGTH] = newWavelength;
    
    setHRot(newHRot);
    setKRot(newKRot);
    setWavelength(newWavelength);
    
    double rSplit = boost::get<3>(results[0]);
    
    std::cout << getFilename() << "\t" << hRot << "\t" << kRot << "\t" << wavelength << "\t" << rSplit << std::endl;

    this->refreshPartialities(params);
    
    this->setRefCorrelation(correlation());

    this->writeToFile("ref-" + getFilename());
}
