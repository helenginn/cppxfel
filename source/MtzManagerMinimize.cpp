#include "StatisticsManager.h"
#include <cmath>
#include "Vector.h"
#include <algorithm>
#include <fstream>

#include "MtzRefiner.h"
#include "FileParser.h"
#include "MtzManager.h"
#include "parameters.h"
#include "GraphDrawer.h"
#include "Reflection.h"
#include "Miller.h"

#include "NelderMead.h"
#include "RefinementStepSearch.h"

#define DEFAULT_RESOLUTION 3.4
typedef boost::tuple<double, double, double, double, double> ResultTuple;

double MtzManager::scoreNelderMead(void *object)
{
    static_cast<MtzManager *>(object)->refreshCurrentPartialities();
    return static_cast<MtzManager *>(object)->exclusionScoreWrapper(object, 0, 0);
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
    MtzManager *mtz = static_cast<MtzManager *>(object);
    
    if (scoreType == ScoreTypeCorrelation
        || scoreType == ScoreTypeCorrelationLog)
    {
        return mtz->exclusionScore(lowRes, highRes,  scoreType);
    }
    else if (scoreType == ScoreTypeMinimizeRSplit)
    {
        return mtz->rSplit(lowRes, highRes);
    }
    else if (scoreType == ScoreTypeRSplitIntensity)
    {
        double rSplit = mtz->rSplit(lowRes, highRes);
        double logIntensity = log(mtz->averageIntensity());
        
        return rSplit / logIntensity;
    }
    else if (scoreType == ScoreTypeMinimizeRSplitLog)
    {
        return mtz->rSplit(lowRes, highRes);
    }
    else if (scoreType == ScoreTypeSymmetry)
    {
        return mtz->rFactorWithManager(RFactorTypeMeas);
    }
    else if (scoreType == ScoreTypePartialityCorrelation
             || scoreType == ScoreTypePartialityLeastSquares)
    {
        double value = mtz->partialityFunction();
        
        return  value;
    }
    else if (scoreType == ScoreTypeMinimizeRMeas)
    {
        double rSplit = mtz->rSplit(lowRes, highRes);
        return rSplit;
    }
    else
    {
        return mtz->exclusionScore(lowRes, highRes, ScoreTypeCorrelation);
    }
}

bool partialGreaterThanPartial(Partial a, Partial b)
{
    return (a.wavelength > b.wavelength);
}

double MtzManager::rSplit(double low, double high)
{
    bool reverse = FileParser::getKey("SMOOTH_FUNCTION", false);
    
    if (referenceManager == NULL)
    {
        return 0;
    }
    
    double scale = this->gradientAgainstManager(referenceManager);
    applyScaleFactor(scale);
    
    double sum_numerator = 0;
    double sum_denominator = 0;
    int count = 0;
    double weights = 0;
    
    vector<ReflectionPtr> reflections1;
    vector<ReflectionPtr> reflections2;
    
    this->findCommonReflections(referenceManager, reflections1, reflections2, NULL, true);
    
    for (int i = 0; i < reflections1.size(); i++)
    {
        ReflectionPtr reflection = reflections1[i];
        ReflectionPtr reflection2 = reflections2[i];
        
        if (!reverse && reflection->acceptedCount() == 0)
            continue;
        
        if (reflection2->millerCount() == 0)
            continue;
        
        if (reflection2->miller(0)->isFree())
            continue;
        
        if (!reflection->betweenResolutions(low, high))
            continue;
        
        for (int j = 0; j < reflection->millerCount(); j++)
        {
            double int1 = 0;
            double int2 = 0;
            double weight = 0;
            
            if (!reverse)
            {
                int1 = reflection->meanIntensity();
                int2 = reflection2->meanIntensityWithExclusion(&filename);
                weight = reflection->meanPartiality();
            }
            else
            {
                int1 = reflection->miller(j)->getRawIntensity();
                int2 = reflection2->meanIntensity() * reflection->miller(j)->getPartiality();
                weight = 1;
            }
            
            if (int1 == 0 || weight == 0 || weight != weight)
            {
                continue;
                
            }
            
            if (int1 != int1 || int2 != int2)
                continue;
            
            if (int1 + int2 < 0)
                continue;
            
            count++;
            weights += weight;
            
            sum_numerator += fabs(int1 - int2) * weight;
            sum_denominator += (int1 + int2) * weight / 2;
        }
        

    }
    
  //  if (!withCutoff)
  //      std::cout << count << std::endl;
    
    double r_split = sum_numerator / (sum_denominator * sqrt(2));
    lastRSplit = r_split;
    
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
    vector<double> intensities;
    
    for (int i = 0; i < reflections.size(); i++)
    {
        ReflectionPtr imageReflection = reflections[i];
        ReflectionPtr refReflection;
        int reflid = (int)imageReflection->getReflId();
        
        referenceManager->findReflectionWithId(reflid, &refReflection);
        
        if (refReflection)
        {
            if (refReflection->meanIntensity() < REFERENCE_WEAK_REFLECTION)
                continue;
            
            if (!refReflection->betweenResolutions(low, high))
                continue;
            
            Partial partial;
            partial.miller = imageReflection->miller(0);
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
        intensities.push_back(partials[i].miller->getRawestIntensity());
        
     //   std::cout << partials[i].partiality << "\t" << partials[i].percentage << std::endl;
    }
    
    double correl = 0;
    
    correl = correlation_through_origin(&partialities, &percentages, &intensities);
    
    return correl;
}

double MtzManager::partialityFunction(double low, double high)
{
    std::vector<double> calcPartials, obsPartials, rawIntensities;
    double calcPartialSum = 0;
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        ReflectionPtr reflection = reflections[i];
        if (!reflection->betweenResolutions(low, high))
            continue;
        
        for (int j = 0; j < reflection->millerCount(); j++)
        {
            MillerPtr miller = reflection->miller(j);
            
            double calcPartial = miller->getPartiality();
            double obsPartial = miller->observedPartiality(MtzManager::getReferenceManager());
            double rawIntensity = miller->getRawIntensity();
            
            if (rawIntensity != rawIntensity || obsPartial != obsPartial || calcPartial != calcPartial)
                continue;
            
            calcPartialSum += calcPartial;
            
            calcPartials.push_back(calcPartial);
            obsPartials.push_back(obsPartial);
            rawIntensities.push_back(rawIntensity);
        }
    }
    
    double score = 0;

    for (int i = 0; i < calcPartials.size(); i++)
    {
        if (obsPartials[i] < 0 || calcPartials[i] <= 0)
            continue;
        
        double addition = obsPartials[i] * calcPartials[i] * rawIntensities[i];
        addition /= calcPartialSum;
        score -= addition;
    }
    
    return score;
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

double MtzManager::refineParameterScore(void *object)
{
    MtzManager *me = static_cast<MtzManager *>(object);
    
    me->refreshCurrentPartialities();
    return exclusionScoreWrapper(me);
    
    
    me->updateLatestMatrix();
    
    for (int i = 0; i < me->reflectionCount(); i++)
    {
        for (int j = 0; j < me->reflection(i)->millerCount(); j++)
        {
            MillerPtr miller = me->reflection(i)->miller(j);
            miller->recalculateBetterPartiality();
        }
    }
    
    me->setScoreType(me->defaultScoreType);
    
    return exclusionScoreWrapper(me);
}

double MtzManager::minimize()
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
    
    RefinementStrategyPtr refiner = RefinementStrategy::userChosenStrategy();
    
    refiner->setJobName("Refining " + getFilename());
    
    if (this->wavelength == 0)
        this->wavelength = wavelength;
    
    double stepSizeUnitCellA = FileParser::getKey("STEP_SIZE_UNIT_CELL_A", 0.5);
    double stepSizeUnitCellB = FileParser::getKey("STEP_SIZE_UNIT_CELL_B", 0.5);
    double stepSizeUnitCellC = FileParser::getKey("STEP_SIZE_UNIT_CELL_C", 0.5);
    
    bool optimisingUnitCellA = FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
    bool optimisingUnitCellB = FileParser::getKey("OPTIMISING_UNIT_CELL_B", false);
    bool optimisingUnitCellC = FileParser::getKey("OPTIMISING_UNIT_CELL_C", false);
    
    int verbosity = FileParser::getKey("VERBOSITY_LEVEL", 0);
    
    refiner->setEvaluationFunction(scoreNelderMead, this);
    refiner->setVerbose((verbosity > 0));
    
    if (optimisingWavelength)
    {
        refiner->addParameter(this, getWavelengthStatic, setWavelengthStatic, stepSizeWavelength, 0, "wavelength");
    }
    
    if (optimisingOrientation)
    {
        refiner->addParameter(this, getHRotStatic, setHRotStatic, stepSizeOrientation, 0, "hRot");
        refiner->addCoupledParameter(this, getKRotStatic, setKRotStatic, stepSizeOrientation, 0, "kRot");
    }
    
    if (optimisingMosaicity)
    {
        refiner->addParameter(this, getMosaicityStatic, setMosaicityStatic, stepSizeMosaicity, 0, "mosaicity");
    }
    
    if (optimisingRlpSize)
    {
        refiner->addParameter(this, getSpotSizeStatic, setSpotSizeStatic, stepSizeRlpSize, 0, "rlpSize");
    }
    
    if (optimisingBandwidth)
    {
        refiner->addParameter(this, getBandwidthStatic, setBandwidthStatic, stepSizeBandwidth, 0, "bandwidth");
    }
    
    if (optimisingExponent)
    {
        refiner->addParameter(this, getExponentStatic, setExponentStatic, stepSizeExponent, 0, "exponent");
    }
    
    if (optimisingUnitCellA)
    {
        refiner->addParameter(this, getUnitCellAStatic, setUnitCellAStatic, stepSizeUnitCellA, 0, "unitCellA");
    }
    
    if (optimisingUnitCellB)
    {
        refiner->addParameter(this, getUnitCellBStatic, setUnitCellBStatic, stepSizeUnitCellB, 0, "unitCellB");
    }
    
    if (optimisingUnitCellC)
    {
        refiner->addParameter(this, getUnitCellCStatic, setUnitCellCStatic, stepSizeUnitCellC, 0, "unitCellC");
    }
    
    refiner->refine();
    
    double correl = correlation(true, minResolution, maxResolutionAll);
    this->refCorrelation = correl;
    
    logged << "Returning correl: " << correl << std::endl;
    
    return correl;
}


void MtzManager::gridSearch(bool silent)
{
    scoreType = defaultScoreType;
    chooseAppropriateTarget();
    std::string scoreDescription = this->describeScoreType();
    
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
            double correl = 0;
            
            correl = minimize();
            
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
    
    scale = this->gradientAgainstManager(this->getReferenceManager());
    applyScaleFactor(scale);
    
    double newCorrel = (scoreType == ScoreTypeSymmetry) ? 0 : correlation(true);
    this->setFinalised(true);
    
    if (newCorrel == 1 || newCorrel == -1)
        setFinalised(false);
    
    bool partialityRejection = FileParser::getKey("PARTIALITY_REJECTION", false);
    
    if (partialityRejection && scoreType != ScoreTypeSymmetry)
        this->excludePartialityOutliers();
    
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
    double rMerge = 0;
    
    logged << filename << "\t" << scoreDescription << "\t" << "\t"
    << newerCorrel << "\t" << rSplitValue << "\t"
    << partCorrel << "\t" << bFactor << "\t" << hits << std::endl;
    
    this->sendLog(silent ? LogLevelDetailed : LogLevelNormal);
    
    delete[] firstParams;
    
    writeToFile(std::string("ref-") + filename);
    writeToDat("ref-");
}

void MtzManager::excludeFromLogCorrelation()
{
    bool correlationRejection = FileParser::getKey("CORRELATION_REJECTION", true);
    
    if (!correlationRejection)
    {
        return;
    }
    
    MtzManager &image1 = *this;
    MtzManager &image2 = *(MtzManager::referenceManager);
    
    double lowCut = 0;
    double highCut = 1 / maxResolutionAll;
    
    vector<double> refIntensities;
    vector<double> imageIntensities;
    vector<double> weights;
    vector<ReflectionPtr> imgReflections;
    
    for (int i = 0; i < image1.reflectionCount(); i++)
    {
        ReflectionPtr reflection = image1.reflection(i);
        ReflectionPtr reflection2;
        
        if (reflection->getResolution() < lowCut
            || reflection->getResolution() > highCut)
            continue;
        
        int refl = (int)reflection->getReflId();
        
        image2.findReflectionWithId(refl, &reflection2);
        
        if (!reflection2)
            continue;
        
        double int1 = (reflection->meanIntensity());
        
        double int2 = (reflection2->meanIntensity());
        
        double weight = reflection->meanWeight();
        
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
                                                                                  RejectReasonCorrelation, true);
            count++;
        }
    }
}

double MtzManager::partialityRatio(ReflectionPtr imgReflection, ReflectionPtr refReflection)
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
        reflection(i)->miller(0)->setRejected(RejectReasonPartiality, false);
    }
}

void MtzManager::excludePartialityOutliers()
{
    vector<ReflectionPtr> refReflections, imgReflections;
    
    applyScaleFactor(this->gradientAgainstManager(referenceManager));
    
    this->findCommonReflections(referenceManager, imgReflections, refReflections,
                                NULL);
    
    vector<double> ratios;
    
    for (int i = 0; i < refReflections.size(); i++)
    {
        ReflectionPtr imgReflection = imgReflections[i];
        ReflectionPtr refReflection = refReflections[i];
        
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
        ReflectionPtr imgReflection = imgReflections[i];
        ReflectionPtr refReflection = refReflections[i];
        
        double ratio = partialityRatio(imgReflection, refReflection);
        
        if (ratio > upperBound || ratio < lowerBound)
        {
            if (!imgReflection->miller(0)->accepted())
                continue;
            
            if (imgReflection->betweenResolutions(1.7, 0))
                continue;
            
            imgReflection->miller(0)->setRejected(RejectReasonPartiality, true);
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

void MtzManager::findSteps(int param1, int param2, std::string csvName)
{
    std::ofstream csvStream;
    std::cout << "Writing to CSV stream " << csvName << std::endl;
    csvStream.open(csvName);
    
    params = new double[PARAM_NUM];
    double *goodParams = new double[PARAM_NUM];
    
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
    double divisions = FileParser::getKey("LANDSCAPE_DIVISIONS", 100);
    
    for (int i = 0; i < 2; i++)
    {
        int chosen = (i == 0) ? param1 : param2;
        
        if (chosen == PARAM_HROT || chosen == PARAM_KROT)
            ranges[i] = 0.2;
        else if (chosen == PARAM_WAVELENGTH)
            ranges[i] = 0.05;
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
    
    vector<ResultTuple> results;
    double bestResult = rSplit(0, maxResolutionAll);
    double lowestResult = exclusionScoreWrapper(this);
    
    for (int i = 0; i < divisions; i++)
    {
        double iParam = iMinParam + iStep * (double)i;
        
        for (int j = 0; j < divisions; j++)
        {
            double jParam = jMinParam + jStep * (double)j;
            
            if (ranges[0] > 0)
                params[param1] = iParam;
            if (ranges[1] > 0)
                params[param2] = jParam;
            
            this->refreshPartialities(params);
            double score = exclusionScoreWrapper(this);
            double rSplitValue = rSplit(0, maxResolutionAll);
            
            if (score < lowestResult)
            {
                lowestResult = score;
                getParams(&goodParams);
            }
            
            ResultTuple result = boost::make_tuple<>(iParam, jParam, 0, rSplitValue, score);
            results.push_back(result);
            
            csvStream << iParam << "," << jParam << "," << score << std::endl;
            sendLog(LogLevelNormal);
            
            if (bestResult > score)
            {
                //     logged << "newmin\t" << iParam << "\t" << jParam << "\t" << kParam << "\t" << rSplitValue << std::endl;
                //     sendLog(LogLevelDetailed);
                bestResult = score;
            }
        }
    }
    
    logged << "New lowest result: " << lowestResult << std::endl;
    sendLog();
    
    refreshPartialities(goodParams);
    writeToFile(std::string("ref-") + filename);
    
    delete [] params;
    delete [] goodParams;
    
    csvStream.close();
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
}

void MtzManager::resetDefaultParameters()
{
    if (!setInitialValues)
    {
        wavelength = FileParser::getKey("INITIAL_WAVELENGTH", wavelength);
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
