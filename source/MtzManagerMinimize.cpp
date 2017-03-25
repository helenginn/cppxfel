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

typedef boost::tuple<double, double, double, double, double> ResultTuple;

double MtzManager::exclusionScoreWrapper(void *object, double lowRes,
                                         double highRes)
{
    ScoreType scoreType = static_cast<MtzManager *>(object)->scoreType;
    MtzManager *mtz = static_cast<MtzManager *>(object);
    
    if (scoreType == ScoreTypeCorrelation)
    {
		return 1 - mtz->correlation(true, lowRes, highRes);
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
    else if (scoreType == ScoreTypeMinimizeRMeas)
    {
        double rSplit = mtz->rSplit(lowRes, highRes);
        return rSplit;
    }
    else
    {
        return mtz->rSplit(lowRes, highRes);
    }
}

double MtzManager::rSplit(double low, double high)
{
    bool reverse = FileParser::getKey("SMOOTH_FUNCTION", false);
    
    if (referenceManager == NULL)
    {
        return 0;
    }
    
    this->scaleToMtz(referenceManager);

    double sum_numerator = 0;
    double sum_denominator = 0;
    int count = 0;
    double weights = 0;
    
    vector<ReflectionPtr> referenceRefs;
    vector<ReflectionPtr> imageRefs;
    
    this->findCommonReflections(referenceManager, imageRefs, referenceRefs, NULL, true, true);
    
    for (int i = 0; i < referenceRefs.size(); i++)
    {
        ReflectionPtr referenceRef = referenceRefs[i];
        ReflectionPtr imageRef = imageRefs[i];
        
        if (!reverse && imageRef->acceptedCount() == 0)
            continue;
        
        if (referenceRef->millerCount() == 0)
            continue;
        
        if (imageRef->miller(0)->isFree())
            continue;
        
        if (!referenceRef->betweenResolutions(low, high))
            continue;
        
        for (int j = 0; j < imageRef->millerCount(); j++)
        {
            if (!imageRef->miller(j)->accepted())
            {
                continue;
            }
            
            double int1 = 0;
            double int2 = 0;
            double weight = 0;
            
            if (!reverse)
            {
                int1 =  imageRef->miller(j)->intensity();
                int2 = referenceRef->meanIntensity();
                weight = imageRef->meanPartiality();
            }
            else
            {
                int1 = imageRef->miller(j)->getRawIntensity();
                int2 = referenceRef->meanIntensity() * imageRef->miller(j)->getPartiality();
                weight = imageRef->meanPartiality();
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
    
    double r_split = sum_numerator / (sum_denominator * sqrt(2));
    lastRSplit = r_split;
    
    return r_split;
    
}

double MtzManager::leastSquaresPartiality(double low, double high)
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
            / refReflection->meanIntensity();
            partial.resolution = imageReflection->getResolution();
            
            partials.push_back(partial);
        }
    }
    
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
    }
    
    double correl = 0;
    
    correl = correlation_through_origin(&partialities, &percentages, &intensities);
    
    return correl;
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
      //      miller->recalculatePartiality();
            // fix me!
        }
    }
    
    me->setScoreType(me->defaultScoreType);
    
    return exclusionScoreWrapper(me);
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
        
        if (it->first > 0.15)
        {
            imgReflections[correlationResults[it->first]]->miller(0)->setRejected(
                                                                                  RejectReasonCorrelation, true);
            count++;
        }
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
