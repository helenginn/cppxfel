/*
 * MtzManagerRefine.cpp
 *
 *  Created on: 28 Aug 2014
 *      Author: helenginn
 */

#include "StatisticsManager.h"
#include <cmath>
#include "Vector.h"
#include <algorithm>
#include "MtzManager.h"
#include "FileParser.h"
#include "GaussianBeam.h"
#include "SpectrumBeam.h"
#include "RefinementStepSearch.h"
#include "Reflection.h"
#include "Miller.h"

void MtzManager::applyUnrefinedPartiality()
{
    if (wavelength == 0)
    {
        wavelength = bestWavelength();
    }
    
    double mosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);
    double rlpSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
    double bandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
    double exponent = FileParser::getKey("INITIAL_EXPONENT", INITIAL_EXPONENT);

    refreshPartialities(0, 0, mosaicity, rlpSize, wavelength, bandwidth, exponent);
}

void MtzManager::addParameters(RefinementStrategyPtr refiner)
{
	double stepSizeUnitCellA = FileParser::getKey("STEP_SIZE_UNIT_CELL_A", 0.5);
	double stepSizeUnitCellB = FileParser::getKey("STEP_SIZE_UNIT_CELL_B", 0.5);
	double stepSizeUnitCellC = FileParser::getKey("STEP_SIZE_UNIT_CELL_C", 0.5);

	bool optimisingUnitCellA = FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
	bool optimisingUnitCellB = FileParser::getKey("OPTIMISING_UNIT_CELL_B", false);
	bool optimisingUnitCellC = FileParser::getKey("OPTIMISING_UNIT_CELL_C", false);


    if (optimisingOrientation)
    {
        refiner->addParameter(this, getHRot, setHRot, stepSizeOrientation, toleranceOrientation, "hRot");
        refiner->addParameter(this, getKRot, setKRot, stepSizeOrientation, toleranceOrientation, "kRot");
    }
    
    if (optimisingRlpSize)
	{
        refiner->addParameter(this, getSpotSizeStatic, setSpotSizeStatic, stepSizeRlpSize, toleranceRlpSize, "rlpSize");
	}

    if (optimisingMosaicity)
	{
        refiner->addParameter(this, getMosaicityStatic, setMosaicityStatic, stepSizeMosaicity, toleranceMosaicity, "mosaicity");
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


    // delete me later

	if (optimisingWavelength)
	{
		refiner->addParameter(this, getWavelengthStatic, setWavelengthStatic, stepSizeWavelength, toleranceWavelength, "wavelength");
	}

	if (optimisingBandwidth)
	{
		refiner->addParameter(this, getBandwidthStatic, setBandwidthStatic, stepSizeBandwidth, 0, "bandwidth");
	}

	if (optimisingExponent)
	{
		refiner->addParameter(this, getExponentStatic, setExponentStatic, stepSizeExponent, 0, "exponent");
	}
}

double MtzManager::refinePartialitiesOrientation(int ambiguity, bool reset)
{
    this->setActiveAmbiguity(ambiguity);
	scoreType = defaultScoreType;

	refineParameterScore(this);
    
	RefinementStrategyPtr refinementMap = RefinementStrategy::userChosenStrategy();

	if (wavelength == 0)
	{
		wavelength = bestWavelength();
	}

	/*
    {
        beam = GaussianBeamPtr(new GaussianBeam(wavelength, bandwidth, exponent));
        
        for (int i = 0; i < reflectionCount(); i++)
        {
            for (int j = 0; j < reflection(i)->millerCount(); j++)
            {
                reflection(i)->miller(j)->setBeam(beam);
            }
        }
    }

    beam->addParameters(refinementMap);
*/

	addParameters(refinementMap);

	refinementMap->setJobName("Refining " + getFilename());
    refinementMap->setEvaluationFunction(refineParameterScore, this);
    
    refinementMap->refine();
    
    double correl = correlation();

	if (reset)
	{
		refinementMap->resetToInitialParameters();
	}

	return correl;
}

void MtzManager::refinePartialities()
{
    std::vector<double> correlations;
    double maxCorrel = -1;
    int bestAmbiguity = 0;
    
    for (int i = 0; i < ambiguityCount(); i++)
    {
        correlations.push_back(refinePartialitiesOrientation(i));
    }
    
    for (int i = 0; i < ambiguityCount(); i++)
    {
        if (correlations[i] > maxCorrel)
        {
            bestAmbiguity = i;
            maxCorrel = correlations[i];
        }
    }
    
    refinePartialitiesOrientation(bestAmbiguity, false);
    double partCorrel = leastSquaresPartiality();
    setRefPartCorrel(partCorrel);
    
    setRefCorrelation(correlation());
    
    double rSplitValue = rSplit(0, 0);
    
    logged << getFilename() << "\t" << describeScoreType() << "\t\t"
	<< refCorrelation << "\t" << rSplitValue << "\t" << refPartCorrel
	<< "\t" << accepted() << std::endl;
    sendLog();
    
    writeToFile(std::string("ref-") + getFilename());
}

void MtzManager::refreshCurrentPartialities()
{
    if (externalScale != -1)
        this->applyScaleFactor(externalScale, 0, 0, true);
    
    applyBFactor(bFactor);
    refreshPartialities(this->hRot,
                        this->kRot,
                        this->mosaicity,
                        this->spotSize,
                        this->wavelength,
                        this->bandwidth,
                        this->exponent);
}

void MtzManager::replaceBeamWithSpectrum()
{
    GaussianBeamPtr gaussian = boost::static_pointer_cast<GaussianBeam>(beam);
    
    if (gaussian)
    {
        SpectrumBeamPtr spectrum = SpectrumBeamPtr(new SpectrumBeam(gaussian, shared_from_this()));
        
        beam = spectrum;
        
        for (int i = 0; i < reflectionCount(); i++)
        {
            for (int j = 0; j < reflection(i)->millerCount(); j++)
            {
                reflection(i)->miller(j)->setBeam(spectrum);
            }
        }
    }
}

void MtzManager::refreshPartialities(double hRot, double kRot, double mosaicity,
		double spotSize, double wavelength, double bandwidth, double exponent)
{
    this->makeSuperGaussianLookupTable(exponent);

    if (!matrix)
        return;
    
    
	lockUnitCellDimensions();

    this->matrix->changeOrientationMatrixDimensions(getUnitCell());

    MatrixPtr newMatrix = MatrixPtr();
    Miller::rotateMatrixHKL(hRot, kRot, 0, matrix, &newMatrix);
    
	for (int i = 0; i < reflections.size(); i++)
	{
		for (int j = 0; j < reflections[i]->millerCount(); j++)
		{
			MillerPtr miller = reflections[i]->miller(j);
			miller->recalculatePartiality(newMatrix, mosaicity, spotSize,
					wavelength, bandwidth, exponent);
		}
	}
    
    logged << "Refreshed partialities with wavelength " << wavelength << ", spot size " << spotSize << " to generate " << accepted() << " accepted reflections." << std::endl;

    sendLog(LogLevelDebug);
}

static bool greaterThan(double num1, double num2)
{
    return (num1 > num2);
}

double MtzManager::medianWavelength(double lowRes, double highRes)
{
    vector<double> wavelengths;
    double refinementIntensityThreshold = FileParser::getKey("REFINEMENT_INTENSITY_THRESHOLD", 200.0);
    
    vector<double> wavelengthRange = FileParser::getKey("WAVELENGTH_RANGE", vector<double>());
    
    for (int i = 0; i < reflectionCount(); i++)
    {
        for (int j = 0; j < reflection(i)->millerCount(); j++)
        {
            double wavelength = reflection(i)->miller(0)->getWavelength();
            
            if (wavelengthRange.size())
            {
                if (wavelength < wavelengthRange[0] || wavelength > wavelengthRange[1])
                    continue;
            }
            
            double isigi = reflection(i)->miller(j)->getRawestIntensity();
            
            if (isigi > refinementIntensityThreshold && isigi == isigi)
            {
                wavelengths.push_back(wavelength);
            }
        }
    }
    
    std::sort(wavelengths.begin(), wavelengths.end(), greaterThan);

    if (wavelengths.size() < 2)
        return 0;
    
    return wavelengths[int(wavelengths.size() / 2)];
}

double MtzManager::bestWavelength(double lowRes, double highRes, bool usingReference)
{
    bool median = FileParser::getKey("MEDIAN_WAVELENGTH", false);
    
    if (median)
        return medianWavelength(lowRes, highRes);
    
    vector<double> wavelengthRange = FileParser::getKey("WAVELENGTH_RANGE", vector<double>());
    
	double totalWavelength = 0;
	double count = 0;

    double low, high;
	StatisticsManager::convertResolutions(lowRes, highRes, &low, &high);

    double refinementIntensityThreshold = FileParser::getKey("REFINEMENT_INTENSITY_THRESHOLD", 200.0);
    
    if (usingReference)
        refinementIntensityThreshold = 0.3;
    
	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			if (reflection(i)->getResolution() < low)
				continue;

			if (reflection(i)->getResolution() > high)
				continue;
            
            double wavelength = reflection(i)->miller(0)->getWavelength();
            
            if (wavelengthRange.size())
            {
                if (wavelength < wavelengthRange[0] || wavelength > wavelengthRange[1])
                {
                    continue;
                }
            }
            
			double weight = 1;
			double isigi = reflection(i)->miller(j)->getRawestIntensity();

            MtzManager *reference = MtzManager::getReferenceManager();
            ReflectionPtr refReflection = ReflectionPtr();
            
            if (reference != NULL)
            {
                reference->findReflectionWithId(reflection(i)->getReflId(), &refReflection);
                
                if (refReflection && refReflection->meanIntensity() < REFERENCE_WEAK_REFLECTION)
                    continue;
            }
            
            if (usingReference)
            {
                if (!refReflection)
                    continue;
                
                double imgIntensity = reflection(i)->miller(j)->getRawestIntensity();
                double refIntensity = refReflection->meanIntensity();
                
                double proportion = imgIntensity / refIntensity;
                isigi = proportion;
            }
            
			if (isigi > refinementIntensityThreshold && isigi == isigi)
			{
				totalWavelength += wavelength
						* weight;
				count += weight;
			}
		}
	}

    if (count == 0)
    {
        this->setRejected(true);
    }

    double newWavelength = totalWavelength / count;
    
    logged << "Best wavelength " << newWavelength << " Å from " << count << " reflections." << std::endl;
    sendLog(LogLevelDetailed);
    return newWavelength;

}

double MtzManager::correlation(bool silent, double lowResolution,
		double highResolution)
{
	if (highResolution == -1)
		highResolution = maxResolutionAll;

	double correlation = StatisticsManager::cc_pearson(this,
			MtzManager::referenceManager, true, NULL,
			NULL, lowResolution, highResolution, false);

	return correlation;
}

double MtzManager::rSplitWithManager(MtzManager *otherManager, bool printHits,
		bool silent, double lowRes, double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog, bool freeOnly)
{
	StatisticsFunction *function = StatisticsManager::r_split;

	return statisticsWithManager(otherManager, function, printHits, silent,
			lowRes, highRes, bins, correlations, shouldLog, freeOnly);
}

double MtzManager::correlationWithManager(MtzManager *otherManager,
		bool printHits, bool silent, double lowRes, double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog, bool freeOnly)
{
	StatisticsFunction *function = StatisticsManager::cc_pearson;

	return statisticsWithManager(otherManager, function, printHits, silent,
			lowRes, highRes, bins, correlations, shouldLog, freeOnly);
}

double MtzManager::statisticsWithManager(MtzManager *otherManager,
		StatisticsFunction *function, bool printHits, bool silent, double lowRes,
		double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog, bool freeOnly)
{
	return statisticsWithManager(otherManager, function, NULL, RFactorNone,
			printHits, silent, lowRes, highRes, bins, correlations, shouldLog, freeOnly);
}

double MtzManager::rFactorWithManager(RFactorType rFactor, bool printHits,
		bool silent, double lowRes, double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	RFactorFunction *function = StatisticsManager::r_factor;
	return statisticsWithManager(NULL, NULL, function, rFactor, printHits,
			silent, lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::statisticsWithManager(MtzManager *otherManager,
		StatisticsFunction *function, RFactorFunction *rFactorFunction,
		RFactorType rFactor, bool printHits, bool silent, double lowRes,
		double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog, bool freeOnly)
{
	vector<double> shells;

	if (bins == 0)
		bins = 10;

	double maxResolution = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		if (reflection(i)->getResolution() > maxResolution
				&& reflection(i)->getResolution() < 1 / 1.2)
		{
            ReflectionPtr bestReflection = reflection(i);
            
			maxResolution = bestReflection->getResolution();
		}
	}

	if (highRes == 0 || highRes < 1 / maxResolution)
		highRes = 1 / maxResolution;

	int hits = 0;
	double multiplicity = 0;

	double statistic = 0;

	if (rFactor == RFactorNone)
		statistic = function(this, otherManager, !printHits, &hits,
				&multiplicity, lowRes, highRes, false, freeOnly);
	else
		statistic = rFactorFunction(rFactor, this, &hits, &multiplicity, lowRes,
				highRes, freeOnly);

    if (!silent)
        logged << "N: " << "lowRes\thighRes\tValue\tHits\tMultiplicity" << std::endl;
    
	if (bins > 1 || !silent)
	{
		StatisticsManager::generateResolutionBins(lowRes, highRes, bins,
				&shells);

		for (int i = 0; i < shells.size() - 1; i++)
		{
			double statistic = 0;

			if (rFactor == RFactorNone)
				statistic = function(this, otherManager, 1, &hits,
						&multiplicity, shells[i], shells[i + 1], false, freeOnly);
			else
				statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
						shells[i], shells[i + 1], freeOnly);

			if (!silent)
			{
				logged << "N: " << shells[i] << "\t" << shells[i + 1] << "\t"
						<< statistic << "\t" << hits << "\t" << multiplicity
						<< std::endl;
			}

			if (correlations != NULL)
			{
				boost::tuple<double, double, double, int> result = boost::make_tuple(
						shells[i], shells[i + 1], statistic, hits);
				correlations->push_back(result);
			}
		}
	}

	if (!silent)
	{
		double statistic = 0;

		if (rFactor == RFactorNone)
			statistic = function(this, otherManager, 1,  &hits,
					&multiplicity, lowRes, highRes, shouldLog, freeOnly);
		else
			statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
					lowRes, highRes, freeOnly);

		logged << "N: *** Overall ***" << std::endl;
		logged << "N: " << lowRes << "\t" << highRes << "\t" << statistic
				<< "\t" << hits << "\t" << multiplicity << std::endl;
	}

    sendLog();
    
	return statistic;
}


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

	hRot = 0;
	kRot = 0;
	allowTrust = FileParser::getKey("ALLOW_TRUST", true);
	bool alwaysTrust = FileParser::getKey("TRUST_INDEXING_SOLUTION", false);

	if (alwaysTrust)
		trust = TrustLevelGood;
}

