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

    refreshPartialities(0, 0, mosaicity, rlpSize, wavelength, bandwidth, exponent, cellDim[0], cellDim[1], cellDim[2]);
}

void MtzManager::setParams(double parameters[], int paramCount)
{
    wavelength = parameters[PARAM_WAVELENGTH];
	bandwidth = parameters[PARAM_BANDWIDTH];
	mosaicity = parameters[PARAM_MOS];
	spotSize = parameters[PARAM_SPOT_SIZE];
	hRot = parameters[PARAM_HROT];
	kRot = parameters[PARAM_KROT];
    exponent = parameters[PARAM_EXPONENT];
    bFactor = parameters[PARAM_B_FACTOR];
    
    applyBFactor(bFactor);
    applyScaleFactor(parameters[PARAM_SCALE_FACTOR], 0, 0, true);
    
    cellDim[0] = parameters[PARAM_UNIT_CELL_A];
    cellDim[1] = parameters[PARAM_UNIT_CELL_B];
    cellDim[2] = parameters[PARAM_UNIT_CELL_C];
}

void MtzManager::getParams(double *parameters[], int paramCount)
{
	(*parameters)[PARAM_WAVELENGTH] = wavelength;
	(*parameters)[PARAM_BANDWIDTH] = bandwidth;
	(*parameters)[PARAM_MOS] = mosaicity;
	(*parameters)[PARAM_SPOT_SIZE] = spotSize;
	(*parameters)[PARAM_HROT] = hRot;
	(*parameters)[PARAM_KROT] = kRot;
    (*parameters)[PARAM_EXPONENT] = exponent;
    (*parameters)[PARAM_B_FACTOR] = bFactor;
    (*parameters)[PARAM_SCALE_FACTOR] = scale;
    (*parameters)[PARAM_UNIT_CELL_A] = cellDim[0];
    (*parameters)[PARAM_UNIT_CELL_B] = cellDim[1];
    (*parameters)[PARAM_UNIT_CELL_C] = cellDim[2];
}

void MtzManager::addParameters(RefinementStrategyPtr map)
{
    if (optimisingOrientation)
    {
        map->addParameter(this, getHRotStatic, setHRotStatic, stepSizeOrientation, toleranceOrientation);
        map->addParameter(this, getKRotStatic, setKRotStatic, stepSizeOrientation, toleranceOrientation);
    }
    
    if (optimisingRlpSize)
        map->addParameter(this, getSpotSizeStatic, setSpotSizeStatic, stepSizeRlpSize, toleranceRlpSize);
    
    if (optimisingMosaicity)
        map->addParameter(this, getMosaicityStatic, setMosaicityStatic, stepSizeMosaicity, toleranceMosaicity);
    
    // delete me later
 
    map->addParameter(this, getWavelengthStatic, setWavelengthStatic, stepSizeWavelength, toleranceWavelength);
    
}

double MtzManager::refinePartialitiesOrientation(int ambiguity, int cycles)
{
  //  logged << "New refinement, ambiguity " << ambiguity << std::endl;
 //  sendLog();
    
    this->setActiveAmbiguity(ambiguity);
    
	RefinementStrategyPtr refinementMap = RefinementStrategy::userChosenStrategy();

    wavelength = bestWavelength();
    
    if (!beam)
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
    addParameters(refinementMap);

    refinementMap->setEvaluationFunction(refineParameterScore, this);
    refinementMap->setCycles(cycles);
    
    refinementMap->refine();
    
    return correlation();
}

void MtzManager::refinePartialities()
{
    std::vector<double> correlations;
    double maxCorrel = -1;
    int bestAmbiguity = 0;
    
    for (int i = 0; i < ambiguityCount(); i++)
    {
        correlations.push_back(refinePartialitiesOrientation(i, 10));
    }
    
    for (int i = 0; i < ambiguityCount(); i++)
    {
        if (correlations[i] > maxCorrel)
        {
            bestAmbiguity = i;
            maxCorrel = correlations[i];
        }
    }
    
    refinePartialitiesOrientation(bestAmbiguity);
    double partCorrel = leastSquaresPartiality();
    setRefPartCorrel(partCorrel);
    
    setRefCorrelation(correlation());
    
    double rSplitValue = rSplit(0, 0);
    
    logged << getFilename() << "\t" << describeScoreType() << "\t\t" << refCorrelation << "\t" << rSplitValue << "\t" << refPartCorrel << "\t" << accepted() << std::endl;
    sendLog();
    
    writeToFile(std::string("ref-") + getFilename());
}

void MtzManager::refreshPartialities(double parameters[])
{
    refreshPartialities(parameters[PARAM_HROT],
                        parameters[PARAM_KROT],
                        parameters[PARAM_MOS],
                        parameters[PARAM_SPOT_SIZE],
                        parameters[PARAM_WAVELENGTH],
                        parameters[PARAM_BANDWIDTH],
                        parameters[PARAM_EXPONENT],
                        parameters[PARAM_UNIT_CELL_A],
                        parameters[PARAM_UNIT_CELL_B],
                        parameters[PARAM_UNIT_CELL_C]);
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
                        this->exponent,
                        this->cellDim[0],
                        this->cellDim[1],
                        this->cellDim[2]);
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
		double spotSize, double wavelength, double bandwidth, double exponent,
                                     double a, double b, double c)
{
    this->makeSuperGaussianLookupTable(exponent);

    if (!matrix)
        return;
    
    
    double spgNum = getLowGroup()->spg_num;
    if (spgNum >= 75 && spgNum <= 194)
    {
        b = a;
    }
    if (spgNum >= 195)
    {
        b = a;
        c = a;
    }

    this->matrix->changeOrientationMatrixDimensions(a, b, c, cellAngles[0], cellAngles[1], cellAngles[2]);

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
