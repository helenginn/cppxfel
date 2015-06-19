/*
 * MtzManagerRefine.cpp
 *
 *  Created on: 28 Aug 2014
 *      Author: helenginn
 */

#include <liblbfgs.h>
#include "StatisticsManager.h"
#include <cmath>
#include "Vector.h"
#include <algorithm>
#include "MtzManager.h"
#include "FileParser.h"

MtzManager *MtzManager::currentManager;

void MtzManager::applyUnrefinedPartiality()
{
    double wavelength = bestWavelength();
    double mosaicity = getMosaicity();
    double rlpSize = getSpotSize();
    double bandwidth = getBandwidth();
    double exponent = getExponent();

    refreshPartialities(0, 0, mosaicity, rlpSize, wavelength, bandwidth, exponent);
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
}

void MtzManager::refreshPartialities(double parameters[])
{
    refreshPartialities(parameters[0], parameters[1], parameters[2],
			parameters[3], parameters[4], parameters[5], parameters[6]);
}

void MtzManager::refreshPartialities(double hRot, double kRot, double mosaicity,
		double spotSize, double wavelength, double bandwidth, double exponent)
{

    this->makeSuperGaussianLookupTable(exponent);

	for (int i = 0; i < holders.size(); i++)
	{
		for (int j = 0; j < holders[i]->millerCount(); j++)
		{
			MillerPtr miller = holders[i]->miller(j);
			miller->recalculatePartiality(hRot, kRot, mosaicity, spotSize,
					wavelength, bandwidth, exponent);
		}
	}

}

double MtzManager::bestWavelength(double lowRes, double highRes, bool usingReference)
{
	double totalWavelength = 0;
	double count = 0;

    double minWavelength = 1.23;
    double maxWavelength = 1.32;
    
	double low, high;
	StatisticsManager::convertResolutions(lowRes, highRes, &low, &high);

    double refinementIntensityThreshold = FileParser::getKey("REFINEMENT_INTENSITY_THRESHOLD", 200.0);
    
    if (usingReference)
        refinementIntensityThreshold = 0.3;
    
	for (int i = 0; i < holderCount(); i++)
	{
		for (int j = 0; j < holder(i)->millerCount(); j++)
		{
			if (holder(i)->getResolution() < low)
				continue;

			if (holder(i)->getResolution() > high)
				continue;
            
            double wavelength = holder(i)->miller(0)->getWavelength();
            if (wavelength < minWavelength || wavelength > maxWavelength)
                continue;

			double weight = 1;
			double isigi = holder(i)->miller(j)->getRawestIntensity();

            if (usingReference)
            {
                MtzManager *reference = MtzManager::getReferenceManager();
                Holder *refHolder = NULL;
                reference->findHolderWithId(holder(i)->getReflId(), &refHolder);
                
                if (refHolder == NULL)
                    continue;
                
                double imgIntensity = holder(i)->miller(j)->getRawestIntensity();
                double refIntensity = refHolder->meanIntensity();
                
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
    
	return totalWavelength / count;
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
		bool shouldLog)
{
	StatisticsFunction *function = StatisticsManager::r_split;

	return statisticsWithManager(otherManager, function, printHits, silent,
			lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::correlationWithManager(MtzManager *otherManager,
		bool printHits, bool silent, double lowRes, double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	StatisticsFunction *function = StatisticsManager::cc_pearson;

	return statisticsWithManager(otherManager, function, printHits, silent,
			lowRes, highRes, bins, correlations, shouldLog);
}

double MtzManager::statisticsWithManager(MtzManager *otherManager,
		StatisticsFunction *function, bool printHits, bool silent, double lowRes,
		double highRes, int bins,
		vector<boost::tuple<double, double, double, int> > *correlations,
		bool shouldLog)
{
	return statisticsWithManager(otherManager, function, NULL, RFactorNone,
			printHits, silent, lowRes, highRes, bins, correlations, shouldLog);
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
		bool shouldLog)
{
	vector<double> shells;

	if (bins == 0)
		bins = 10;

	double maxResolution = 0;

	for (int i = 0; i < holderCount(); i++)
	{
		if (holder(i)->getResolution() > maxResolution
				&& holder(i)->getResolution() < 1 / 1.4)
		{
            Holder *bestHolder = holder(i);
            
			maxResolution = bestHolder->getResolution();
		}
	}

	if (highRes == 0 || highRes < 1 / maxResolution)
		highRes = 1 / maxResolution;

	int hits = 0;
	double multiplicity = 0;

	double statistic = 0;

	if (rFactor == RFactorNone)
		statistic = function(this, otherManager, !printHits, &hits,
				&multiplicity, lowRes, highRes, shouldLog);
	else
		statistic = rFactorFunction(rFactor, this, &hits, &multiplicity, lowRes,
				highRes);

	if (bins > 1 || !silent)
	{
		StatisticsManager::generateResolutionBins(lowRes, highRes, bins,
				&shells);

		for (int i = 0; i < shells.size() - 1; i++)
		{
			double statistic = 0;

			if (rFactor == RFactorNone)
				statistic = function(this, otherManager, 1, &hits,
						&multiplicity, shells[i], shells[i + 1], shouldLog);
			else
				statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
						shells[i], shells[i + 1]);

			if (!silent)
			{
				std::cout << "N: " << shells[i] << "\t" << shells[i + 1] << "\t"
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
					&multiplicity, lowRes, highRes, shouldLog);
		else
			statistic = rFactorFunction(rFactor, this, &hits, &multiplicity,
					lowRes, highRes);

		std::cout << "N: *** Overall ***" << std::endl;
		std::cout << "N: " << lowRes << "\t" << highRes << "\t" << statistic
				<< "\t" << hits << "\t" << multiplicity << std::endl;
	}

	return statistic;
}

double MtzManager::wavelengthStandardDeviation()
{
    double minResolution = maxResolutionRlpSize;
    double maxResolution = maxResolutionAll;
    int reflectionCount = 100;
    
    vector<Holder *>refHolders, imageHolders;
    std::map<double, double> wavelengths;
    
    this->findCommonReflections(getReferenceManager(), imageHolders, refHolders);
    
    for (int i = 0; i < imageHolders.size(); i++)
    {
        double refIntensity = refHolders[i]->meanIntensity();
        if (refIntensity < 1000)
            continue;
        
        if (!imageHolders[i]->betweenResolutions(minResolution, maxResolution))
            continue;
        
        for (int j = 0; j < imageHolders[i]->millerCount(); j++)
        {
            double intensity = imageHolders[i]->miller(j)->getRawestIntensity();
            
            double wavelength = imageHolders[i]->miller(0)->getWavelength();
            
            if (wavelength < 1.25 || wavelength > 1.33)
                continue;
            
            double proportion = intensity / refIntensity;
            
            wavelengths[proportion] = wavelength;
        }
    }
    
    if (wavelengths.size() == 0)
        return FLT_MAX;
    
    int count = 0;
    vector<double> bigWavelengths;
    
    std::map<double, double>::iterator it = wavelengths.end();
    it--;
    
    for (; it != wavelengths.begin() && count < reflectionCount; --it)
    {
        bigWavelengths.push_back(it->second);
        count++;
    }
    
    
    double stDev = standard_deviation(&bigWavelengths);
    this->wavelength = median(&bigWavelengths);
    
//    std::cout << "St dev: " << stDev << ", median: " << wavelength << std::endl;
    
    return stDev;
}