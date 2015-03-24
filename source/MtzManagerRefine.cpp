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

double MtzManager::bestWavelength(double lowRes, double highRes)
{
	double totalWavelength = 0;
	double count = 0;

	double low, high;
	StatisticsManager::convertResolutions(lowRes, highRes, &low, &high);

    double refinementIntensityThreshold = FileParser::getKey("REFINEMENT_INTENSITY_THRESHOLD", 200.0);
    
	for (int i = 0; i < holderCount(); i++)
	{
		for (int j = 0; j < holder(i)->millerCount(); j++)
		{
			if (holder(i)->getResolution() < low)
				continue;

			if (holder(i)->getResolution() > high)
				continue;

			double weight = 1;
			double isigi = holder(i)->miller(j)->getRawIntensity();

			if (isigi > refinementIntensityThreshold && isigi == isigi)
			{
				totalWavelength += holder(i)->miller(j)->getWavelength()
						* weight;
				count += weight;
			}
		}
	}

	return totalWavelength / count;
}

double MtzManager::correlation(bool silent, double lowResolution,
		double highResolution)
{
	if (highResolution == -1)
		highResolution = maxResolutionAll;

	double correlation = StatisticsManager::cc_pearson(this,
			MtzManager::referenceManager, true, isInverse(), NULL,
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
			maxResolution = holder(i)->getResolution();
		}
	}

	if (highRes == 0 || highRes < 1 / maxResolution)
		highRes = 1 / maxResolution;

	int hits = 0;
	bool inverse = this->isInverse();
	double multiplicity = 0;

	double statistic = 0;

	if (rFactor == RFactorNone)
		statistic = function(this, otherManager, !printHits, inverse, &hits,
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
				statistic = function(this, otherManager, 1, inverse, &hits,
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
			statistic = function(this, otherManager, 1, inverse, &hits,
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
