#include <string>
#include <iostream>
#include <cmath>
#include "misc.h"
#include "StatisticsManager.h"
#include <new>
#include <fstream>
#include "FileParser.h"
#include "CSV.h"
#include "Reflection.h"
#include "Miller.h"

void StatisticsManager::setMtzs(vector<MtzPtr> newMtzs)
{
	mtzs = newMtzs;
    mtz_num = (int)newMtzs.size();
}


void StatisticsManager::loadFiles(char **filenames, int filenum, int partiality)
{
	std::cout << "Loading " << filenum << " files..." << std::endl;

	for (int i = 0; i < filenum; i++)
	{
		std::cout << "Loading file " << filenames[i] << std::endl;
		std::string filename_string = filenames[i];

		MtzPtr mtz = MtzPtr(new MtzManager());
		mtzs.push_back(mtz);

		mtzs[i]->setFilename(filename_string);
		mtzs[i]->loadReflections(partiality);
	}

	mtz_num = filenum;

	std::cout << "Loaded " << filenum << " files." << std::endl;
}

double StatisticsManager::midPointBetweenResolutions(double minD, double maxD)
{
	vector<double> bins;
	generateResolutionBins(minD, maxD, 2, &bins);

	return bins[1];
}

void StatisticsManager::generateResolutionBins(double minD, double maxD,
		int binCount, vector<double> *bins)
{
	double minRadius = (minD == 0) ? 0 : 1 / minD;
	double maxRadius = 1 / maxD;

	if (maxD <= 0)
	{
		std::cout << "Warning: maximum resolution set to 0. Ignoring.";
        return;
	}

	double maxVolume = pow(maxRadius, 3);
	double minVolume = pow(minRadius, 3);
	double totalVolume = maxVolume - minVolume;

	double eachVolume = totalVolume / binCount;

	double r1 = minRadius;
	double r2 = 0;

	bins->push_back(1 / r1);

	for (int i = 0; i < binCount; i++)
	{
		double r2_cubed = pow(r1, 3) + eachVolume;
		r2 = cbrt(r2_cubed);

		bins->push_back(1 / r2);

		r1 = r2;
	}
}

double StatisticsManager::cc_pearson(int num1, int num2, int silent,
		int *hits, double *multiplicity, double lowResolution,
		double highResolution, bool shouldLog, bool freeOnly)
{
	MtzManager *shot1 = &*(mtzs[num1]);
	MtzManager *shot2 = &*(mtzs[num2]);

	return cc_pearson(shot1, shot2, silent, hits, multiplicity,
			lowResolution, highResolution, shouldLog, freeOnly);
}

void StatisticsManager::convertResolutions(double lowAngstroms,
		double highAngstroms, double *lowReciprocal, double *highReciprocal)
{
	*lowReciprocal = (lowAngstroms == 0) ? 0 : 1 / lowAngstroms;
	*highReciprocal = (highAngstroms == 0) ? FLT_MAX : 1 / highAngstroms;
}

double StatisticsManager::cc_pearson(MtzManager *shot1, MtzManager *shot2,
		int silent, int *hits, double *multiplicity,
		double lowResolution, double highResolution, bool shouldLog, bool freeOnly)
{
	vector<ReflectionPtr> reflections1;
	vector<ReflectionPtr> reflections2;
	int num = 0;
    
	shot1->findCommonReflections(shot2, reflections1, reflections2, &num, true);

    if (reflections1.size() <= 2)
	{
		return -1;
	}

	double invHigh = 1 / highResolution;
	double invLow = 1 / lowResolution;

	if (highResolution == 0)
		invHigh = FLT_MAX;

	if (lowResolution == 0)
		invLow = 0;

	if (!silent || shouldLog)
	{
        std::ostringstream logged;
        
        if (!silent)
            logged << "h k l\tFirst intensity\tSecond intensity\tResolution" << std::endl;
        
        CSVPtr csv = CSVPtr(new CSV(6, "h", "k", "l", "First intensity", "Second intensity", "Resolution"));
        
		for (int i = 0; i < num; i++)
		{
            int h = reflections1[i]->miller(0)->getH();
            int k = reflections1[i]->miller(0)->getK();
            int l = reflections1[i]->miller(0)->getL();
            
            if (freeOnly && !reflections2[i]->miller(0)->isFree())
                continue;

			double int1 = reflections1[i]->meanIntensity();
			double int2 = reflections2[i]->meanIntensity();

			if (int1 != int1 || int2 != int2)
				continue;

			if (!(reflections1[i]->betweenResolutions(lowResolution, highResolution)))
				continue;
            
            double resolution = 1 / reflections1[i]->getResolution();

            csv->addEntry(6, (double)h, (double)k, (double)l, int1, int2, resolution);
			
            if (!silent)
                logged << h << " " << k << " " << l << "\t" << int1 << "\t" << int2 << "\t" << resolution << std::endl;
		}

        if (!silent)
        {
            std::map<std::string, std::string> plotMap;
            plotMap["filename"] = "correlation";
            plotMap["xHeader0"] = "First intensity";
            plotMap["yHeader0"] = "Second intensity";
            plotMap["style0"] = "scatter";
            csv->plotPNG(plotMap);
            csv->writeToFile("correlation.csv");
        }
        csv->plotColumns(3, 4);
        
        logged << "Data logged to correlation.csv" << std::endl;
        Logger::mainLogger->addStream(&logged);
        
	}

	double sum_x = 0;
	double sum_y = 0;
	double weight_counted = 0;
	int num_counted = 0;

	bool useSigmaOne = true;
	bool useSigmaTwo = false;

	for (int i = 0; i < num; i++)
	{
		if (reflections1[i]->miller(0)->getPartiality() < 1)
		{
			useSigmaOne = false;
		}
		if (reflections2[i]->miller(0)->getPartiality() < 1)
		{
			useSigmaTwo = false;
			break;
		}
	}

	for (int i = 0; i < num; i++)
	{
        if (freeOnly && !reflections2[i]->miller(0)->isFree())
            continue;

		if (!(reflections1[i]->getResolution() > invLow
				&& reflections1[i]->getResolution() < invHigh))
			continue;
        
        double weight =
        useSigmaOne ?
        1 / reflections1[i]->meanSigma() :
        reflections1[i]->meanPartiality();
        
        weight *=
        useSigmaTwo ?
        1 / reflections2[i]->meanSigma() :
        reflections2[i]->meanPartiality();

        if (weight < 0)
            continue;
        
		double mean1 = reflections1[i]->meanIntensity();
        double mean2 = reflections2[i]->meanIntensity();

		if (mean1 != mean1 || mean2 != mean2 || weight != weight)
			continue;

		num_counted++;
		sum_x += mean1 * weight;
		sum_y += mean2 * weight;

		weight_counted += weight;
	}

	double mean_x = sum_x / weight_counted;
	double mean_y = sum_y / weight_counted;

	double sum_x_y_minus_mean_x_y = 0;
	double sum_x_minus_mean_x_sq = 0;
	double sum_y_minus_mean_y_sq = 0;

	for (int i = 0; i < num; i++)
	{
        if (freeOnly && !reflections2[i]->miller(0)->isFree())
            continue;

		if (!(reflections1[i]->getResolution() > invLow
				&& reflections1[i]->getResolution() < invHigh))
			continue;

		double amp_x = reflections1[i]->meanIntensity();
		double amp_y = reflections2[i]->meanIntensity();

		double weight =
				useSigmaOne ?
						1 / reflections1[i]->meanSigma() :
						reflections1[i]->meanPartiality();
		weight *=
				useSigmaTwo ?
						1 / reflections2[i]->meanSigma() :
						reflections2[i]->meanPartiality();

        if (weight < 0)
            continue;

		if (amp_x != amp_x || amp_y != amp_y)
			continue;

		if (weight != weight)
			continue;

		sum_x_y_minus_mean_x_y += weight * (amp_x - mean_x) * (amp_y - mean_y);
		sum_x_minus_mean_x_sq += weight * pow(amp_x - mean_x, 2);
		sum_y_minus_mean_y_sq += weight * pow(amp_y - mean_y, 2);
	}

	sum_x_y_minus_mean_x_y /= weight_counted;
	sum_x_minus_mean_x_sq /= weight_counted;
	sum_y_minus_mean_y_sq /= weight_counted;

    if (hits != NULL)
		*hits = num_counted;

	double r = sum_x_y_minus_mean_x_y
			/ (sqrt(sum_x_minus_mean_x_sq * sum_y_minus_mean_y_sq));

	if (r < 0)
		r = 0;
	if (r != r)
		r = -1;

	return r;
}

double StatisticsManager::r_factor(RFactorType rFactor, MtzManager *shot1, int *hits,
		double *multiplicity, double lowResolution, double highResolution, bool freeOnly)
{
	double rmerge_numerator = 0;
	double rmerge_denominator = 0;

	int all_reflections_used = 0;
	int all_millers_used = 0;
    double sumMultiplicity = 0;
    double multiplicityDivide = 0;
	double mult = 0;

	double dMin = 0;
	double dMax = 0;

	convertResolutions(lowResolution, highResolution, &dMin, &dMax);

	for (int i = 0; i < shot1->reflectionCount(); i++)
	{
		ReflectionPtr reflection = shot1->reflection(i);

        if (!reflection->anyAccepted())
            continue;
        
		MillerPtr firstMiller = reflection->miller(0);

		double res = reflection->getResolution();

		if (res < dMin || res > dMax)
			continue;

		double sum_numerator = 0;
		double sum_denominator = 0;
		double n = reflection->millerCount();

        sumMultiplicity += reflection->acceptedCount();
        multiplicityDivide++;
        
		if (reflection->acceptedCount() == 1)
        {
            continue;
        }
        
        all_reflections_used++;

		for (int j = 0; j < n; j++)
		{
			if (!reflection->miller(j)->accepted())
				continue;

			double mean_i = reflection->meanIntensity();
			double weight = reflection->meanPartiality() / reflection->meanSigma();

            weight = 1;
            
			if (weight != weight || mean_i != mean_i)
				continue;

			all_millers_used++;

			double one_i = reflection->miller(j)->intensity();
            
			sum_numerator += fabs(mean_i - one_i) * weight;
			sum_denominator += fabs(one_i) * weight;
		}

		double one_multiple = (double) all_millers_used
				/ (double) all_reflections_used;
		mult += one_multiple;

		double scalingTerm = 1;

		if (rFactor == RFactorTypePim)
		{
			scalingTerm = sqrt(1 / (n - 1));
		}
		else if (rFactor == RFactorTypeMeas)
		{
			scalingTerm = sqrt(n / (n - 1));
		}

		rmerge_numerator += scalingTerm * sum_numerator;
		rmerge_denominator += sum_denominator;
	}

	mult /= all_reflections_used;

	if (multiplicity != NULL)
    {
		*multiplicity = mult;
    }

	if (hits != NULL)
		*hits = all_reflections_used;

	double r_merge = rmerge_numerator / rmerge_denominator;

	return r_merge;
}

double StatisticsManager::r_split(MtzManager *shot1, MtzManager *shot2,
		int silent, int *hits, double *multiplicity,
		double lowResolution, double highResolution, bool shouldLog, bool freeOnly)
{
	double sum_numerator = 0;
	double sum_denominator = 0;
	double mult = 0;
	int count = 0;

	double dMin = 0;
	double dMax = 0;

	convertResolutions(lowResolution, highResolution, &dMin, &dMax);
    std::string filename = shot1->getFilename();

	for (int i = 0; i < shot1->reflectionCount(); i++)
	{
		ReflectionPtr reflection = shot1->reflection(i);
        int reflid = (int)reflection->getReflId();

		ReflectionPtr reflection2;
		shot2->findReflectionWithId(reflid, &reflection2);

		if (!reflection2)
			continue;
        
        if (freeOnly && !reflection2->miller(0)->isFree())
            continue;
        
		double int1 = reflection->meanIntensity();
		double int2 = reflection2->meanIntensity();
        
        double weight = reflection->meanWeight();

        if (int1 == 0 || weight == 0 || weight != weight)
            continue;
        
        if (int1 != int1 || int2 != int2)
            continue;
        
        if (int1 + int2 < 0)
            continue;

		double res = reflection->getResolution();

		if (res < dMin || res > dMax)
			continue;

		count++;
		sum_numerator += fabs(int1 - int2) * weight;
		sum_denominator += fabs((int1 + int2) / 2) * weight;

		mult += reflection->acceptedCount() + reflection2->acceptedCount();
	}

	mult /= count;

	if (hits != NULL)
		*hits = count;

	if (multiplicity != NULL)
		*multiplicity = mult;

	double r_split = sum_numerator / (sum_denominator * sqrt(2));

	return r_split;
}

double StatisticsManager::cc_through_origin(int num1, int num2, int silent,
		int inverted, int *hits)
{
	return cc_through_origin(num1, num2, silent, inverted, hits, 0, 0, 0);
}

double StatisticsManager::cc_through_origin(int num1, int num2, int silent,
		int inverted, int *hits, double lowResolution, double highResolution,
		bool shouldLog)
{
	MtzManager *shot1 = &*(mtzs[num1]);
	MtzManager *shot2 = &*(mtzs[num2]);

	vector<ReflectionPtr> reflections1;
	vector<ReflectionPtr> reflections2;
	int num = 0;

	double highRes = 1 / highResolution;
	double lowRes = 1 / lowResolution;

	if (highResolution == 0)
		highRes = FLT_MAX;

	if (lowResolution == 0)
		lowRes = 0;

	shot1->findCommonReflections(shot2, reflections1, reflections2, &num);
	(*hits) = num;

	if (num <= 1)
		return -1;

	if (!silent)
	{
		for (int i = 0; i < num; i++)
		{
			if (reflections1[i]->getResolution() < lowRes
					|| reflections1[i]->getResolution() > highRes)
			{
				std::cout << "Out of resolution range" << std::endl;
				continue;
			}

			double int1 = reflections1[i]->meanIntensity();
			double int2 = reflections2[i]->meanIntensity();

			if (int1 != int1 || int2 != int2)
				continue;

			std::cout << int1 << "\t" << int2 << std::endl;
		}
	}

	double x_squared = 0;
	double x_y = 0;
	double sum_y = 0;

	double residuals_squared = 0;
	double denominator = 0;

	for (int i = 0; i < num; i++)
	{
		if (reflections1[i]->getResolution() < lowRes
				|| reflections1[i]->getResolution() > highRes)
		{
			continue;
		}

		double mean1 = reflections1[i]->meanIntensity();
		double mean2 = reflections2[i]->meanIntensity();

		if (mean1 != mean1 || mean2 != mean2)
			continue;

		x_squared += pow(mean1, 2);
		x_y += mean1 * mean2;
		sum_y += mean2;
	}

	double mean_y = sum_y / num;
	double grad = x_y / x_squared;

	for (int i = 0; i < num; i++)
	{
		if (reflections1[i]->getResolution() < lowRes
				|| reflections1[i]->getResolution() > highRes)
		{
			continue;
		}

		double mean1 = reflections1[i]->meanIntensity();
		double mean2 = reflections2[i]->meanIntensity();

		if (mean1 != mean1 || mean2 != mean2)
			continue;

		residuals_squared += pow(mean2 - grad * mean1, 2);
		denominator += pow(mean2 - mean_y, 2);
	}

	double R_squared = 1 - residuals_squared / denominator;
    
    if (R_squared < 0)
		R_squared = 0;
	double R = sqrt(R_squared);
    
	if (grad < 0)
		R = 0;
	if (R != R)
		R = -1;

	return R;
}

StatisticsManager::StatisticsManager(void)
{
	mtz_num = 0;
	mtzs.clear();
}

StatisticsManager::~StatisticsManager(void)
{
	mtzs.clear();
	std::vector<MtzPtr>().swap(mtzs);
}

