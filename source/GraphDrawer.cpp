/*
 * GraphDrawer.cpp
 *
 *  Created on: 22 Oct 2014
 *      Author: helenginn
 */

#include "GraphDrawer.h"
#include "misc.h"
#include <vector>
#include<iostream>
#include <algorithm>
#include <cmath>
#include <boost/variant.hpp>
#include <map>
#include <sstream>
#include <tuple>
#include "Panel.h"
#include <fstream>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/space_group.h>
#include "Reflection.h"
#include "Miller.h"
#include "UnitCellLattice.h"
#include "FileParser.h"
#include "PNGFile.h"
#include "Image.h"

GraphDrawer::GraphDrawer(MtzManager *mtz)
{
	// TODO Auto-generated constructor stub
	this->mtz = mtz;
}

bool sortByWavelength(Partial x, Partial y)
{
	return (x.wavelength < y.wavelength);
}

bool sortByPercentage(Partial x, Partial y)
{
	return (x.percentage < y.percentage);
}



std::string GraphDrawer::generateFilename(std::string stem)
{
	return generateFilename(stem, "png");
}

std::string GraphDrawer::generateFilename(std::string stem, std::string ext)
{
	time_t cputime;
	time(&cputime);

	std::ostringstream timeString;
	timeString << stem;

	timeString << "." << ext;

	return timeString.str();
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<vector<double> > xs,
		vector<vector<double> > ys)
{
	vector<double> x2, y2;

	return plot(filename, properties, xs, ys, x2, y2);
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<double> x, vector<double> y)
{
	vector<double> x2, y2;

	return plot(filename, properties, x, y, x2, y2);
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<double> x, vector<double> y, vector<double> x2,
		vector<double> y2)
{
	vector<vector<double> > xs;
	vector<vector<double> > ys;

    xs.push_back(x); xs.push_back(x2);
    ys.push_back(y); ys.push_back(y2);

	return plot(filename, properties, xs, ys, x2, y2);
}

std::string GraphDrawer::plot(std::string filename, GraphMap properties,
		vector<vector<double> > xs,
		vector<vector<double> > ys, vector<double> x2,
		vector<double> y2)
{
    std::string csvName = generateFilename(filename, "csv");
    
    std::ofstream csv;
    csv.open(csvName);
    
    int biggest = 0;
    
    for (int i = 0; i < xs.size(); i++)
    {
        if (xs[i].size() > biggest)
            biggest = (int)xs[i].size();
    }
    
    for (int i = 0; i < biggest; i++)
    {
        for (int j = 0; j < xs.size(); j++)
        {
            if (i < xs[j].size())
                csv << xs[j][i] << "," << ys[j][i] << ",";
            else csv << ",,";
            
        }
        
        csv << std::endl;
    }
    
    csv.close();
    
    return csvName;
}

void GraphDrawer::correlationPlot(std::string filename, double xMax, double yMax)
{
	GraphMap properties;

	properties["title"] = "Correlation plot";
	properties["xTitle"] = "Intensity (image)";
	properties["yTitle"] = "Intensity (reference MTZ)";
	properties["plotType"] = "point";

	if (xMax > 0)
		properties["xMax"] = xMax;

	if (yMax > 0)
		properties["yMax"] = yMax;

	vector<double> x;
	vector<double> y;

	MtzManager *shot1 = this->mtz;
	MtzManager *shot2 = mtz->getReferenceManager();
//	double scale = mtz->gradientAgainstManager(*shot2);
//	mtz->applyScaleFactor(scale);

	vector<ReflectionPtr> reflections1;
	vector<ReflectionPtr> reflections2;
	int num = 0;

	shot1->findCommonReflections(shot2, reflections1, reflections2, &num);

	for (int i = 0; i < num; i++)
	{
		double int1 = reflections1[i]->meanIntensity();
		std::string filename = shot1->getFilename();
		double int2 = reflections2[i]->meanIntensity();

		if (int1 != int1 || int2 != int2)
			continue;

		x.push_back(int1);
		y.push_back(int2);
	}

	std::string fullName = this->plot(filename, properties, x, y);
	std::cout << "N: plot " << fullName << std::endl;
}

void GraphDrawer::resolutionStatsPlot(vector<MtzManager *>& managers,
		std::string filename, GraphMap properties, bool intensityBins, bool image)
{
	if (!intensityBins)
		properties["title"] = "Intensity agreement vs resolution";
	else
		properties["title"] = "Intensity agreement vs strength of signal";

	if (!intensityBins)
	{
		properties["xTitle"] = "Resolution (1 / d^2)";
	}
	else
	{
		if (image)
		{
			properties["xTitle"] = "Log(image intensity)";
		}
		else
		{
			properties["xTitle"] = "Reference intensity";
		}
	}

	properties["yTitle"] =
			"Percentage of merged data set (individual reflections)";
	properties["plotType"] = "point";

	properties["xMin"] = 0;
//	properties["xMax"] = 0.5;
	properties["yMin"] = 0;
	properties["yMax"] = 4;

	properties["colour_0"] = 10;
	properties["colour_1"] = 4;

	vector<double> bins;
	StatisticsManager::generateResolutionBins(50, 1.6, 18, &bins);

	MtzManager *referenceManager = MtzManager::getReferenceManager();
	double grad = 1000 / referenceManager->averageIntensity();
	referenceManager->applyScaleFactor(grad);

	vector<double> blueX, blueY, redX, redY;
    
    double percentSum = 0;
    double sum = 0;
    
    for (int i = 0; i < managers.size(); i++)
	{
		mtz = managers[i];

	//	mtz->excludeFromLogCorrelation();
        
        double correl = mtz->correlation(true);
        double invCorrel = correl;
        
        if (correl < 0.85)
            continue;
        
        if (MtzManager::getReferenceManager()->ambiguityCount() == 2)
        {
            mtz->setActiveAmbiguity(1);
            invCorrel = mtz->correlation(true);
            mtz->setActiveAmbiguity(0);
            
            if (invCorrel > correl)
                mtz->setActiveAmbiguity(1);
        }
        double newCorrel = mtz->correlation(true);
        mtz->setRefCorrelation(newCorrel);
        
        std::cout << mtz->getFilename() << "\t" << correl << "\t"
        << invCorrel << std::endl;

		double scale = mtz->gradientAgainstManager(referenceManager);
        std::cout << "Scale: " << scale << std::endl;
        mtz->applyScaleFactor(scale);
     
		vector<ReflectionPtr> refReflections, imgReflections;

		mtz->findCommonReflections(referenceManager, imgReflections, refReflections,
		NULL);

		for (int j = 0; j < refReflections.size(); j++)
		{
			if (!imgReflections[j]->anyAccepted())
				continue;

			double imgIntensity = imgReflections[j]->meanIntensity();
			double refIntensity = refReflections[j]->meanIntensity();
            
			double percent = imgIntensity / refIntensity;

			vector<double> *chosenX = &redX;
			vector<double> *chosenY = &redY;

			if (log(imgIntensity) > 7)
			{
				chosenX = &blueX;
				chosenY = &blueY;
			}
            
			if (!intensityBins)
			{
				double resolution = refReflections[j]->getResolution();
				double res_squared = pow(resolution, 2);
				chosenX->push_back(res_squared);
			}
			else
			{
				double logIntensity = 0;
				if (image)
					logIntensity = log(imgIntensity);
				else
					logIntensity = log(refIntensity);

				chosenX->push_back(logIntensity);
			}
			chosenY->push_back(percent);

        }
    }
    
    double aveSum = percentSum / (double)sum;
    
    std::cout << "Average sum = " << aveSum << std::endl;
    
	if (redX.size() == 0)
		properties["colour_0"] = 1;

	vector<vector<double> > xs, ys;
	xs.push_back(blueX);
	xs.push_back(redX);
	ys.push_back(blueY);
	ys.push_back(redY);

	plot(filename, properties, xs, ys);
}

void GraphDrawer::resolutionStatsCSV(std::vector<MtzManager *>& managers)
{
    double maxRes = this->mtz->maxResolution();
    std::vector<std::vector<double> > intensitiesPerBinList;
    vector<double> bins;
    StatisticsManager::generateResolutionBins(50, 1.7, 6, &bins);

    for (int i = 0; i < managers.size(); i++)
    {
        MtzManager *myMtz = managers[i];
        
/*        if (myMtz->accepted() < 130)
            continue;
  */
        std::vector<ReflectionPtr> imgRefls, refRefls;
        this->mtz->findCommonReflections(myMtz, imgRefls, refRefls);
        
        double correl = myMtz->correlation(true);
        
        if (MtzManager::getReferenceManager()->ambiguityCount() == 2)
        {
            myMtz->setActiveAmbiguity(1);
            double invCorrel = myMtz->correlation(true);
            
            if (invCorrel < correl)
                myMtz->setActiveAmbiguity(0);
        }

        double scale = myMtz->gradientAgainstManager(this->mtz);
        myMtz->applyScaleFactor(scale);
 
      /*
        double scale, bFactor;
        myMtz->bFactorAndScale(&scale, &bFactor);
        myMtz->applyScaleFactor(scale);
        myMtz->applyBFactor(bFactor);*/
        
        std::vector<double> intensitiesPerBin;
        
        for (int j = 0; j < bins.size() - 1; j++)
        {
            double lowResCut = bins[j];
            double highResCut = bins[j + 1];
            
            double imgSum = 0;
            double refSum = 0;
            
            for (int k = 0; k < imgRefls.size(); k++)
            {
                ReflectionPtr imgRefl = imgRefls[k];
                
                if (!imgRefl->betweenResolutions(lowResCut, highResCut))
                {
                    continue;
                }
                
                if (!imgRefl->anyAccepted())
                    continue;
                
                ReflectionPtr refRefl = refRefls[k];
                
                double imgIntensity = imgRefl->meanIntensity();
                double refIntensity = refRefl->meanIntensity();
                
                if (imgIntensity != imgIntensity || refIntensity != refIntensity)
                {
                    continue;
                }
                
                imgSum += imgIntensity;
                refSum += refIntensity;
            }
            
            double ratio = imgSum / refSum;
            
            if (refSum == 0)
            {
                ratio = 0;
            }
            
            intensitiesPerBin.push_back(ratio);
        }
        
        intensitiesPerBinList.push_back(intensitiesPerBin);
    }
  /*
    if (!intensitiesPerBinList.size())
        return;
*/
    for (int i = 0; i < intensitiesPerBinList[0].size(); i++)
    {
        for (int j = 0; j < intensitiesPerBinList.size(); j++)
        {
            std::cout << 1 / pow(bins[i], 2) << ", ";
            std::cout << intensitiesPerBinList[j][i] << ", ";
        }
        
        std::cout << std::endl;
    }
    
    std::cout << "Behaved" << std::endl;
    
    for (int j = 0; j < intensitiesPerBinList.size(); j++)
    {
        if (intensitiesPerBinList[j][0] < 1.15)
        {
            std::cout << "good ";
        }
        else
        {
            std::cout << "bad ";
        }
        
        std::cout << managers[j]->getFilename();
        
        double correl = managers[j]->correlation();
        double scale = managers[j]->getScale();
        double bFactor = managers[j]->bFactor;
        
        std::cout << " " << correl << " " << scale << " " << bFactor << std::endl;
    }
}

void GraphDrawer::bFactorPlot(vector<MtzManager *>& managers, std::string filename,
		GraphMap properties)
{
	vector<vector<std::string> > filenames;
	time_t cputime;
	time(&cputime);

	for (int i = 0; i < managers.size(); i++)
	{
		filenames.push_back(vector<std::string>());

		for (int j = 0; j < 3; j++)
		{
			std::ostringstream stream;
			stream << cputime << "_bfactor_" << i << "_" << j;
			filenames[i].push_back(stream.str());
		}

		cputime++;
	}

	std::ostringstream results;

	properties["title"] = "Scale and B factor plot";
	properties["xTitle"] = "Resolution (1 / d^2)";
	properties["yTitle"] = "Ratio of intensities (image over reference)";
	properties["plotType"] = "line";

	properties["xMin"] = 0;
//	properties["xMax"] = 0.3;
	properties["yMin"] = 0;
	properties["yMax"] = 4;

	vector<double> bins;
	StatisticsManager::generateResolutionBins(50, 2.1, 18, &bins);

	vector<vector<double> > xs;
	vector<vector<double> > ys;

	MtzManager *referenceManager = MtzManager::getReferenceManager();
	/*
	double grad = 1000 / referenceManager->averageIntensity();
	referenceManager->applyScaleFactor(grad);
*/
	for (int i = 0; i < managers.size(); i++)
	{
		mtz = managers[i];
		mtz->excludeFromLogCorrelation();

		double correl = managers[i]->correlation(true);
		if (correl < 0.85)
			continue;

		 double gradient = managers[i]->gradientAgainstManager(
		 MtzManager::getReferenceManager());

	//	 double bFactor = 0;

	//	 managers[i]->bFactorAndScale(&gradient, &bFactor);
		 //	gradient = 1000 / managers[i]->averageIntensity();

		 	managers[i]->applyScaleFactor(gradient);


		vector<double> bins;
		StatisticsManager::generateResolutionBins(50, 1.6, 10, &bins);

		vector<double> x, y;

		for (int shell = 0; shell < bins.size() - 3; shell++)
		{
			double low = bins[shell];
			double high = bins[shell + 3];

			vector<ReflectionPtr> refReflections, imgReflections;

			mtz->findCommonReflections(referenceManager, imgReflections, refReflections,
			NULL);

            double weights = 0;
			double refMean = 0;
			double imgMean = 0;
			int count = 0;

			for (int i = 0; i < imgReflections.size(); i++)
			{
				if (!imgReflections[i]->anyAccepted())
					continue;

				if (imgReflections[i]->betweenResolutions(low, high))
				{
					weights += imgReflections[i]->meanPartiality();
					refMean += refReflections[i]->meanIntensity()
							* imgReflections[i]->meanPartiality();
					imgMean += imgReflections[i]->meanIntensity()
							* imgReflections[i]->meanPartiality();
					count++;
				}
			}

			refMean /= weights;
			imgMean /= weights;

			double ratio = refMean / imgMean;
			ratio = 1 / ratio;

			if (ratio != ratio)
				continue;

			x.push_back(1 / pow(low, 2));
			y.push_back(ratio);
		}

		xs.push_back(x);
		ys.push_back(y);
	}

	this->plot(filename, properties, xs, ys);

	xs.clear();
	ys.clear();
}

void GraphDrawer::plotPolarisation(vector<MtzPtr> mtzs)
{
    int count = 36;
    int divide = 360 / count;
    
    std::map<int, double> histogram;
    std::map<int, int> counts;
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        for (int j = 0; j < mtzs[i]->reflectionCount(); j++)
        {
            if (!mtzs[i]->reflection(j)->betweenResolutions(1.8, 1.4))
                continue;
            
            for (int k = 0; k < mtzs[i]->reflection(j)->millerCount(); k++)
            {
                MillerPtr miller = mtzs[i]->reflection(j)->miller(k);
                
                if (miller->getRawIntensity() < 0)
                    continue;
                
                vec hkl = new_vector(miller->getH(), miller->getK(), miller->getL());
                mtzs[i]->getMatrix()->multiplyVector(&hkl);
                
                double angle = cartesian_to_angle(hkl.h, hkl.k);
                double degrees = angle * 180 / M_PI;
                
                int category = (int)(degrees / divide) * divide;
                
                if (histogram.count(category) == 0)
                {
                    histogram[category] = 0;
                    counts[category] = 0;
                }
                
                histogram[category] += miller->getRawIntensity();
                counts[category]++;
            }
        }
    }
    
    vector<double> xs, ys;
    GraphMap graphMap;
    graphMap["yMin"] = 0;
    graphMap["yMax"] = 200;
    graphMap["title"] = "Average raw intensity vs angle on detector";
    graphMap["xTitle"] = "Angle on detector";
    graphMap["yTitle"] = "Average raw intensity";
    
    int num = 0;
    
    for (int i = 0; i <= 0; i++)
    {
        for (std::map<int, double>::iterator it = histogram.begin(); it != histogram.end(); it++)
        {
            if (i == 0)
                histogram[it->first] /= counts[it->first];
            
            xs.push_back(it->first + (i * 360));
            ys.push_back(histogram[it->first]);
            num++;
        }
    }
    
    std::cout << "Number of X values: " << num << std::endl;
    
    plot("polarisation", graphMap, xs, ys);
}

void GraphDrawer::partialityPlot(std::string filename, GraphMap properties, double maxRes)
{
    double correl = mtz->correlation();
    mtz->setActiveAmbiguity(1);
    double invCorrel = mtz->correlation();
    
    std::cout << "Ambiguity 0: " << correl << ", ambiguity 1: " << invCorrel << std::endl;
    
    if (correl > invCorrel)
        mtz->setActiveAmbiguity(0);
    
    double gradient = mtz->gradientAgainstManager(MtzManager::getReferenceManager());
    mtz->applyScaleFactor(gradient);
    
    properties["title"] = "Partiality plot " + mtz->getFilename();
	properties["xTitle"] = "Wavelength (Å)";
	properties["yTitle"] = "Fraction of merged data set (%)";
	properties["plotType"] = "fill";

	vector<double> resolutions;
	StatisticsManager::generateResolutionBins(0, maxRes, 4, &resolutions);

	vector<std::string> files;

	vector<Partial> partials;
	StatisticsManager::twoImagePartialityStatsWritten(&partials,
			&MtzManager::getReferenceManager(), &mtz);

	std::cout << "Total number of reflections in MTZ: " << partials.size() << std::endl;

	double maxPercentage = 1000;
	std::sort(partials.begin(), partials.end(), sortByPercentage);

	maxPercentage = partials[partials.size() - 10].percentage;

    if (maxPercentage < 500 || !std::isfinite(maxPercentage))
		maxPercentage = 500;

	maxPercentage = 500;

	std::sort(partials.begin(), partials.end(), sortByWavelength);

	double minX = partials[50].wavelength;
	double maxX = partials[partials.size() - 50].wavelength;
    double middle = partials[partials.size() / 2].wavelength;

    double wantedWidth = 0.08;
    
  //  std::cout << "minX: " << minX << ", maxX: " << maxX << ", middle: " << middle << std::endl;
    
	if (maxX - minX > wantedWidth)
	{
        maxX = middle + wantedWidth / 2;
		minX = middle - wantedWidth / 2;
	}

	int resCount = (int)resolutions.size();

	properties["xMin"] = minX;
	properties["xMax"] = maxX;
	properties["yMin2"] = 0;
	properties["yMax2"] = 5.0;
	properties["yMax"] = maxPercentage;

	vector<double> xs, xs2, xs3, xs4, ys, ys2, ys3, ys4, scatterX, scatterY;

    std::cout << std::setw(12) << "Low res" << std::setw(12) << "High res" << std::setw(12) << "Num refl." << std::endl;
    
    
	for (int shell = 0; shell < resCount - 1; shell++)
	{
		xs.clear();
		ys.clear();
        xs2.clear();
		ys2.clear();
        xs3.clear();
        ys3.clear();

		double lowRes = 1 / resolutions[shell];
		double highRes = 1 / resolutions[shell + 1];

        std::ofstream partLog;
        partLog.open("partiality_" + i_to_str(shell) + ".csv");
        partLog << "h,k,l,wavelength,partiality,percentage,intensity,resolution" << std::endl;
        
        int count = 0;
        
		for (int i = 0; i < partials.size(); i++)
		{
			if (partials[i].resolution > lowRes
					&& partials[i].resolution < highRes)
			{
                if (partials[i].wavelength != partials[i].wavelength
                    || partials[i].percentage != partials[i].percentage
                    || partials[i].partiality != partials[i].partiality)
                    continue;
                
                if (partials[i].percentage > maxPercentage)
                    continue;

                double h = partials[i].miller->getH();
                double k = partials[i].miller->getK();
                double l = partials[i].miller->getL();
                double wavelength = partials[i].wavelength;
                double partiality = partials[i].partiality;
                double percentage = partials[i].percentage;
                double intensity = partials[i].miller->getRawestIntensity();
                double resolution = partials[i].resolution;

                partLog << h << "," << k << "," << l << "," << wavelength << "," << partiality << "," <<
                percentage << "," << intensity << "," << resolution << std::endl;
                count++;
            }
		}
        
        partLog.close();

		std::cout << std::setw(12) << 1 / lowRes << std::setw(12) <<  1 / highRes << std::setw(12) << count
				<< std::endl;
	}
}

void GraphDrawer::plotPartialityStats(int h, int k, int l)
{
 //   double threshold = 100;
    GraphMap map;
    map["xMin"] = 0.0;
    map["xMax"] = 1;
    map["yMax"] = 3;
    map["xTitle"] = "Calculated partiality";
    map["yTitle"] = "Proportion of merged value";
    map["plotType"] = "point";
    
    vector<double> xs, ys;
    
    bool checkingHKL = !(h == 0 && k == 0 && l == 0);
    double reflId = Reflection::reflectionIdForCoordinates(h, k, l);
    
    for (int i = 0; i < mtz->reflectionCount(); i++)
    {
        bool okay = false;
        
        if (checkingHKL && mtz->reflection(i)->getReflId() == reflId)
        {
            okay = true;
            std::ostringstream logged;
            logged << "Found reflection " << h << " " << k << " " << l << std::endl;
            Logger::mainLogger->addStream(&logged);
        }
        
        if (okay == false)
            continue;
    
        for (int j = 0; j < mtz->reflection(i)->millerCount(); j++)
        {
            double partiality = mtz->reflection(i)->miller(j)->getPartiality();
            double percentage = mtz->reflection(i)->miller(j)->getRawIntensity();
            {
                xs.push_back(partiality);
                ys.push_back(percentage);
            }
        }
    }
    
    plot("partialityStats", map, xs, ys);
}

void GraphDrawer::plotOrientationStats(vector<MtzPtr> mtzs)
{
    CCP4SPG *spaceGroup = mtzs[0]->getLowGroup();
    UnitCellLattice lattice = UnitCellLattice(100, 100, 100, 90, 90, 90, spaceGroup->spg_num);
    
    std::ofstream pdbLog;
    pdbLog.open("plot.pdb");
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        MatrixPtr matrix = mtzs[i]->getMatrix()->getRotation();
        
        for (int j = 0; j < lattice.symOperatorCount(); j++)
        {
            MatrixPtr op = lattice.symOperator(j);
            
            MatrixPtr mat2 = matrix->copy();
            vec test = new_vector(1, 0, 0);
            mat2->preMultiply(*op);
            std::cout << mtzs[i]->getFilename() << "\t" << mat2->summary() << std::endl;
            mat2->multiplyVector(&test);
            scale_vector_to_distance(&test, 1);
            
            pdbLog << "HETATM";
            pdbLog << std::fixed;
            pdbLog << std::setw(5) << j << "                   ";
            pdbLog << std::setw(8) << std::setprecision(2) << test.h * 50;
            pdbLog << std::setw(8) << std::setprecision(2) << test.k * 50;
            pdbLog << std::setw(8) << std::setprecision(2) << test.l * 50;
            pdbLog << "                       O" << std::endl;
            
        }
    }
    
    pdbLog.close();

}

void GraphDrawer::plotSingleMillerFromMtzs(std::vector<MtzPtr> mtzs, int h, int k, int l)
{
    int index = Reflection::reflectionIdForCoordinates(h, k, l);
    
    std::ofstream csv;
    csv.open("miller_" + i_to_str(h) + "_" + i_to_str(k) + "_" + i_to_str(l) + ".csv");
    csv << "filename,delay,intensity,countingSigma,scale,partiality,phase,resolution" << std::endl;
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        MtzPtr mtz = mtzs[i];
        
        for (int j = 0; j < mtz->reflectionCount(); j++)
        {
            ReflectionPtr refl = mtz->reflection(j);
            
            long unsigned int thisID = refl->getReflId();
            
            if (thisID != index)
            {
                continue;
            }
            
            if (!refl->anyAccepted())
                continue;
            
            for (int k = 0; k < refl->millerCount(); k++)
            {
                MillerPtr miller = refl->miller(k);
                
                if (!miller->accepted())
                    continue;
                
                csv << mtz->getFilename() << "\t" << mtz->getTimeDelay() << "\t" << miller->intensity() << "\t" << miller->getCountingSigma() << "\t" << miller->getScale() << "\t" << miller->getPartiality() << "\t" << miller->getPhase() << "\t" << 1 / refl->getResolution() << std::endl;
            }
            
        }
    }

    csv.close();
}

void GraphDrawer::plotReflectionFromMtzs(std::vector<MtzPtr> mtzs, int h, int k, int l)
{
    if (!(h == 0 && k == 0 && l == 0))
    {
        plotSingleMillerFromMtzs(mtzs, h, k, l);
        return;
    }

    // if not, we need to plot all out to resolution
    
    double resol = FileParser::getKey("MAX_INTEGRATED_RESOLUTION", 3.0);
    std::vector<double> unitCell = mtzs[0]->getUnitCell();
    
    UnitCellLattice lattice = UnitCellLattice(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5], mtzs[0]->getLowGroup()->spg_num);
    
    int hMax, kMax, lMax;
    
    lattice.getMaxMillerIndicesForResolution(resol, &hMax, &kMax, &lMax);
    
    for (int i = -hMax; i < hMax; i++)
    {
        for (int j = -kMax; j < kMax; j++)
        {
            for (int k = -lMax; k < lMax; k++)
            {
                bool asu = ccp4spg_is_in_asu(mtzs[0]->getLowGroup(), i, j, k);
                
                if (asu)
                {
                    plotSingleMillerFromMtzs(mtzs, i, j, k);
                }
            }
        }
    }
}

void GraphDrawer::cutoutIntegrationAreas(std::vector<MtzPtr> mtzs, int h, int k, int l)
{
    int index = Reflection::reflectionIdForCoordinates(h, k, l);
    int windowPadding = 10;
    double threshold = FileParser::getKey("PNG_THRESHOLD", 1000.);
    
    PNGFile png = PNGFile("windows.png", 1000, 1000);
    
    int xGrid = 0;
    int yGrid = 0;
    
    std::ostringstream logged;
    
    logged << "Checking " << mtzs.size() << " mtz files." << std::endl;
    Logger::log(logged);
    
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        MtzPtr mtz = mtzs[i];
        ImagePtr image = mtz->getImagePtr();
        
        for (int j = 0; j < mtz->reflectionCount(); j++)
        {
            ReflectionPtr refl = mtz->reflection(j);
            
            long unsigned int thisID = refl->getReflId();
            
            if (thisID != index)
            {
                continue;
            }
            
            
            if (!refl->anyAccepted())
                continue;
            
            
            for (int k = 0; k < refl->millerCount(); k++)
            {
                MillerPtr miller = refl->miller(k);
                
                if (!miller->accepted())
                    continue;
                
                
                int xStart = xGrid * (windowPadding * 2 + 3) + 2;
                int yStart = yGrid * (windowPadding * 2 + 3) + 2;
                
                PanelPtr panel = Panel::panelForMiller(&*miller);
                Coord bestShift = panel->shiftForMiller(&*miller);
                
                double shiftedX = miller->getLastX() + bestShift.first;
                double shiftedY = miller->getLastY() + bestShift.second;
                
                logged << "Found Miller on " << mtz->getFilename() << " - will be plotting at (" << shiftedX << ", " << shiftedY << ")" << std::endl;
                Logger::log(logged);
                
                for (int s = -windowPadding; s < windowPadding + 1; s++)
                {
                    for (int t = -windowPadding; t < windowPadding + 1; t++)
                    {
                        int pngX = xStart + s + windowPadding;
                        int pngY = yStart + t + windowPadding;
                        
                        int localX = shiftedX + s;
                        int localY = shiftedY + t;
                        
                        double value = image->valueAt(localX, localY);
                        if (value < 0)
                            value = 0;
                        
                        unsigned char pixelValue = std::min(value, threshold) * 255 / threshold;
                        pixelValue = 255 - pixelValue;
                        
                        png.setPixelColour(pngX, pngY, pixelValue, pixelValue, pixelValue);
                    }
                }
                
                xGrid++;
                
                if (xGrid > 10)
                {
                    xGrid = 0;
                    yGrid++;
                }
            }
            
        }
    }
    
    png.writeImageOutput();
    
    logged << "Plotted Miller " << h << ", " << k << ", " << l << std::endl;
    Logger::log(logged);
}

GraphDrawer::~GraphDrawer()
{
// TODO Auto-generated destructor stub
}

