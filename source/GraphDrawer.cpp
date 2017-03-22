/*
 * GraphDrawer.cpp
 *
 *  Created on: 22 Oct 2014
 *      Author: helenginn
 */

#include "GraphDrawer.h"
#include "misc.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <boost/variant.hpp>
#include <map>
#include <sstream>
#include <tuple>
#include "Detector.h"
#include <fstream>
#include "Reflection.h"
#include "Miller.h"
#include "UnitCellLattice.h"
#include "FileParser.h"
#include "PNGFile.h"
#include "CSV.h"
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

void GraphDrawer::partialityPNGResolutionShell(std::string filename, double meanWavelength,
                                               std::vector<ReflectionPtr> refRefls,
                                               std::vector<ReflectionPtr> imageRefls,
                                               double minRes, double maxRes)
{
    CSVPtr csv = CSVPtr(new CSV(4, "wavelength", "intensity", "partiality", "percentage"));

    for (int i = 0; i < refRefls.size(); i++)
    {
        ReflectionPtr imageRefl = imageRefls[i];
        ReflectionPtr refRefl = refRefls[i];
        
        if (!(refRefl->betweenResolutions(minRes, maxRes)))
        {
            continue;
        }
        
        for (int j = 0; j < imageRefl->millerCount(); j++)
        {
            MillerPtr miller = imageRefl->miller(j);
            double wavelength = miller->getWavelength();
            
            if (fabs(wavelength - meanWavelength) > meanWavelength * 0.05)
            {
                continue;
            }
            
            double intensity = miller->getRawIntensity();
            double partiality = miller->getPartiality();
            double percentage = intensity / refRefl->meanIntensity() * 100;
            
            csv->addEntry(0, wavelength, intensity, partiality, percentage);
        }
    }
    
    std::ostringstream logged;
    logged << csv->entryCount() << " reflections between " << minRes << " and " << maxRes << std::endl;
    Logger::log(logged);
    
    std::string extendedFilename = filename + "_" + f_to_str(minRes, 3) + "_to_" + f_to_str(maxRes, 3) + "_partiality";
    
    std::map<std::string, std::string> plotMap;
    plotMap["filename"] = extendedFilename;
    plotMap["xHeader0"] = "wavelength";
    plotMap["yHeader0"] = "percentage";
    plotMap["xMax0"] = f_to_str(meanWavelength * 0.95);
    plotMap["xMin0"] = f_to_str(meanWavelength * 1.05);
    plotMap["yMax0"] = "250";
    plotMap["yMin0"] = "0";
    plotMap["xTitle0"] = "Ewald sphere wavelength (Ang)";
    plotMap["style0"] = "line";
    
    plotMap["xHeader1"] = "wavelength";
    plotMap["yHeader1"] = "partiality";
    plotMap["xMax1"] = f_to_str(meanWavelength * 0.95);
    plotMap["xMin1"] = f_to_str(meanWavelength * 1.05);
    plotMap["yMax1"] = "2.5";
    plotMap["yMin1"] = "0";
    plotMap["style1"] = "line";
    plotMap["colour1"] = "blue";
    
    csv->writeToFile(extendedFilename + ".csv");
    csv->plotPNG(plotMap);
}

void GraphDrawer::partialityPNG(MtzPtr mtz, double maxRes)
{
    MtzManager *reference = MtzManager::getReferenceManager();
    double wavelength = mtz->bestWavelength();
    
    if (maxRes == 0)
    {
        maxRes = 1 / mtz->maxResolution();
    }
    
    double correl = mtz->correlation();
    mtz->setActiveAmbiguity(1);
    double invCorrel = mtz->correlation();
    
    std::cout << "Ambiguity 0: " << correl << ", ambiguity 1: " << invCorrel << std::endl;
    
    if (correl > invCorrel)
        mtz->setActiveAmbiguity(0);
    
    double gradient = mtz->gradientAgainstManager(MtzManager::getReferenceManager());
    mtz->applyScaleFactor(gradient);
    
    const int binCount = 4;
    vector<double> resolutions;
    StatisticsManager::generateResolutionBins(0, maxRes, binCount, &resolutions);
    
    std::vector<ReflectionPtr> refRefls, imageRefls;
    mtz->findCommonReflections(reference, imageRefls, refRefls);
    
    for (int i = 0; i < binCount; i++)
    {
        partialityPNGResolutionShell(mtz->getFilename(), wavelength, refRefls, imageRefls, resolutions[i], resolutions[i + 1]);
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
    UnitCellLatticePtr lattice = UnitCellLattice::getMainLattice();
    
    std::ofstream pdbLog;
    pdbLog.open("plot.pdb");
    
    for (int i = 0; i < mtzs.size(); i++)
    {
        MatrixPtr matrix = mtzs[i]->getMatrix()->getRotation();
        
        for (int j = 0; j < lattice->symOperatorCount(); j++)
        {
            MatrixPtr op = lattice->symOperator(j);
            
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

    std::ostringstream logged;
    logged << "Plotted Miller " << h << ", " << k << ", " << l << std::endl;
    Logger::log(logged);
    
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
    
    UnitCellLatticePtr lattice = UnitCellLattice::getMainLattice();
    
    int hMax, kMax, lMax;
    
    lattice->getMaxMillerIndicesForResolution(resol, &hMax, &kMax, &lMax);
    
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
    int windowPadding = 20;
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
                
                double shiftedX, shiftedY;
                DetectorPtr detector = Detector::getMaster()->spotCoordForMiller(miller, &shiftedX, &shiftedY);
                
                if (!detector)
                {
                    continue;
                }
                
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

