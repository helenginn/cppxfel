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

void GraphDrawer::resolutionStatsCSV(std::vector<MtzManager *>& managers)
{
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

        myMtz->scaleToMtz(this->mtz);

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
            
            if (fabs(wavelength - meanWavelength) > meanWavelength * 0.025)
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
    plotMap["xMax0"] = f_to_str(meanWavelength * 0.975);
    plotMap["xMin0"] = f_to_str(meanWavelength * 1.025);
    plotMap["yMax0"] = "250";
    plotMap["yMin0"] = "0";
    plotMap["xTitle0"] = "Ewald sphere wavelength (Ang)";
    plotMap["style0"] = "line";
    
    plotMap["xHeader1"] = "wavelength";
    plotMap["yHeader1"] = "partiality";
    plotMap["xMax1"] = f_to_str(meanWavelength * 0.975);
    plotMap["xMin1"] = f_to_str(meanWavelength * 1.025);
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
    
    mtz->scaleToMtz(MtzManager::getReferenceManager());
    
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

void GraphDrawer::plotOrientationStats(vector<MtzPtr> mtzs)
{
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
                bool asu = ccp4spg_is_in_asu(mtzs[0]->getSpaceGroup(), i, j, k);
                
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

