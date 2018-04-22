//
//  ProfileFit.cpp
//  cppxfel
//
//  Created by Helen Mary Elizabeth Duyvesteyn on 04/03/2018.
//  Copyright © 2018 Structural Biology. All rights reserved.
//
#include "ProfileFit.h"
//#include "parameters.h"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
//#include <boost/thread/thread.hpp>
#include "Image.h"
#include "Vector.h"
//#include "Miller.h"
#include "Logger.h"
#include "FileParser.h"
//#include "Shoebox.h"
//#include "Hdf5ManagerProcessing.h"

#include "misc.h"
#include "CSV.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>




void ProfileFit::makeShoebox(ShoeboxPtr inBox)
{
    
    int foregroundLength = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",
                                              SHOEBOX_FOREGROUND_PADDING);
    int neitherLength = FileParser::getKey("SHOEBOX_NEITHER_PADDING",
                                           SHOEBOX_NEITHER_PADDING);
    int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",
                                              SHOEBOX_BACKGROUND_PADDING);
    bool shoeboxEven = FileParser::getKey("SHOEBOX_MAKE_EVEN", false);
    
    inBox->simpleShoebox(foregroundLength, neitherLength, backgroundLength, shoeboxEven);
}



double ProfileFit::estimateFWHM(double stdevI)
{
    double fWHM = 0;
    fWHM = 2 * sqrt(2 * log(2)) * stdevI;
    return fWHM;
}


void ProfileFit::addShoeboxes(ShoeboxPtr shoeTmp)
{
    int slowSide1 = 0;
    int slowSide2 = 0;
    int fastSide1 = 0;
    int fastSide2 = 0;
    int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",SHOEBOX_BACKGROUND_PADDING);
    int limitB = (backgroundLength * 2 + 1);
    
    std::cout << "Attempting to add together shoeboxes!" << std::endl;
    //shoeMaster->sideLengths(&slowSide1, &fastSide1); 22
    //shoeTmp->sideLengths(&slowSide2, &fastSide2); 23
    //if (( slowSide1 == slowSide2) && ( fastSide1 == fastSide2))
    //{
    for (int i = 0; i < limitB; i++)
    {
        for (int j = 0; j < limitB; j++)
        {
            double val = (*shoeTmp)[i][j];
            if (val > 1)
            {
                (*this)[i][j] += val;
            }
            else
            {
                (*this)[i][j] = (*this)[i][j];
            }
            
        }
    }
    //}
    //else
    //{
    //std::cout << "Error comparing shoeboxes." << slowSide1 << "rather than " << slowSide2 << std::endl;
    //}
}

void ProfileFit::subtractfromShoebox(double toSubtract)
{
    //int slowSide = 0;
    //int fastSide = 0;
    int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",
                                              SHOEBOX_BACKGROUND_PADDING);
    int limitB = (backgroundLength * 2);
    //this->sideLengths(&slowSide, &fastSide);
    std::cout << "Attempting subtract value from shoebox!" << std::endl;
    for (int i = 0; i < limitB; i++)
    {
        for (int j = 0; j < limitB; j++)
        {
            if ((*this)[i][j] > 1)
            {
                std::cout << "Attempting to subtract:" << toSubtract << std::endl;
                (*this)[i][j] -= toSubtract ;
            }
	}        
    }
}

std::vector<double> listBox(ShoeboxPtr inBox)
{
    std::vector<double> outList;
    int slowSide = 0;
    int fastSide = 0;
    inBox ->sideLengths(&slowSide, &fastSide);
    for (int i = 0; i < slowSide; i++)
    {
        for (int j = 0; j < fastSide; j++)
        {
            outList.push_back((*inBox)[i][j]);
        }
    }
    
    return outList;
}


//

ShoeboxPtr ProfileFit::intensitySearchGrid(int xCoord, int yCoord, ImagePtr im, std::string fName, int spotNum, double *backGround)
{
    //Returns vector with intensities.
    int foregroundLength = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",
                                              SHOEBOX_FOREGROUND_PADDING);
    int neitherLength = FileParser::getKey("SHOEBOX_NEITHER_PADDING",
                                           SHOEBOX_NEITHER_PADDING);
    int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",
                                              SHOEBOX_BACKGROUND_PADDING);
    bool shoeboxEven = FileParser::getKey("SHOEBOX_MAKE_EVEN", false);
    
    double startX = (xCoord - backgroundLength);
    double startY = (yCoord - backgroundLength);
    
    ShoeboxPtr intensityBox = ShoeboxPtr(new Shoebox(MillerPtr()));
    //double intensityBox;
    double bgSummation = 0;
    int bgCount = 0;
    int x = 0;
    int y = 0;
    intensityBox->simpleShoebox(foregroundLength, neitherLength, backgroundLength, shoeboxEven);
    intensityBox->centre(&xCoord, &yCoord);
    std::cout << "Checking intensity box." << std::endl;
    //intensityBox->printShoebox();
    
    int slowSide = 0;
    int fastSide = 0;
    
    intensityBox->sideLengths(&slowSide, &fastSide);
    
    float counting = 0;
    int count = 0;
    
    CSVPtr csvspot = CSVPtr(new CSV(3,"x", "y", "Intensity"));
    for (int i = 0; i < slowSide; i++)
    {
        int panelPosX = (i - x) + startX ;
        for (int j = 0; j < fastSide; j++)
        {
            int panelPosY = (j - y) + startY ;
            double intensity = im->valueAt(panelPosX, panelPosY); //WRONG!
            if ((*intensityBox)[i][j] == 1)
            {
                (*intensityBox)[i][j] = intensity;
                count += 1;
                csvspot->addEntry(3, (double)i, (double)j, (double)intensity);
            }
            else if ((*intensityBox)[i][j] == -1)
            {
                bgSummation += intensity;
                bgCount += 1;
            }
            
        }
    }
    
    *backGround += (double)(bgSummation / bgCount);
    double tmpBg = (bgSummation / bgCount);
    std::cout<<"Searched grid. Writing to file:" << "\t" << fName << std::endl;
    std::string rootFilename = getBaseFilename(fName);
    
    csvspot->writeToFile("profs_" + rootFilename + "_s" + i_to_str(spotNum) + ".csv");
    return intensityBox;
}


//Complete spot.


void ProfileFit::correctSpotCentring(double *xCoord, double *yCoord, ImagePtr img)
{
    bool even = false;
    int shoeBoxF = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",SHOEBOX_FOREGROUND_PADDING);
    if (shoeBoxF % 2 == 0)
        even = true;
    
    std::cout << "orig:" << xCoord << "," << yCoord  << std::endl;
    
    //get Mss
    int mss = FileParser::getKey("METROLOGY_SEARCH_SIZE", -1);
    std::cout << "Metrology:" << mss << std::endl;
    img->focusOnAverageMax(xCoord, yCoord, mss, 0, even);
    
    //double correctedX = xCoord;
    //double correctedY = yCoord;
    std::cout << "after:" << xCoord << "," << yCoord  << std::endl;
    
}


void ProfileFit::calculateImageProfile(ImagePtr img)
{
    std::cout << "Hello! Calculating image profile." << std::endl;
    std::vector<double> imageSpotStdevs;
    std::vector<double> imageSpotMeans;
    std::vector<double> imageSpotFWHMs;
    std::vector<int> imageSpotNums;
    std::vector<double> imageSpotIntensities;
    std::vector<int> imageSpotXs;
    std::vector<int> imageSpotYs;
    std::vector<double> imageTotalVector;
    
    
    int foregroundLength = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",
                                              SHOEBOX_FOREGROUND_PADDING);
    int neitherLength = FileParser::getKey("SHOEBOX_NEITHER_PADDING",
                                           SHOEBOX_NEITHER_PADDING);
    int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",
                                              SHOEBOX_BACKGROUND_PADDING);
    bool shoeboxEven = FileParser::getKey("SHOEBOX_MAKE_EVEN", false);
    int limitB = (backgroundLength * 2);
    
    simpleShoebox(foregroundLength, neitherLength, backgroundLength, shoeboxEven);
    
    CSVPtr csv = CSVPtr(new CSV(7,"SpotNum", "x", "y", "Intensity", "MeanIntensity","Stdeviation","FWHM"));
    std::string fName = img->getFilename();
    std::string rootFilename = getBaseFilename(fName);
    std::cout << "Setting filename to: " << "/t" << fName << std::endl;
    
    std::cout << "Spot count for image:" << img->spotCount() << std::endl;
    
    double backGround = 0;
    int spotCount = 0;
    
    for (int j = 0; j < img->spotCount(); j++) //
    {
        //ReflectionPtr ref = mtz->reflection(j);
        std::cout << "Looking at spot number #" << j << std::endl;
        //int xCoord = mtz->miller->getCorrectedX();
        //int yCoord = mtz->miller->getCorrectedY();
        double xCoord = img->spot(j)->getRawX(); //\no geom, as opp to getX.
        double yCoord = img->spot(j)->getRawY();
        std::cout << "calcorig:" << xCoord << "," << yCoord  << std::endl;
        correctSpotCentring(&xCoord, &yCoord, img); //Correct centring.
        std::cout << "calcafter:" << xCoord << "," << yCoord  << std::endl;
        if (img->accepted(xCoord, yCoord))
        {
            
            std::vector<double> intensityList;
            //ShoeboxPtr intensityBox ;
            ShoeboxPtr intensityBox = ShoeboxPtr(new Shoebox(MillerPtr()));
            intensityBox = intensitySearchGrid(xCoord, yCoord, img, fName, j, &backGround);
            std::cout << "IntensityBox now completed." << std::endl;
            
            addShoeboxes(intensityBox);
            
            std::cout << "Finished search grid!" << std::endl;
            double imageSpotStdev = 0;
            double imageSpotMean = 0;
            double imageSpotFWHM = 0;
            int imageSpotIntensity = 0;
            
            intensityList = listBox(intensityBox);
            
            imageSpotIntensity = std::accumulate(intensityList.rbegin(), intensityList.rend(), 0);
            
            imageSpotStdev = standard_deviation(&intensityList);
            std::cout << "Raw standard deviation of" << imageSpotStdev << std::endl;
            imageSpotMean = weighted_mean(&intensityList);
            imageSpotFWHM = estimateFWHM(imageSpotStdev);
            std::cout << "Reading out standard deviations:" << imageSpotStdev << std::endl;
            std::cout << "Reading out fwhms:" << imageSpotFWHM << std::endl;
            imageSpotIntensities.push_back(imageSpotIntensity);
            imageSpotNums.push_back(j);
            //Summary of indiv spot profiles for one image.
            csv->addEntry(7, (double)j, (double)xCoord, (double)yCoord, (double)imageSpotIntensity,(double)imageSpotMean,(double)imageSpotStdev,(double)imageSpotFWHM);
            std::cout << "Background:" << backGround << std::endl;
            spotCount +=1;
            
        }
    }
    
    csv->writeToFile("profsTot_" + rootFilename + ".csv");
    
    double totalAverageBg = (backGround / spotCount);
    
    CSVPtr csv2 = CSVPtr(new CSV(3,"RelativeX","RelativeY","SummedIntensity"));
    
    subtractfromShoebox(totalAverageBg);
    
   printShoebox();
    
    for (int i = 0; i < limitB; i++)
    {
        for (int j = 0; j < limitB; j++)
        {
            double val = (*this)[i][j];
            csv2->addEntry(3, (double)i, (double)j, (double)val);
        }
    }
    
    csv2->writeToFile("spotsTotIntensity_" + rootFilename + ".csv");
    std::cout << "Calculated average background of:" << "\t" << totalAverageBg << std::endl;
    
    std::cout << "Printing intensities for file" << fName << std::endl;
}



///////////////////

/////DESTRUCTOR

//////////////////

//ProfileFit::~ProfileFit()




/*

ProfileFit::ProfileFit()
{
    int logInt = FileParser::getKey("VERBOSITY_LEVEL", 0);
    Logger::mainLogger->changePriorityLevel((LogLevel)logInt);
    ProfileFit::shoeBoxF = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",SHOEBOX_FOREGROUND_PADDING);
    ProfileFit::shoeBoxB = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",SHOEBOX_FOREGROUND_PADDING);
    ProfileFit::shoeBoxN = FileParser::getKey("SHOEBOX_NEITHER_PADDING",SHOEBOX_FOREGROUND_PADDING);
    //loadedImage = false;
    //hasFitted = false;
}
{
    std::cout << "Deallocating ProfileFitter." << std::endl;
    
}
 */

