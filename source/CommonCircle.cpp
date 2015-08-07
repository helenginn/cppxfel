//
//  CommonCircle.cpp
//  cppxfel
//
//  Created by Helen Ginn on 04/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "CommonCircle.h"
#include "Image.h"
#include "Vector.h"
#include "FileParser.h"
#include "Miller.h"
#include "Logger.h"

void CommonCircle::signature(std::vector<double> *angles)
{
    for (int i = 0; i < spots.size(); i++)
    {
        SpotPtr spot = spots[i];
        double angle = spot->angleFromSpotToCentre(x, y);
        angles->push_back(angle);
    }
}

vec CommonCircle::beamToCentre()
{
    vec beam = new_vector(parentImage->getBeamX(), parentImage->getBeamY(), 0);
    vec centre = new_vector(x, y, 0);
    vec beamToCentre = vector_between_vectors(beam, centre);
    
    return beamToCentre;
}

double CommonCircle::swivelToHorizontalAxis(bool negative)
{
    vec thisBeamToCentre = beamToCentre();
    vec thatBeamToCentre = new_vector(negative ? -1 : 1, 0, 0);
    
    double angle = angleBetweenVectors(thisBeamToCentre, thatBeamToCentre);
    
    return angle;
}

std::vector<double> CommonCircle::radiusRangeForTwoThetaRange(double degStart, double degEnd, Image *image)
{
    double radStart = degStart * M_PI / 180;
    double radEnd = degEnd * M_PI / 180;
    double wavelength = image->getWavelength();
    
    double rRadiusStart = twoThetaToReciprocalRadius(radStart, wavelength);
    double rRadiusEnd = twoThetaToReciprocalRadius(radEnd, wavelength);
    
    std::vector<double> range;
    range.push_back(reciprocalRadiusToPixelForImage(image, rRadiusEnd));
    range.push_back(reciprocalRadiusToPixelForImage(image, rRadiusStart));
    
    return range;
}

double CommonCircle::twoThetaToReciprocalRadius(double twoTheta, double wavelength)
{
    double theta = twoTheta / 2;
    double cosTheta = cos(theta);
    double radius = cosTheta / wavelength;
    
    return radius;
}

double CommonCircle::getAngleSeparation()
{
    logged << "Reciprocal radius: " << rRadius << std::endl;
    double twoTheta = reciprocalRadiusToTwoTheta(rRadius, this->parentImage->getWavelength());
    sendLog(LogLevelDebug);
    return twoTheta;
}

double CommonCircle::reciprocalRadiusToTwoTheta(double repRadius, double wavelength)
{
    double ewaldRadius = 1 / wavelength;
    double theta = acos(repRadius / (2 * ewaldRadius));
    
    return 2 * theta;
}

double CommonCircle::reciprocalRadiusToPixelForImage(Image *image, double reciprocalRadius)
{
    double wavelength = image->getWavelength();
    double distance = image->getDetectorDistance();
    double mmPerPixel = image->getMmPerPixel();
    
    double newH = 0;
    double newK = sqrt((4 * pow(reciprocalRadius, 2) - pow(reciprocalRadius, 4) * pow(wavelength, 2)) / 4);
    double newL = 0 - pow(reciprocalRadius, 2) * wavelength / 2;
    
    vec hkl = new_vector(newH, newK, newL);
    
    double x_mm = (hkl.k * distance / (1 / wavelength + hkl.l));
    double y_mm = (hkl.h * distance / (1 / wavelength + hkl.l));
    
    double x_coord = x_mm / mmPerPixel;
    double y_coord = y_mm / mmPerPixel;
    
    return sqrt(pow(x_coord, 2) + pow(y_coord, 2));
}


CommonCircle::CommonCircle(Image *image, double xPix, double yPix)
{
    parentImage = image;
    pixTolerance = FileParser::getKey("PIXEL_TOLERANCE", 3.0);
    maxUnitCell = FileParser::getKey("MAX_UNIT_CELL", 500.0);
    
    double rMaxUnitCell = 1 / maxUnitCell;
    pixMaxUnitCell = reciprocalRadiusToPixelForImage(image, rMaxUnitCell);
    
    x = xPix;
    y = yPix;
    
    double beamX = image->getBeamX();
    double beamY = image->getBeamY();
    
    double diamX = 2 * (xPix - x) + x;
    double diamY = 2 * (yPix - y) + y;
    
    double detector_distance = image->getDetectorDistance();
    double wavelength = image->getWavelength();
    
    vec crystalRealSpace = new_vector(beamX, beamY, -detector_distance);
    vec crystalReciprocal = new_vector(0, 0, -1 / wavelength);
    vec detectorPos = new_vector(diamX, diamY, 0);
    vec scatterVec = vector_between_vectors(crystalRealSpace, detectorPos);
    
    scale_vector_to_distance(&scatterVec, 1 / wavelength);
    vec originToPoint = copy_vector(crystalReciprocal);
    add_vector_to_vector(&originToPoint, scatterVec);
    
    rRadius = length_of_vector(originToPoint) / 2;
    pixRadius = sqrt(pow(xPix - beamX, 2) + pow(yPix - beamY, 2));
    pixRadSqr = pow(pixRadius, 2);
    pixRadSqrTolerance = pow(pixRadius + pixTolerance, 2) - pow(pixRadius, 2);
    
    coarseOuterMinX = x - pixRadius - pixTolerance * 2;
    coarseOuterMaxX = x + pixRadius + pixTolerance * 2;
    
    coarseOuterMinY = y - pixRadius - pixTolerance * 2;
    coarseOuterMaxY = y + pixRadius + pixTolerance * 2;
    
    coarseMinDistance = pixRadius / sqrt(2) - pixTolerance * 2;
    
}

bool CommonCircle::spotIsOnCommonCircle(SpotPtr spot)
{
    double xSpot = spot->getX();

    if (xSpot > coarseOuterMaxX || xSpot < coarseOuterMinX)
        return false;
    
    if (fabs(xSpot - x) < coarseMinDistance)
        return false;
    
    double ySpot = spot->getY();

    if (ySpot > coarseOuterMaxY || ySpot < coarseOuterMinY)
        return false;

    if (fabs(ySpot - y) < coarseMinDistance)
        return false;
    
    double dSqr = pow(xSpot - x, 2) + pow(ySpot - y, 2);
    
    return (fabs(dSqr - pixRadSqr) < pixRadSqrTolerance);
}

void CommonCircle::addSpotToCommonCircle(SpotPtr spot)
{
    spots.push_back(spot);
}

bool CommonCircle::spotIsTooCloseToOtherSpots(SpotPtr spot)
{
    for (int i = 0; i < spots.size(); i++)
    {
        double x = spots[i]->getX();
        double y = spots[i]->getY();
        double xSpot = spot->getX();
        double ySpot = spot->getY();
        
        double distance = sqrt(pow(x - xSpot, 2) + pow(y - ySpot, 2));
        if (distance < pixMaxUnitCell)
            return true;
    }
    
    return false;
}

bool CommonCircle::addSpotIfOnCommonCircle(SpotPtr spot)
{
    if (spotIsOnCommonCircle(spot) && !spotIsTooCloseToOtherSpots(spot))
    {
        addSpotToCommonCircle(spot);
        return true;
    }
    
    return false;
}

int CommonCircle::addSpotsIfOnCommonCircle(std::vector<SpotPtr> newSpots)
{
    int total = 0;
    
    for (int i = 0; i < newSpots.size(); i++)
    {
        total += addSpotIfOnCommonCircle(newSpots[i]);
    }
    
    return total;
}

void CommonCircle::printSpots()
{
    logged << "Common circle pixel radius " << pixRadius << " for position (" << x << ", " << y << ")" << std::endl;
    logged << "Image: " << parentImage->getFilename() << std::endl;
    
    for (int i = 0; i < spots.size(); i++)
    {
        logged << spots[i]->getX() << "\t" << spots[i]->getY() << "\t" << spots[i]->angleFromSpotToCentre(x, y) << std::endl;
    }
    
    logged << std::endl;
    
    sendLog(LogLevelNormal);
}

void CommonCircle::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

static bool spotIsGreaterAngleThanSpot(std::pair<SpotPtr, double> one, std::pair<SpotPtr, double> two)
{
    return (one.second > two.second);
}

double CommonCircle::angleForSpot(SpotPtr spot)
{
    return spot->angleFromSpotToCentre(x, y);
}

// rewrite this
void CommonCircle::orderSpotsInPairs(CommonCirclePtr otherCircle, std::vector<std::pair<SpotPtr, SpotPtr> > *orderedSpots)
{
    std::vector<std::pair<SpotPtr, double> > theseAnglesAndSpots, thoseAnglesAndSpots;
    
    for (int i = 0; i < spots.size(); i++)
        theseAnglesAndSpots.push_back(std::make_pair(spots[i], angleForSpot(spots[i])));
    
    for (int i = 0; i < otherCircle->spots.size(); i++)
        thoseAnglesAndSpots.push_back(std::make_pair(otherCircle->spots[i], otherCircle->angleForSpot(otherCircle->spots[i])));
    
    std::sort(theseAnglesAndSpots.begin(), theseAnglesAndSpots.end(), spotIsGreaterAngleThanSpot);
    std::sort(thoseAnglesAndSpots.rbegin(), thoseAnglesAndSpots.rend(), spotIsGreaterAngleThanSpot);
 
    if (theseAnglesAndSpots.size() != thoseAnglesAndSpots.size())
    {
        double thisGreaterThanOther = (theseAnglesAndSpots.size() > thoseAnglesAndSpots.size());
        
        std::vector<std::pair<SpotPtr, double> > bigger = thisGreaterThanOther ? theseAnglesAndSpots : thoseAnglesAndSpots;
        std::vector<std::pair<SpotPtr, double> > smaller = thisGreaterThanOther ? thoseAnglesAndSpots : theseAnglesAndSpots;
        std::vector<std::pair<SpotPtr, double> > reducedBigger;
        
        for (int i = 0; i < smaller.size(); i++)
        {
            std::pair<SpotPtr, double> smallerOne = smaller[i];
            double minDifference = 1000;
            int minJ = -1;
            
            for (int j = 0; j < bigger.size(); j++)
            {
                double diff = (bigger[j].second + smallerOne.second); // smaller one is reversed
                if (minDifference > diff)
                {
                    minDifference = diff;
                    minJ = j;
                }
            }
            
            if (minJ != -1)
            {
                reducedBigger.push_back(bigger[minJ]);
                bigger.erase(bigger.begin() + minJ);
            }
        }
        
        theseAnglesAndSpots = thisGreaterThanOther ? reducedBigger : smaller;
        thoseAnglesAndSpots = thisGreaterThanOther ? smaller : reducedBigger;
    }
    
    for (int i = 0; i < theseAnglesAndSpots.size(); i++)
    {
        orderedSpots->push_back(std::make_pair(theseAnglesAndSpots[i].first, thoseAnglesAndSpots[i].first));
    }
}

double CommonCircle::similarityToCommonCircle(CommonCirclePtr otherCircle)
{
    if (fabs(this->pixRadius - otherCircle->pixRadius) > pixTolerance)
        return 1000;
    
    std::vector<double> thisSignature, otherSignature;
    this->signature(&thisSignature);
    otherCircle->signature(&otherSignature);
    bool thisGreaterThanOther = false;
    
    std::sort(thisSignature.begin(), thisSignature.end());
    std::sort(otherSignature.rbegin(), otherSignature.rend()); // reverse
    
    if (thisSignature.size() != otherSignature.size())
    {
        thisGreaterThanOther = (thisSignature.size() > otherSignature.size());
        
        std::vector<double> bigger = thisGreaterThanOther ? thisSignature : otherSignature;
        std::vector<double> smaller = thisGreaterThanOther ? otherSignature : thisSignature;
        std::vector<double> reducedBigger;
        
        for (int i = 0; i < smaller.size(); i++)
        {
            double smallerOne = smaller[i];
            double minDifference = 1000;
            int minJ = -1;
            
            for (int j = 0; j < bigger.size(); j++)
            {
                double diff = (bigger[j] + smallerOne); // smaller one is reversed
                if (minDifference > diff)
                {
                    minDifference = diff;
                    minJ = j;
                }
            }
            
            if (minJ != -1)
            {
                reducedBigger.push_back(bigger[minJ]);
                bigger.erase(bigger.begin() + minJ);
            }
        }
        
        thisSignature = thisGreaterThanOther ? reducedBigger : smaller;
        otherSignature = thisGreaterThanOther ? smaller : reducedBigger;
    }
    
    double sumSquares = 0;
    
    for (int i = 0; i < thisSignature.size(); i++)
    {
        double thisAngle = thisSignature[i];
        double otherAngle = otherSignature[i];
        
        double squareDifference = pow(thisAngle + otherAngle, 2); // otherangle is reversed
        sumSquares += squareDifference;
    }
    
    sumSquares /= thisSignature.size();
    
    return sqrt(sumSquares);
}

std::string CommonCircle::printResult()
{
    std::ostringstream result;
    result << "circle " << x << " " << y << std::endl;
    
    return result.str();
}

int CommonCircle::makeSuccessful()
{
    int total = 0;
    
    for (int i = 0; i < spots.size(); i++)
    {
        int currentSuccessful = spots[i]->successfulLineCount();
        
        spots[i]->addSuccessfulCommonLineCount();
        
        if (currentSuccessful == 0)
            total++;
    }
    
    return total;
}

void CommonCircle::removeSuccess()
{
    for (int i = 0; i < spots.size(); i++)
    {
        spots[i]->deleteSuccessfulCommonLineCount();
    }
}
