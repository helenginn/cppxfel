/*
 * Spot.cpp
 *
 *  Created on: 30 Dec 2014
 *      Author: helenginn
 */

#include "Spot.h"
#include <iostream>
#include <cmath>
#include <cfloat>
#include "Image.h"
#include <algorithm>
#include "Vector.h"
#include <fstream>
#include "FileParser.h"
#include "CSV.h"
#include "Shoebox.h"
#include "Detector.h"

vector<vector<double> > Spot::probe;

double Spot::maxResolution = 0;
double Spot::minIntensity = 0;
double Spot::minCorrelation = -3;
bool Spot::checkRes = false;

Spot::Spot(ImagePtr image)
{
        // TODO Auto-generated constructor stub
        parentImage = image;
    angleDetectorPlane = 0;
    setAngle = false;
    checked = false;
    successfulCommonLines = 0;
    correctedX = -1;
    correctedY = -1;
    _isBeamCentre = false;
    rejected = false;
    _isFake = false;
        _estimatedVec = new_vector(0, 0, 0);
    x = 0;
    y = 0;
    height = FileParser::getKey("IMAGE_SPOT_PROBE_HEIGHT", 100);
    background = FileParser::getKey("IMAGE_SPOT_PROBE_BACKGROUND", 0);
    int probePadding = FileParser::getKey("IMAGE_SPOT_PROBE_PADDING", 1);
    length = probePadding * 2 + 1;
    backgroundPadding = FileParser::getKey("IMAGE_SPOT_PROBE_BG_PADDING", probePadding) * 2 + 1;
    storedRadius = 1 / image->getWavelength();

        makeProbe(height, background, length, backgroundPadding);
}

Spot::~Spot()
{

}

bool Spot::isAcceptable(ImagePtr image)
{
        int length = (int)probe.size();
        int tolerance = (length - 1) / 2;

        for (int i = x - tolerance; i < x + tolerance; i++)
        {
                for (int j = y - tolerance; j < y + tolerance; j++)
                {
                        if (!image->accepted(i, j))
                                return false;
                }
        }

        return true;
}

double Spot::focusOnNearbySpot(double maxShift, double trialX, double trialY, int round)
{
    double focusedX = trialX;
    double focusedY = trialY;

    if (maxShift > 0)
    {
        this->getParentImage()->focusOnAverageMax(&focusedX, &focusedY, maxShift);
    }

    setXY(focusedX, focusedY);

    double pixelIntensity = this->getParentImage()->valueAt(focusedX, focusedY);

    if (pixelIntensity < minIntensity)
    {
        return false;
    }

    if (checkRes)
    {
        double resol = this->resolution();

        if (resol > (1. / maxResolution)) return false;
        if (resol < 0) return false;
    }

    std::vector<double> probeIntensities, realIntensities;
    probeIntensities.reserve(backgroundPadding * backgroundPadding);
    realIntensities.reserve(backgroundPadding * backgroundPadding);

    int padding = (backgroundPadding - 1) / 2;

    for (int j = -padding; j < padding + 1; j++)
    {
        for (int i = -padding; i < padding + 1; i++)
        {
            if (!(this->getParentImage()->accepted(focusedX + i, focusedY + j)))
            {
                return false;
            }

            int probeX = i + padding;
            int probeY = j + padding;

            double probeIntensity = probe[probeX][probeY];
            double realIntensity = this->getParentImage()->valueAt(focusedX + i, focusedY + j);

            if (!(probeIntensity == probeIntensity && realIntensity == realIntensity))
                continue;

            probeIntensities.push_back(probeIntensity);
            realIntensities.push_back(realIntensity);
        }
    }

    double correlation = correlation_between_vectors(&probeIntensities, &realIntensities);

    return correlation;
}

void Spot::makeProbe(int height, int background, int length, int backPadding)
{
    if (probe.size() > 0)
    {
        return;
    }

    probeMutex.lock();

    if (probe.size() > 0)
    {
        probeMutex.unlock();
        return;
    }

        logged << "Setting global spot parameters..." << std::endl;
        sendLog();

        checkRes = FileParser::hasKey("MAX_INTEGRATED_RESOLUTION");
        maxResolution = FileParser::getKey("MAX_INTEGRATED_RESOLUTION", 2.0);
        minIntensity = FileParser::getKey("IMAGE_MIN_SPOT_INTENSITY", 100.);
        minCorrelation = FileParser::getKey("IMAGE_MIN_CORRELATION", 0.7);


    if (backPadding < length)
        backPadding = length;

        if (backPadding % 2 == 0)
                backPadding++;

    if (length % 2 == 0)
        length++;

    int backgroundThickness = (backPadding - length) / 2;

        int size = (backPadding - 1) / 2;

        int centre = size;

        for (int i = 0; i < backPadding; i++)
        {
                probe.push_back(vector<double>());

                for (int j = 0; j < backPadding; j++)
                {
                        probe[i].push_back(background);
                }
        }

        for (int i = backgroundThickness; i < length + backgroundThickness; i++)
        {
                for (int j = backgroundThickness; j < length + backgroundThickness; j++)
                {
                        int from_centre_i = i - centre;
                        int from_centre_j = j - centre;

                        double distance_from_centre = sqrt(
                                        from_centre_i * from_centre_i
                                                        + from_centre_j * from_centre_j);

                        double fraction = distance_from_centre / (double) (size + 1);

            // further out the fraction, the lower the value
                        if (fraction > 1)
            {
                                probe[i][j] = background;
            }
                        else
            {
                                probe[i][j] = (1 - fraction) * height + background;
            }
                }
        }

    probeMutex.unlock();
}

void Spot::setXY(double x, double y)
{
        this->x = x;
        this->y = y;
}

double Spot::resolution()
{
    vec spotVec = this->estimatedVector();

    return length_of_vector(spotVec);
}

double Spot::angleFromSpotToCentre(double centreX, double centreY)
{
    double beamX = getParentImage()->getBeamX();
    double beamY = getParentImage()->getBeamY();

    vec beamCentre = new_vector(beamX, beamY, 0);
    vec circleCentre = new_vector(centreX, centreY, 0);
    vec circleCentreToBeam = vector_between_vectors(circleCentre, beamCentre);

    return angleInPlaneOfDetector(centreX, centreY, circleCentreToBeam);
}

double Spot::angleInPlaneOfDetector(double centreX, double centreY, vec upBeam)
{
    if (centreX == 0 && centreY == 0)
    {
        centreX = getParentImage()->getBeamX();
        centreY = getParentImage()->getBeamY();
    }

    vec centre = new_vector(centreX, centreY, 0);
    vec spotVec = new_vector(getX(), getY(), 0);
    vec spotVecFromCentre = vector_between_vectors(centre, spotVec);

    angleDetectorPlane = angleBetweenVectors(upBeam, spotVecFromCentre);

    return angleDetectorPlane;
}

void Spot::setUpdate()
{
    getX(true);
    getY(true);
}

Coord Spot::getXY()
{
    if (isBeamCentre())
    {
        x = getParentImage()->getBeamX();
        y = getParentImage()->getBeamY();

        return std::make_pair(x, y);
    }


    vec arrangedPos = new_vector(FLT_MAX, FLT_MAX, FLT_MAX);

    DetectorPtr detector = Detector::getMaster()->spotToAbsoluteVec(shared_from_this(), &arrangedPos);

    return std::make_pair(arrangedPos.h, arrangedPos.k);
}

double Spot::getX(bool update)
{
    if (update || correctedX == -1)
    {
        Coord shift = getXY();
        correctedX = shift.first;
    }

    return correctedX;
}

double Spot::getY(bool update)
{
    if (update || correctedY == -1)
    {
        Coord shift = getXY();
        correctedY = shift.second;
    }

    return correctedY;
}

Coord Spot::getRawXY()
{
    return std::make_pair(x, y);
}

void Spot::recentreInWindow(int windowPadding)
{
        if (windowPadding == 0)
        {
                windowPadding = backgroundPadding + 4;
        }

        ImagePtr thisImage = getParentImage();
        double x = getRawX();
        double y = getRawY();

        recentreInWindow(thisImage, &x, &y, windowPadding);

        setXY(x, y);
}

void Spot::recentreInWindow(ImagePtr thisImage, double *x, double *y, int windowPadding)
{
    std::vector<double> xPositions, yPositions, weights;

    for (int i = -windowPadding; i < windowPadding + 1; i++)
    {
        for (int j = -windowPadding; j < windowPadding + 1; j++)
        {
            int myX = *x + i;
            int myY = *y + j;

            if (thisImage->accepted(myX, myY))
                continue;

            int pixelIntensity = thisImage->valueAt(myX, myY);

            if (pixelIntensity <= 10)
                continue;

            xPositions.push_back(myX);
            yPositions.push_back(myY);
            weights.push_back(pixelIntensity);
        }
    }

    double newX = weighted_mean(&xPositions, &weights);
    double newY = weighted_mean(&yPositions, &weights);

    if (weights.size() > 0)
    {
                *x = newX;
                *y = newY;
    }
}

bool Spot::isOnSameLineAsSpot(SpotPtr otherSpot, double toleranceDegrees)
{
    double tolerance = toleranceDegrees * M_PI / 180;

    double thisAngle = angleInPlaneOfDetector();
    double otherAngle = otherSpot->angleInPlaneOfDetector();

    if (thisAngle < 0)
        thisAngle += M_PI;

    if (otherAngle < 0)
        otherAngle += M_PI;

    return (thisAngle > otherAngle - tolerance && thisAngle < otherAngle + tolerance);
}

void Spot::writeDatFromSpots(std::string filename, std::vector<SpotPtr> spots)
{
    bool enableCSVs = FileParser::getKey("ENABLE_IMAGE_CSVS", false);

    if (!enableCSVs)
        return;

    int count = 0;

    CSV csv = CSV(9, "h", "k", "l", "angle", "1", "1", "x", "y", "1");

    for (int i = 0; i < spots.size(); i++)
    {
        vec estimated = spots[i]->estimatedVector();

        csv.addEntry(1000, estimated.h, estimated.k, estimated.l, spots[i]->angleInPlaneOfDetector(), 1., 1., spots[i]->getRawXY().first, spots[i]->getRawXY().second, 1.);

    //    dat << "0\t0\t0\t" << spots[i]->angleInPlaneOfDetector() << "\t1\t1\t"
    //    << spots[i]->x << "\t" << spots[i]->y
    //    << "\t1" << std::endl;

        count++;
    }

    csv.writeToFile(filename);
    std::ostringstream logged;
    logged << "Written " << spots.size() << " spots to " << filename << std::endl;
    Logger::log(logged);

    //csv.plotColumns(6, 7);

//    dat.close();
}

vec Spot::estimatedVector()
{
    if (_isBeamCentre)
    {
        vec reciprocalPos = new_vector(0, 0, 0);

        return reciprocalPos;
    }

    if (hasDetector())
    {
        getDetector()->spotCoordToAbsoluteVec(x, y, &_estimatedVec);
    }
    else
    {
                DetectorPtr det  = Detector::getMaster()->spotToAbsoluteVec(shared_from_this(), &_estimatedVec);

        if (!det)
        {
            return new_vector(0, 0, 0);
        }
    }

    scale_vector_to_distance(&_estimatedVec, storedRadius);
    _estimatedVec.l -= storedRadius;

    return _estimatedVec;
}

double Spot::integrate()
{
    int foregroundLength = FileParser::getKey("SHOEBOX_FOREGROUND_PADDING",
                                              SHOEBOX_FOREGROUND_PADDING);
    int neitherLength = FileParser::getKey("SHOEBOX_NEITHER_PADDING",
                                           SHOEBOX_NEITHER_PADDING);
    int backgroundLength = FileParser::getKey("SHOEBOX_BACKGROUND_PADDING",
                                              SHOEBOX_BACKGROUND_PADDING);
    bool shoeboxEven = FileParser::getKey("SHOEBOX_MAKE_EVEN", false);

    ShoeboxPtr shoebox = ShoeboxPtr(new Shoebox(MillerPtr()));
    shoebox->simpleShoebox(foregroundLength, neitherLength, backgroundLength, shoeboxEven);

    float counting = 0;
    Coord rawXY = getRawXY();
    intensity = getParentImage()->intensityAt(rawXY.first, rawXY.second, shoebox, &counting);

    if (intensity != intensity)
    {
        intensity = 0;
    }

    return intensity;
}

std::string Spot::spotLine()
{
    if (_isBeamCentre)
    {
        return std::string("");
    }

    std::ostringstream line;

    std::string fakeFlag = isFake() ? "fake" : "";

    line << x << "\t" << y << "\t" << fakeFlag << std::endl;

    return line.str();
}

bool Spot::isSameAs(SpotPtr spot2)
{
    return (spot2->getX() == getX() && spot2->getY() == getY());
}

double Spot::closeToSecondSpot(SpotPtr spot2, double squareMinDistance)
{
    if (this == &*spot2)
        return false;

    vec reciprocalSpot1 = estimatedVector();
    vec reciprocalSpot2 = spot2->estimatedVector();

    take_vector_away_from_vector(reciprocalSpot1, &reciprocalSpot2);

    double square = length_of_vector_squared(reciprocalSpot2);


    return (square < squareMinDistance);
}

void Spot::spotsAndVectorsToResolution(double lowRes, double highRes, std::vector<SpotPtr> spots, std::vector<SpotVectorPtr> spotVectors, std::vector<SpotPtr> *lowResSpots, std::vector<SpotVectorPtr> *lowResSpotVectors)
{
    double minResol = lowRes == 0 ? 0 : (1 / lowRes);
    double maxResol = highRes == 0 ? FLT_MAX : (1 / highRes);

    for (int i = 0; i < spots.size(); i++)
    {
        if (spots[i]->resolution() < maxResol && spots[i]->resolution() > minResol)
        {
            lowResSpots->push_back(spots[i]);
        }
    }

    for (int i = 0; i < spotVectors.size(); i++)
    {
        SpotPtr firstSpot = spotVectors[i]->getFirstSpot();
        SpotPtr secondSpot = spotVectors[i]->getSecondSpot();

        if ((firstSpot->resolution() < maxResol && secondSpot->resolution() < maxResol)
            && (firstSpot->resolution() > minResol && secondSpot->resolution() > minResol))
        {
            lowResSpotVectors->push_back(spotVectors[i]);
        }
    }
}

void Spot::addToMask(int *mask, int width, int height)
{
    int padding = (backgroundPadding - 1) / 2 + 2;

    for (int j = -padding; j < padding + 1; j++)
    {
        for (int i = -padding; i < padding + 1; i++)
        {
            double xPos = getRawXY().first;
            double yPos = getRawXY().second;

            int xShifted = xPos + i;
            int yShifted = yPos + j;

            int position = width * yShifted + xShifted;

            if (position > width * height || position < 0)
            {
                continue;
            }

            mask[position] = 1;
        }
    }
}
