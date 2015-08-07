//
//  CommonCircle.h
//  cppxfel
//
//  Created by Helen Ginn on 04/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__CommonCircle__
#define __cppxfel__CommonCircle__

#include <stdio.h>
#include "parameters.h"
#include "Logger.h"
#include "Vector.h"

class Image;

class CommonCircle
{
private:
    Image *parentImage;
    double x; double y;
    double maxUnitCell;
    double pixMaxUnitCell;
    
    double coarseOuterMinX;
    double coarseOuterMinY;
    double coarseOuterMaxX;
    double coarseOuterMaxY;
    double coarseMinDistance;
    
    double rRadius;
    double pixRadius;
    double pixRadSqr;
    double pixTolerance;
    double pixRadSqrTolerance;
    std::vector<SpotPtr> spots;
    std::ostringstream logged;
    
    void signature(std::vector<double> *angles);
    void sendLog(LogLevel priority = LogLevelNormal);
public:
    CommonCircle(Image *image, double xPix, double yPix);
    bool spotIsOnCommonCircle(SpotPtr spot);
    bool addSpotIfOnCommonCircle(SpotPtr spot);
    void addSpotToCommonCircle(SpotPtr spot);
    int addSpotsIfOnCommonCircle(std::vector<SpotPtr> spots);
    bool spotIsTooCloseToOtherSpots(SpotPtr spot);
    vec beamToCentre();
    int makeSuccessful();
    void removeSuccess();
    
    double angleForSpot(SpotPtr spot);
    double swivelToHorizontalAxis(bool negative);
    static double reciprocalRadiusToPixelForImage(Image *image, double reciprocalRadius);
    static double twoThetaToReciprocalRadius(double twoTheta, double wavelength);
    static std::vector<double> radiusRangeForTwoThetaRange(double degStart, double degEnd, Image *image);
    static std::vector<double> coordsForPixRadiusAndAngle(double angle, double pixRadius);
    double getAngleSeparation();
    double reciprocalRadiusToTwoTheta(double repRadius, double wavelength);

    void orderSpotsInPairs(CommonCirclePtr otherCircle, std::vector<std::pair<SpotPtr, SpotPtr> > *orderedSpots);
    double similarityToCommonCircle(CommonCirclePtr otherCircle);
    
    void printSpots();
    std::string printResult();
    
    Image *getParentImage()
    {
        return parentImage;
    }
    
};

#endif /* defined(__cppxfel__CommonCircle__) */
