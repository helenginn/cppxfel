//
//  Detector.cpp
//  cppxfel
//
//  Created by Helen Ginn on 26/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "Detector.h"
#include "Matrix.h"
#include "Miller.h"
#include "FileParser.h"
#include "Vector.h"
#include "Spot.h"
#include "CSV.h"

double Detector::mmPerPixel = 0;
DetectorPtr Detector::masterPanel = DetectorPtr();
int Detector::detectorActive = -1;

// MARK: initialisation and constructors

void Detector::initialiseZeros()
{
    unarrangedTopLeftX = 0;
    unarrangedTopLeftY = 0;
    unarrangedBottomRightX = 0;
    unarrangedBottomRightY = 0;
    unarrangedMidPointX = 0;
    unarrangedMidPointY = 0;
    
    memset(&arrangedMidPoint.h, 0, sizeof(arrangedMidPoint.h) * 3);
    memset(&slowDirection.h, 0, sizeof(slowDirection.h) * 3);
    memset(&fastDirection.h, 0, sizeof(fastDirection.h) * 3);
    memset(&rotationAngles.h, 0, sizeof(rotationAngles.h) * 3);
    
    parent = DetectorPtr();
    rotMat = MatrixPtr(new Matrix());
    changeOfBasisMat = MatrixPtr(new Matrix());
    
    if (mmPerPixel == 0)
    {
        mmPerPixel = FileParser::getKey("MM_PER_PIXEL", 0.11);
    }
}

Detector::Detector()
{
    initialiseZeros();
}

Detector::Detector(double distance, double beamX, double beamY)
{
    initialiseZeros();
    
    setSlowDirection(0, 1, 0);
    setFastDirection(1, 0, 0);
    
    double pixDistance = distance / mmPerPixel;
    
    arrangedMidPoint = new_vector(-beamX, -beamY, pixDistance);
}

Detector::Detector(DetectorPtr parent, Coord unarrangedTopLeft, Coord unarrangedBottomRight,
                   vec slowDir, vec fastDir, vec _arrangedTopLeft, bool lastIsMiddle)
{
    initialiseZeros();
    
    setParent(parent);
    setUnarrangedTopLeft(unarrangedTopLeft.first, unarrangedTopLeft.second);
    setUnarrangedBottomRight(unarrangedBottomRight.first, unarrangedBottomRight.second);
    setSlowDirection(slowDir.h, slowDir.k, slowDir.l);
    setFastDirection(fastDir.h, fastDir.k, fastDir.l);
    
    if (lastIsMiddle)
    {
        arrangedMidPoint = _arrangedTopLeft;
    }
    else
    {
        setArrangedTopLeft(_arrangedTopLeft.h, _arrangedTopLeft.k, _arrangedTopLeft.l);
    }
}

// MARK: Housekeeping and additional initialisation

void Detector::updateUnarrangedMidPoint()
{
    unarrangedMidPointX = (double)(unarrangedBottomRightX + unarrangedTopLeftX) / 2;
    unarrangedMidPointY = (double)(unarrangedBottomRightY + unarrangedTopLeftY) / 2;
}

void Detector::setArrangedTopLeft(double newX, double newY, double newZ)
{
    assert(directionSanityCheck());
    
    /* This is where top left offset currently is from mid point */
    vec currentArrangedBottomRightOffset;
    spotCoordToRelativeVec(unarrangedBottomRightX, unarrangedBottomRightY,
                           &currentArrangedBottomRightOffset);
    
    /* Now we create the arranged top left position (absolute) from CrystFEL (say). */
    vec currentArrangedTopLeft = new_vector(newX, newY, newZ);
    
    /* Now we add the movement required to establish the middle position (absolute) */
    add_vector_to_vector(&currentArrangedTopLeft, currentArrangedBottomRightOffset);

    /* Now we need to change basis of vectors to know where this is relative to parent */
    
    getParent()->getChangeOfBasis()->multiplyVector(&currentArrangedTopLeft);

    /* This change of basis should be the relative arranged mid point! maybe. */
    arrangedMidPoint = copy_vector(currentArrangedTopLeft);
}

bool Detector::directionSanityCheck()
{
    double slowLength = length_of_vector(slowDirection);
    double fastLength = length_of_vector(fastDirection);
    
    return (fabs(slowLength - 1.) < 0.05 && fabs(fastLength - 1.) < 0.05);
}

void Detector::updateCurrentRotation()
{
    MatrixPtr myRotation = MatrixPtr(new Matrix());
    myRotation->rotate(rotationAngles.h, rotationAngles.k, rotationAngles.l);
    MatrixPtr parentRot = (!isLUCA() ? getParent()->getRotation() : MatrixPtr(new Matrix()));
    
    rotMat = myRotation;
    rotMat->preMultiply(*parentRot);
    
    slowRotated = copy_vector(slowDirection);
    fastRotated = copy_vector(fastDirection);
    
    rotMat->multiplyVector(&slowRotated);
    rotMat->multiplyVector(&fastRotated);
    
    vec rotatedSlow = getRotatedSlowDirection();
    vec rotatedFast = getRotatedFastDirection();
    cross = cross_product_for_vectors(rotatedFast, rotatedSlow);
    
    double matComponents[9] = {rotatedFast.h, rotatedSlow.h, cross.h,
        rotatedFast.k, rotatedSlow.k, cross.k,
        rotatedFast.l, rotatedSlow.l, cross.l};
                                    
    MatrixPtr mat1 = MatrixPtr(new Matrix(matComponents));
    changeOfBasisMat = mat1->inverse3DMatrix();
    
    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->updateCurrentRotation();
    }
}

// Housekeeping calculations

vec Detector::getRotatedSlowDirection()
{
    return slowRotated;
}

vec Detector::getRotatedFastDirection()
{
    return fastRotated;
}

vec Detector::midPointOffsetFromParent()
{
    if (isLUCA())
    {
        return arrangedMidPoint;
    }
    
    vec parentMidPoint = getParent()->midPointOffsetFromParent();
    
    vec slow = getParent()->getRotatedSlowDirection();
    vec fast = getParent()->getRotatedFastDirection();
    
    multiply_vector(&slow, arrangedMidPoint.k);
    multiply_vector(&fast, arrangedMidPoint.h);
    
    vec modOffset = parentMidPoint;
    add_vector_to_vector(&modOffset, slow);
    add_vector_to_vector(&modOffset, fast);
    modOffset.l += arrangedMidPoint.l;
    
    return modOffset;
}

void Detector::cumulativeMidPoint(vec *cumulative)
{
    if (hasNoChildren())
    {
        *cumulative = new_vector(0, 0, 0);
    }
    
    add_vector_to_vector(cumulative, midPointOffsetFromParent());
}

double Detector::halfSlow()
{
    return (unarrangedBottomRightY - unarrangedTopLeftY) / 2;
}

double Detector::halfFast()
{
    return (unarrangedBottomRightX - unarrangedTopLeftX) / 2;
}

bool Detector::isAncestorOf(DetectorPtr detector)
{
    if (shared_from_this() == detector)
    {
        return true;
    }
    
    for (int i = 0; i < childrenCount(); i++)
    {
        if (getChild(i)->isAncestorOf(detector))
        {
            return true;
        }
    }
    
    return false;
}

// MARK: Spot to physical position where X-rays hit the detector

void Detector::spotCoordToRelativeVec(double unarrangedX, double unarrangedY,
                                      vec *arrangedPos)
{
    assert(hasNoChildren());
    
    double unArrangedOffsetX = unarrangedX - unarrangedMidPointX;
    double unArrangedOffsetY = unarrangedY - unarrangedMidPointY;
    
    vec slowOffset = getRotatedSlowDirection();
    *arrangedPos = getRotatedFastDirection();

    multiply_vector(&slowOffset, unArrangedOffsetY);
    multiply_vector(arrangedPos, unArrangedOffsetX);
    
    add_vector_to_vector(arrangedPos, slowOffset);
}

void Detector::spotCoordToAbsoluteVec(double unarrangedX, double unarrangedY,
                                      vec *arrangedPos)
{
    if (hasNoChildren())
    {
        spotCoordToRelativeVec(unarrangedX, unarrangedY, arrangedPos);
        
        /* Add my modifications onto my children's.
         * This could probably be improved for speed by
         * storing some pre-calculations automatically. */
        
        vec cumulative = new_vector(0, 0, 0);
        cumulativeMidPoint(&cumulative);
        
        add_vector_to_vector(&cumulative, *arrangedPos);
        
        *arrangedPos = copy_vector(cumulative);
    }
}

DetectorPtr Detector::findDetectorPanelForSpotCoord(double xSpot, double ySpot)
{
    if (hasChildren())
    {
        DetectorPtr probe;
        
        for (int i = 0; i < childrenCount(); i++)
        {
            DetectorPtr child = getChild(i);
            probe = child->findDetectorPanelForSpotCoord(xSpot, ySpot);
            
            if (probe)
            {
                break;
            }
        }
        
        return probe;
    }
    else
    {
        if (xSpot >= unarrangedTopLeftX && xSpot <= unarrangedBottomRightX &&
            ySpot >= unarrangedTopLeftY && ySpot <= unarrangedBottomRightY)
        {
            return shared_from_this();
        }
        else
        {
            return DetectorPtr();
        }
    }
}

DetectorPtr Detector::spotToAbsoluteVec(SpotPtr spot, vec *arrangedPos)
{
    if (spot->getDetector())
    {
        spotCoordToAbsoluteVec(spot->getRawXY().first, spot->getRawXY().second, arrangedPos);
        return spot->getDetector();
    }
    
    DetectorPtr detector = findDetectorAndSpotCoordToAbsoluteVec(spot->getRawXY().first, spot->getRawXY().second, arrangedPos);
    spot->setDetector(detector);
    
    return detector;
}

DetectorPtr Detector::findDetectorAndSpotCoordToAbsoluteVec(double xSpot, double ySpot,
                       vec *arrangedPos)
{
    DetectorPtr detector = findDetectorPanelForSpotCoord(xSpot, ySpot);
    
    if (detector)
    {
        detector->spotCoordToAbsoluteVec(xSpot, ySpot, arrangedPos);
    }
    
    return detector;
}

// MARK: ray trace to spot position

void Detector::intersectionToSpotCoord(vec absoluteVec, double *xSpot, double *ySpot)
{
    // make the assumption that absoluteVec is indeed in the plane of the detector
    
    /* Need the cumulative midpoint for this panel */
    
    vec cumulative;
    cumulativeMidPoint(&cumulative);
    
    /* Need offset from centre of the panel */
    take_vector_away_from_vector(cumulative, &absoluteVec);
    
    /* how many slow directions vs how many fast directions...? */
    /* i.e. change of basis... */
    
    changeOfBasisMat->multiplyVector(&absoluteVec);
    
    /* Now to add back the unarranged mid point */
    absoluteVec.h += unarrangedMidPointX;
    absoluteVec.k += unarrangedMidPointY;
    
    *xSpot = absoluteVec.h;
    *ySpot = absoluteVec.k;
}

void Detector::intersectionWithRay(vec ray, vec *intersection)
{
    vec cumulative;
    cumulativeMidPoint(&cumulative);
    
    double t = cumulative.h * cross.h + cumulative.k * cross.k + cumulative.l * cross.l;
    t /= (cross.h * ray.h + cross.k * ray.k + cross.l * ray.l);
    
    *intersection = new_vector(t * ray.h, t * ray.k, t * ray.l);
}

bool Detector::spotCoordWithinBounds(double xSpot, double ySpot)
{
    xSpot -= unarrangedMidPointX;
    ySpot -= unarrangedMidPointY;
    
    bool intersects = (fabs(xSpot) <= halfFast() && fabs(ySpot) <= halfSlow());
    
    return intersects;
}

DetectorPtr Detector::detectorForRayIntersection(vec ray, vec *intersection)
{
    if (hasChildren())
    {
        DetectorPtr probe;
        
        for (int i = 0; i < childrenCount(); i++)
        {
            DetectorPtr child = getChild(i);
            probe = child->detectorForRayIntersection(ray, intersection);
            
            if (probe)
            {
                break;
            }
        }
        
        return probe;
    }
    else
    {
        // Time to find the ray intersection
        // point on plane: arrangedMidPoint
        // plane normal: cross
        
        intersectionWithRay(ray, intersection);
        
        double xSpot, ySpot;
        
        intersectionToSpotCoord(*intersection, &xSpot, &ySpot);
        
        bool intersects = spotCoordWithinBounds(xSpot, ySpot);
        
        if (intersects)
        {
            return shared_from_this();
        }
        else
        {
            return DetectorPtr();
        }
    }
}

DetectorPtr Detector::intersectionForMiller(MillerPtr miller, vec *intersection)
{
    vec hkl = miller->getTransformedHKL();
    double imageWavelength = miller->getImage()->getWavelength();
    double tmp = hkl.k;
    hkl.k = -hkl.h;
    hkl.h = -tmp;
    hkl.l += 1 / imageWavelength;
    
    DetectorPtr probe = miller->getDetector();
    
    if (!probe)
    {
        probe = detectorForRayIntersection(hkl, intersection);
        
        if (probe)
        {
            miller->setDetector(probe);
        }
    }
    else
    {
        probe->intersectionWithRay(hkl, intersection);
    }
    
    return probe;
}

DetectorPtr Detector::spotCoordForMiller(MillerPtr miller, double *xSpot, double *ySpot)
{
    vec intersection;
    DetectorPtr probe = intersectionForMiller(miller, &intersection);
    
    if (probe)
    {
        probe->intersectionToSpotCoord(intersection, xSpot, ySpot);
    }
    
    return probe;
}

DetectorPtr Detector::spotCoordForRayIntersection(vec ray, double *xSpot, double *ySpot)
{
    vec intersection;
    
    DetectorPtr probe = detectorForRayIntersection(ray, &intersection);
    
    if (probe)
    {
        probe->intersectionToSpotCoord(intersection, xSpot, ySpot);
    }
    
    return probe;
}

// MARK: Keeping track of added Millers

void Detector::addMillerCarefully(MillerPtr miller)
{
    millerMutex.lock();
    
    millers.push_back(miller);
    
    millerMutex.unlock();
    
    if (!isLUCA())
    {
        getParent()->addMillerCarefully(miller);
    }
}

// MARK: read/write detector files

std::string indents(int indentCount)
{
    std::string indent = "";
    
    for (int i = 0; i < indentCount; i++)
    {
        indent += "  ";
    }
    
    return indent;
}

void Detector::applyRotations()
{
    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->applyRotations();
    }
    
    updateCurrentRotation();
}

std::string Detector::writeGeometryFile(int indentCount)
{
    std::ostringstream output;
    
    output << indents(indentCount) << "panel " << getTag() << std::endl;
    output << indents(indentCount) << "{" << std::endl;
    output << indents(indentCount + 1) << "min_fs = " << unarrangedTopLeftX << std::endl;
    output << indents(indentCount + 1) << "min_ss = " << unarrangedTopLeftY << std::endl;
    output << indents(indentCount + 1) << "max_fs = " << unarrangedBottomRightX << std::endl;
    output << indents(indentCount + 1) << "max_ss = " << unarrangedBottomRightY << std::endl;
    output << indents(indentCount + 1) << "fs_x = " << fastDirection.h << std::endl;
    output << indents(indentCount + 1) << "fs_y = " << fastDirection.k << std::endl;
    output << indents(indentCount + 1) << "fs_z = " << fastDirection.l << std::endl;
    output << indents(indentCount + 1) << "ss_x = " << slowDirection.h << std::endl;
    output << indents(indentCount + 1) << "ss_y = " << slowDirection.k << std::endl;
    output << indents(indentCount + 1) << "ss_z = " << slowDirection.l << std::endl;
    output << indents(indentCount + 1) << "midpoint_x = " << arrangedMidPoint.h << std::endl;
    output << indents(indentCount + 1) << "midpoint_y = " << arrangedMidPoint.k << std::endl;
    output << indents(indentCount + 1) << "midpoint_z = " << arrangedMidPoint.l << std::endl;
    output << indents(indentCount + 1) << "alpha = " << rotationAngles.h << std::endl;
    output << indents(indentCount + 1) << "beta = " << rotationAngles.k << std::endl;
    output << indents(indentCount + 1) << "gamma = " << rotationAngles.l << std::endl;

    for (int i = 0; i < childrenCount(); i++)
    {
        output << std::endl << getChild(i)->writeGeometryFile(indentCount + 1);
    }
    
    output << indents(indentCount) << "}" << std::endl;
    
    return output.str();
}

// MARK: Generate resolution histogram

CSVPtr Detector::resolutionHistogram()
{
    double minResolution = FileParser::getKey("MIN_INTEGRATED_RESOLUTION", 0.);
    double maxResolution = FileParser::getKey("MAX_INTEGRATED_RESOLUTION", 1.4);
    double approxWavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.);
    
    if (approxWavelength == 0.)
    {
        std::ostringstream logged;
        logged << "Cannot create resolution histogram without an approximate value for INTEGRATION_WAVELENGTH. Please provide!" << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelNormal, true);
    }
    
    std::vector<double> size = FileParser::getKey("DETECTOR_SIZE", std::vector<double>(0, 2));
    int xDim = size[0];
    int yDim = size[1];
    
    CSVPtr csv = CSVPtr(new CSV(0));
    csv->setupHistogram(minResolution, maxResolution, -1000, "Resolution", 1, "numPixels");
    
    DetectorPtr master = Detector::getMaster();
    const vec beamAxis = new_vector(0, 0, 1);
    
    for (int i = 0; i < yDim; i++)
    {
        for (int j = 0; j < xDim; j++)
        {
            vec intersection;
            DetectorPtr aDetector = master->findDetectorAndSpotCoordToAbsoluteVec(j, i, &intersection);
            
            if (!aDetector)
            {
                continue;
            }
            
            double twoTheta = angleBetweenVectors(beamAxis, intersection);
            double theta = twoTheta / 2;
            
            double resolution = approxWavelength / (2 * sin(theta));
            
            csv->addOneToFrequency(resolution, "numPixels");
        }
    }
    
    csv->writeToFile("resolutions.csv");
    return csv;
}

// MARK: Scoring functions

double Detector::millerScoreWrapper(void *object)
{
    return static_cast<Detector *>(object)->millerScore();
}

double Detector::millerScore(bool ascii)
{
    if (!millers.size())
    {
        logged << "No integrated images so cannot use reflection predictions as a scoring function." << std::endl;
        sendLog();
        
        return 0;
    }
    
    Miller::refreshMillerPositions(millers);
    
    if (!xShifts.size())
    {
        yShifts.clear();
        
        for (int i = 0; i < millerCount(); i++)
        {
            MillerPtr myMiller = miller(i);
            
            if (!myMiller || !myMiller->reachesThreshold())
                continue;
            
            float *shiftX = myMiller->getXShiftPointer();
            float *shiftY = myMiller->getYShiftPointer();
            
            xShifts.push_back(shiftX);
            yShifts.push_back(shiftY);
        }
    }

    std::vector<double> xs, ys;
    
    CSV *csv;
    
    if (ascii)
    {
        csv = new CSV(2, "x", "y");
    }
    
    for (int i = 0; i < xShifts.size(); i++)
    {
        xs.push_back(*xShifts[i]);
        ys.push_back(*yShifts[i]);
        
        if (ascii)
        {
            csv->addEntry(0, xs[i], ys[i]);
        }
    }

    if (ascii)
    {
        csv->plotColumns(0, 1);
        delete csv;
    }
    
    double stdevX = standard_deviation(&xs, NULL, 0);
    double stdevY = standard_deviation(&ys, NULL, 0);
    double aScore = stdevX + stdevY;
    
    return aScore;
}

