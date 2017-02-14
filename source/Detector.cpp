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
#include "misc.h"
#include "CSV.h"

double Detector::mmPerPixel = 0;
DetectorPtr Detector::masterPanel = DetectorPtr();
int Detector::detectorActive = -1;
bool Detector::noisy = false;
ImagePtr Detector::drawImage = ImagePtr();
DetectorType Detector::detectorType = DetectorTypeMPCCD;
bool Detector::enabledNudge = false;
int Detector::specialImageCounter = 0;
std::mutex Detector::threadMutex;

// MARK: initialisation and constructors

void Detector::initialiseZeros()
{
    gain = 1;
    
    unarrangedTopLeftX = 0;
    unarrangedTopLeftY = 0;
    unarrangedBottomRightX = 0;
    unarrangedBottomRightY = 0;
    unarrangedMidPointX = 0;
    unarrangedMidPointY = 0;
    enabledNudge = false;
    
    memset(&arrangedMidPoint.h, 0, sizeof(arrangedMidPoint.h) * 3);
    memset(&slowDirection.h, 0, sizeof(slowDirection.h) * 3);
    memset(&fastDirection.h, 0, sizeof(fastDirection.h) * 3);
    memset(&rotationAngles.h, 0, sizeof(rotationAngles.h) * 3);
    memset(&nudgeTranslation.h, 0, sizeof(nudgeTranslation.h) * 3);
    memset(&nudgeRotation.h, 0, sizeof(nudgeRotation.h) * 3);
    memset(&quickMidPoint.h, 0, sizeof(quickMidPoint.h) * 3);
    memset(&quickAngles.h, 0, sizeof(quickAngles.h) * 3);

    parent = DetectorPtr();
    rotMat = MatrixPtr(new Matrix());
    nudgeMat = MatrixPtr(new Matrix());
    changeOfBasisMat = MatrixPtr(new Matrix());
    fixedBasis = MatrixPtr(new Matrix());
    mustUpdateMidPoint = true;
    
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
    double pixDistance = distance / mmPerPixel;
    arrangedMidPoint = new_vector(-beamX, -beamY, pixDistance);

    setSlowDirection(0, 1, 0);
    setFastDirection(1, 0, 0);
}

Detector::Detector(DetectorPtr parent, vec arrangedMiddle, std::string tag)
{
    initialiseZeros();
    
    setSlowDirection(0, 1, 0);
    setFastDirection(1, 0, 0);
    
    arrangedMidPoint = arrangedMiddle;
    setParent(parent);
    setTag(tag);
}

Detector::Detector(DetectorPtr parent, Coord arrangedTopLeft, Coord arrangedBottomRight,
                   double angle, double offsetX, double offsetY, double gain)
{
    initialiseZeros();
    setParent(parent);
    
    double rad = -angle * M_PI / 180;
    
    vec slowAxis = new_vector(0, 1, 0);
    vec fastAxis = new_vector(1, 0, 0);
    
    MatrixPtr mat = MatrixPtr(new Matrix());
    mat->rotate2D(rad);
    
    mat->multiplyVector(&slowAxis);
    mat->multiplyVector(&fastAxis);
    
    setSlowDirection(slowAxis.h, slowAxis.k, slowAxis.l);
    setFastDirection(fastAxis.h, fastAxis.k, fastAxis.l);
    
    updateCurrentRotation();

    arrangedMidPoint = new_vector((arrangedTopLeft.first + arrangedBottomRight.first) / 2,
                                  (arrangedTopLeft.second + arrangedBottomRight.second) / 2, 0);
    
    vec topLeft = new_vector(arrangedTopLeft.first, arrangedTopLeft.second, 0);
    vec bottomRight = new_vector(arrangedBottomRight.first, arrangedBottomRight.second, 0);
    vec reverseOffset = new_vector(offsetX, offsetY, 0);
    
    take_vector_away_from_vector(arrangedMidPoint, &topLeft);
    take_vector_away_from_vector(arrangedMidPoint, &bottomRight);
    
    mat->multiplyVector(&topLeft);
    mat->multiplyVector(&bottomRight);
    
    add_vector_to_vector(&topLeft, arrangedMidPoint);
    add_vector_to_vector(&bottomRight, arrangedMidPoint);
    
    add_vector_to_vector(&topLeft, reverseOffset);
    add_vector_to_vector(&bottomRight, reverseOffset);
    
    setUnarrangedTopLeft(std::min(topLeft.h, bottomRight.h), std::min(topLeft.k, bottomRight.k));
    setUnarrangedBottomRight(std::max(topLeft.h, bottomRight.h), std::max(topLeft.k, bottomRight.k));
}

Detector::Detector(DetectorPtr parent)
{
    initialiseZeros();
    
    setParent(parent);
}

void Detector::initialise(Coord unarrangedTopLeft, Coord unarrangedBottomRight,
                   vec slowDir, vec fastDir, vec _arrangedTopLeft, bool lastIsMiddle, bool ghost)
{
    setUnarrangedTopLeft(unarrangedTopLeft.first, unarrangedTopLeft.second);
    setUnarrangedBottomRight(unarrangedBottomRight.first, unarrangedBottomRight.second);
    setSlowDirection(slowDir.h, slowDir.k, slowDir.l);
    setFastDirection(fastDir.h, fastDir.k, fastDir.l);
    
    if (ghost)
    {
        calculateChangeOfBasis(&fastDir, &slowDir, fixedBasis);
        setFastDirection(1, 0, 0);
        setSlowDirection(0, 1, 0);
    }
    
    if (!isLUCA())
    {
        MatrixPtr parentBasis = getParent()->getChangeOfBasis();
        fixedBasis->multiply(*parentBasis->inverse3DMatrix());
    }
    
    rotateAxisRecursive();
    
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

void Detector::lockNudges()
{
    getMaster()->rotateAxisRecursive(true);
}

void Detector::recalculateChangeOfBasis(vec *fastAxis, vec *slowAxis)
{
    changeOfBasisMat = calculateChangeOfBasis(fastAxis, slowAxis, invBasisMat);
}

MatrixPtr Detector::calculateChangeOfBasis(vec *fastAxis, vec *slowAxis, MatrixPtr inv)
{
    cross = cross_product_for_vectors(*fastAxis, *slowAxis);
    double matComponents[9] = {fastAxis->h, slowAxis->h, cross.h,
        fastAxis->k, slowAxis->k, cross.k,
        fastAxis->l, slowAxis->l, cross.l};
    MatrixPtr mat = MatrixPtr(new Matrix(matComponents));
    memcpy(inv->components, mat->components, 16 * sizeof(double));
    return inv->inverse3DMatrix();
}

void Detector::addToBasisChange(vec angles, MatrixPtr chosenMat)
{
    if (!chosenMat)
    {
        chosenMat = changeOfBasisMat;
    }
    
    vec genFast = new_vector(1, 0, 0);
    vec genSlow = new_vector(0, 1, 0);
    
    rotMat->setIdentity();
    rotMat->rotate(angles.h, angles.k, angles.l);
    
    rotMat->multiplyVector(&genFast);
    rotMat->multiplyVector(&genSlow);
    
    MatrixPtr invBasis = MatrixPtr(new Matrix());
    calculateChangeOfBasis(&genFast, &genSlow, invBasis);
    chosenMat->multiply(*invBasis);
    invBasisMat = chosenMat->inverse3DMatrix();
}


void Detector::rotateAxisRecursive(bool fix)
{
    MatrixPtr tempNewChange = MatrixPtr(new Matrix());
    
    changeOfBasisMat->copyComponents(fixedBasis);
    
    if (fix)
    {
        tempNewChange->copyComponents(changeOfBasisMat);
    }
    
    if (!isLUCA())
    {
        std::lock_guard<std::mutex> lg(threadMutex);
        
        MatrixPtr mat = getParent()->getChangeOfBasis();
        changeOfBasisMat->multiply(*mat);
    }
    
    addToBasisChange(rotationAngles, changeOfBasisMat);
    addToBasisChange(rotationAngles, tempNewChange);
    
    if (enabledNudge)
    {
        vec rebasedNudge = new_vector(0, 0, 0);
        mustUpdateMidPoint = true;
        midPointOffsetFromParent(true, &rebasedNudge, fix);

        addToBasisChange(rebasedNudge);
        addToBasisChange(rebasedNudge, tempNewChange);
    }
    
    if (fix)
    {
        fixedBasis->copyComponents(tempNewChange);
        
        memset(&rotationAngles.h, 0, sizeof(rotationAngles.h) * 3);
        memset(&nudgeTranslation.h, 0, sizeof(nudgeTranslation.h) * 3);
        memset(&nudgeRotation.h, 0, sizeof(nudgeRotation.h) * 3);
    }
    
    slowRotated = slowDirection;
    changeOfBasisMat->multiplyVector(&slowRotated);
    fastRotated = fastDirection;
    changeOfBasisMat->multiplyVector(&fastRotated);
    
    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->rotateAxisRecursive(fix);
    }

}

void Detector::updateCurrentRotation()
{
    rotateAxisRecursive();
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

vec Detector::midPointOffsetFromParent(bool useParent, vec *angles, bool resetNudge)
{
    if (!enabledNudge)
    {
        vec myRawMidPoint = rawMidPointOffsetFromParent(useParent);
        
        return myRawMidPoint;
    }
    
    if (enabledNudge && !mustUpdateMidPoint && useParent && !resetNudge)
    {
        if (angles != NULL)
        {
            *angles = quickAngles;
        }
        
        return quickMidPoint;
    }
    
    vec parentedMidPoint = rawMidPointOffsetFromParent(true);
    
    vec cross = copy_vector(parentedMidPoint);
    scale_vector_to_distance(&cross, 1);
    vec horiz = new_vector(1, 0, -cross.h / cross.l);
    scale_vector_to_distance(&horiz, 1);
    
    vec vert = cross_product_for_vectors(cross, horiz);
    
    double matComponents[9] = {horiz.h, vert.h, cross.h,
        horiz.k, vert.k, cross.k,
        horiz.l, vert.l, cross.l};
    
    MatrixPtr changeToXYZ = MatrixPtr(new Matrix(matComponents));
    MatrixPtr changeToXYZInv = changeToXYZ->inverse3DMatrix();
    
    vec xyzMovement = copy_vector(nudgeTranslation);
    xyzMovement.l = 0;
    
    changeToXYZ->multiplyVector(&xyzMovement);
    double distanceFromOrigin = length_of_vector(parentedMidPoint);
    
    if (angles != NULL)
    {
        vec myAngles = copy_vector(nudgeRotation);
//        myAngles.l = 0;
        
        double tanHoriz = nudgeTranslation.k / distanceFromOrigin;
        double tanVert = nudgeTranslation.h / distanceFromOrigin;
        myAngles.h -= atan(tanHoriz);
        myAngles.k += atan(tanVert);
        /* End not buggy [B] */
        
        changeToXYZ->multiplyVector(&myAngles);
        
        *angles = copy_vector(myAngles);
        quickAngles = *angles;
    }

    vec shifted = copy_vector(xyzMovement);

    add_vector_to_vector(&shifted, parentedMidPoint);
    scale_vector_to_distance(&shifted, distanceFromOrigin);

    // add back in at some point
    /*    MatrixPtr zMat = MatrixPtr(new Matrix());
    zMat->rotate(0, 0, nudgeRotation.l);
    zMat->multiplyVector(&shifted);*/
    
    xyzMovement = new_vector(0, 0, nudgeTranslation.l);
    changeToXYZ->multiplyVector(&xyzMovement);
    add_vector_to_vector(&shifted, xyzMovement);
    
    mustUpdateMidPoint = false;

    if (resetNudge)
    {
        vec shift = copy_vector(shifted);
        
        take_vector_away_from_vector(parentedMidPoint, &shift);
        
        changeOfBasisMat->multiplyVector(&nudgeRotation);
        add_vector_to_vector(&rotationAngles, nudgeRotation);
        add_vector_to_vector(&arrangedMidPoint, shift);
    }
   
    if (!useParent)
    {
        vec myParentMidPoint = new_vector(0, 0, 0);
        
        if (!isLUCA())
        {
            myParentMidPoint = getParent()->midPointOffsetFromParent();
        }
        
        take_vector_away_from_vector(myParentMidPoint, &shifted);
    }
    else
    {
        quickMidPoint = shifted;
    }
    
    return shifted;
}

vec Detector::rawMidPointOffsetFromParent(bool useParent)
{
    if (isLUCA())
    {
        return arrangedMidPoint;
    }
    
    vec parentMidPoint = new_vector(0, 0, 0);
    
    if (useParent)
    {
        parentMidPoint = getParent()->midPointOffsetFromParent(useParent);
    }
    
    vec slow = getParent()->getRotatedSlowDirection();
    vec fast = getParent()->getRotatedFastDirection();
    vec cross = cross_product_for_vectors(fast, slow);
    
    vec midPointAdded = arrangedMidPoint;
    
    multiply_vector(&slow, midPointAdded.k);
    multiply_vector(&fast, midPointAdded.h);
    multiply_vector(&cross, midPointAdded.l);

    add_vector_to_vector(&parentMidPoint, slow);
    add_vector_to_vector(&parentMidPoint, fast);
    add_vector_to_vector(&parentMidPoint, cross);
    
    return parentMidPoint;
}

void Detector::removeMidPointRelativeToParent()
{
    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->removeMidPointRelativeToParent();
    }
    
    vec parentPos = getParent()->midPointOffsetFromParent(false, NULL);
    take_vector_away_from_vector(parentPos, &arrangedMidPoint);
}

void Detector::fullDescription()
{
    Detector::getMaster()->description(true);
}

void Detector::description(bool propogate)
{
    logged << "I am detector " << getTag();
    if (!isLUCA())
        logged << " and my parent is " << getParent()->getTag();
    vec rebasedNudge = new_vector(0, 0, 0);
    vec midpoint = midPointOffsetFromParent(true, &rebasedNudge);
    logged  << std::endl << "My raw midpoint is " << prettyDesc(arrangedMidPoint) << std::endl;
    logged  << std::endl << "My corrected midpoint is " << prettyDesc(midpoint) << std::endl;
    logged << "Distance of midpoint from origin is " << length_of_vector(midpoint) << " z = " << midpoint.l << std::endl;
    logged << "My midpoint relative to my parent is " << prettyDesc(midPointOffsetFromParent(false)) << std::endl;
    logged << "My rotation angles are " << std::setprecision(10) << prettyDesc(rotationAngles) << std::endl;
    logged << "My rotated fast axis is " << std::setprecision(10) << prettyDesc(getRotatedFastDirection()) << std::endl;
    logged << "My rotated slow axis is " << std::setprecision(10) << prettyDesc(getRotatedSlowDirection()) << std::endl;
    
    if (!hasChildren())
    {
        vec intersection;
        spotCoordToAbsoluteVec(unarrangedTopLeftX, unarrangedTopLeftY, &intersection);
        logged << "My top left corner maps to " << prettyDesc(intersection) << std::endl;
    }
    
    if (enabledNudge)
    {
        logged << "Rebased nudge from nudge translation is " << prettyDesc(rebasedNudge) << std::endl;
        logged << "Nudge translation is " << prettyDesc(nudgeTranslation) << " and rotation is " << prettyDesc(nudgeRotation) << std::endl;
        sendLog();
    }
    
    logged << "I have " << childrenCount() << " children." << std::endl << std::endl;
    
    sendLog();
    
    if (propogate && hasChildren())
    {
        logged << "Now for my children (" << getTag() << ")..." << std::endl;
        sendLog();
        
        for (int i = 0; i < childrenCount(); i++)
        {
            getChild(i)->description(true);
        }
        
        logged << "Finished children (" << getTag() << ")..." << std::endl;
        sendLog();
    }
}

double Detector::halfSlow()
{
    return (unarrangedBottomRightY - unarrangedTopLeftY) / 2;
}

double Detector::halfFast()
{
    return (unarrangedBottomRightX - unarrangedTopLeftX) / 2;
}

void Detector::rearrangeCoord(std::pair<float, float> *aShift)
{
    vec slow = getRotatedSlowDirection();
    vec fast = getRotatedFastDirection();
    
    double newX = aShift->first * slow.h + aShift->second * fast.h;
    double newY = aShift->first * slow.k + aShift->second * fast.k;
    
    aShift->first = newX;
    aShift->second = newY;
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

double Detector::spotCoordToResolution(double unarrangedX, double unarrangedY, double wavelength)
{
    vec arrangedPos;
    spotCoordToAbsoluteVec(unarrangedX, unarrangedY, &arrangedPos);
    scale_vector_to_distance(&arrangedPos, 1 / wavelength);
    arrangedPos.l -= 1 / wavelength;
    return length_of_vector(arrangedPos);
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
        
        vec cumulative = midPointOffsetFromParent();
        
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
        spot->getDetector()->spotCoordToAbsoluteVec(spot->getRawXY().first, spot->getRawXY().second, arrangedPos);
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
    
    vec cumulative = midPointOffsetFromParent(true, NULL);
    
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
    vec cumulative = midPointOffsetFromParent(true, NULL);
    
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
    double imageWavelength = miller->getPredictedWavelength();
    
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


std::string Detector::writeGeometryFile(int indentCount)
{
  //  lockNudges();
    std::ostringstream output;
    
    output << indents(indentCount) << "panel " << getTag() << std::endl;
    output << indents(indentCount) << "{" << std::endl;
    output << indents(indentCount + 1) << "min_fs = " << unarrangedTopLeftX << std::endl;
    output << indents(indentCount + 1) << "min_ss = " << unarrangedTopLeftY << std::endl;
    output << indents(indentCount + 1) << "max_fs = " << unarrangedBottomRightX << std::endl;
    output << indents(indentCount + 1) << "max_ss = " << unarrangedBottomRightY << std::endl;
    output << indents(indentCount + 1) << "fs_x = " << fastRotated.h << std::endl;
    output << indents(indentCount + 1) << "fs_y = " << fastRotated.k << std::endl;
    output << indents(indentCount + 1) << "fs_z = " << fastRotated.l << std::endl;
    output << indents(indentCount + 1) << "ss_x = " << slowRotated.h << std::endl;
    output << indents(indentCount + 1) << "ss_y = " << slowRotated.k << std::endl;
    output << indents(indentCount + 1) << "ss_z = " << slowRotated.l << std::endl;
    output << indents(indentCount + 1) << "midpoint_x = " << arrangedMidPoint.h << std::endl;
    output << indents(indentCount + 1) << "midpoint_y = " << arrangedMidPoint.k << std::endl;
    output << indents(indentCount + 1) << "midpoint_z = " << arrangedMidPoint.l << std::endl;
    output << indents(indentCount + 1) << "ghost = " << (hasChildren() ? "true" : "false") << std::endl;

    for (int i = 0; i < childrenCount(); i++)
    {
        output << std::endl << getChild(i)->writeGeometryFile(indentCount + 1);
    }
    
    output << indents(indentCount) << "}" << std::endl;
    
    return output.str();
}

// MARK: Resolution things.

void updateMinMax(double *min, double *max, double value)
{
    if (value < *min)
    {
        *min = value;
    }
    if (value > *max)
    {
        *max = value;
    }
}

void Detector::resolutionLimits(double *min, double *max, double wavelength)
{
    resolutionOrZLimits(min, max, wavelength, 0);
}


void Detector::zLimits(double *min, double *max)
{
    CSVPtr csv = CSVPtr(new CSV(5, "x", "y", "h", "k", "l"));

    resolutionOrZLimits(min, max, 0, 1, csv);
    
    std::string geometryList = FileParser::getKey("DETECTOR_LIST", std::string("error"));
    csv->writeToFile("image_to_detector_map_" + geometryList + ".csv");
}

// type = 0 = resolution, type = 1 = z limits
void Detector::resolutionOrZLimits(double *min, double *max, double wavelength, int type, CSVPtr csv)
{
    if (isLUCA())
    {
        *min = FLT_MAX;
        *max = -FLT_MAX;
    }
    
    if (hasChildren())
    {
        for (int i = 0; i < childrenCount(); i++)
        {
            double myChildMin = FLT_MAX;
            double myChildMax = -FLT_MAX;
            getChild(i)->resolutionOrZLimits(&myChildMin, &myChildMax, wavelength, type, csv);
            updateMinMax(min, max, myChildMin);
            updateMinMax(min, max, myChildMax);
        }
    }
    else if (type == 0)
    {
        double myMin = FLT_MAX;
        double myMax = -FLT_MAX;
        double resol = 0;
    
        resol = spotCoordToResolution(unarrangedTopLeftX, unarrangedTopLeftY, wavelength);
        updateMinMax(&myMin, &myMax, resol);
        
        resol = spotCoordToResolution(unarrangedTopLeftX, unarrangedBottomRightY, wavelength);
        updateMinMax(&myMin, &myMax, resol);

        resol = spotCoordToResolution(unarrangedBottomRightX, unarrangedTopLeftY, wavelength);
        updateMinMax(&myMin, &myMax, resol);

        resol = spotCoordToResolution(unarrangedBottomRightX, unarrangedBottomRightY, wavelength);
        updateMinMax(&myMin, &myMax, resol);

        *min = myMin;
        *max = myMax;
    }
    else if (type == 1)
    {
        double lowZ = FLT_MAX;
        double highZ = 0;
        
        vec arrangedPos = new_vector(0, 0, 0);
        spotCoordToAbsoluteVec(unarrangedTopLeftX, unarrangedTopLeftY, &arrangedPos);
        updateMinMax(&lowZ, &highZ, arrangedPos.l);
        csv->addEntry(5, (double)unarrangedTopLeftX, (double)unarrangedTopLeftY, arrangedPos.h, arrangedPos.k, arrangedPos.l);
        
        spotCoordToAbsoluteVec(unarrangedTopLeftX, unarrangedBottomRightY, &arrangedPos);
        updateMinMax(&lowZ, &highZ, arrangedPos.l);
        csv->addEntry(5, (double)unarrangedTopLeftX, (double)unarrangedBottomRightY, arrangedPos.h, arrangedPos.k, arrangedPos.l);
        
        spotCoordToAbsoluteVec(unarrangedBottomRightX, unarrangedTopLeftY, &arrangedPos);
        updateMinMax(&lowZ, &highZ, arrangedPos.l);
        csv->addEntry(5, (double)unarrangedBottomRightX, (double)unarrangedTopLeftY, arrangedPos.h, arrangedPos.k, arrangedPos.l);

        spotCoordToAbsoluteVec(unarrangedBottomRightX, unarrangedBottomRightY, &arrangedPos);
        updateMinMax(&lowZ, &highZ, arrangedPos.l);
        csv->addEntry(5, (double)unarrangedBottomRightX, (double)unarrangedBottomRightY, arrangedPos.h, arrangedPos.k, arrangedPos.l);
        
        *min = lowZ;
        *max = highZ;
    }
}

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

double Detector::millerStdevScoreWrapper(void *object)
{
    return static_cast<Detector *>(object)->millerScore(false, true);
}

double Detector::millerScoreWrapper(void *object)
{
    return static_cast<Detector *>(object)->millerScore();
}

double Detector::millerScore(bool ascii, bool stdev)
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

    std::vector<double> distances;
    
    CSVPtr csv = CSVPtr(new CSV(2, "x", "y"));
    
    double meanDistance = 0;
    double xSum = 0;
    double ySum = 0;
    
    for (int i = 0; i < xShifts.size(); i++)
    {
        double x = *xShifts[i];
        double y = *yShifts[i];
        
        xSum += x;
        ySum += y;
        
        double distance = sqrt(x * x + y * y);
        
        if (ascii)
        {
            csv->addEntry(0, x, y);
        }
        
        meanDistance += distance;
    }
    
    meanDistance /= xShifts.size();
    xSum /= xShifts.size();
    ySum /= xShifts.size();
    
    double leastSquares = 0;
    
    for (int i = 0; i < xShifts.size(); i++)
    {
        double x = *xShifts[i];
        double y = *yShifts[i];
        
        double distance = sqrt((x - xSum) * (x - xSum) + (y - ySum) * (y - ySum));
        
        if (ascii)
        {
            csv->addEntry(0, x, y);
        }
        
        leastSquares += distance;
    }
    
    leastSquares /= xShifts.size();
    
    if (ascii)
    {
        logged << std::endl << "ASCII plot of coordinates for " << getTag() << std::endl;
        sendLog();
        double edge = FileParser::getKey("METROLOGY_SEARCH_SIZE", 5);
        csv->setMinMaxXY(-edge, -edge, +edge, +edge);
        csv->plotColumns(0, 1);
    }
    
    if (!stdev)
    {
        return meanDistance;
    }
    else
    {
        return leastSquares;
    }
}

void Detector::drawSpecialImage(std::string filename)
{
    if (!drawImage)
        return;
    
    if (filename == "")
    {
        filename = "special_" + i_to_str(specialImageCounter) + ".png";
        specialImageCounter++;
    }
    
    drawImage->drawSpotsOnPNG(filename);
}

void Detector::fixMidpoints()
{
    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->fixMidpoints();
    }
    
    if (!childrenCount())
    {
        return;
    }
    
    if (isLUCA())
    {
        setUpdateMidPoint();
        return;
    }
    
    vec aveMidpoint = new_vector(0, 0, 0);
    
    for (int i = 0; i < childrenCount(); i++)
    {
        vec myMidpoint = getChild(i)->midPointOffsetFromParent(false);
        getChangeOfBasis()->multiplyVector(&myMidpoint);
        add_vector_to_vector(&aveMidpoint, myMidpoint);
    }
    
    multiply_vector(&aveMidpoint, 1 / (double)childrenCount());
    
    setArrangedMidPoint(aveMidpoint);
    
    for (int i = 0; i < childrenCount(); i++)
    {
        vec childPoint = getChild(i)->getArrangedMidPoint();
        take_vector_away_from_vector(aveMidpoint, &childPoint);
        getChild(i)->setArrangedMidPoint(childPoint);
    }
}