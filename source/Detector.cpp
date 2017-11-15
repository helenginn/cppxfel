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
#include <vector>
#include <iomanip>

double Detector::mmPerPixel = 0;
DetectorPtr Detector::masterPanel = DetectorPtr();
int Detector::detectorActive = -1;
bool Detector::noisy = false;
ImagePtr Detector::drawImage = ImagePtr();
DetectorType Detector::detectorType = DetectorTypeMPCCD;
bool Detector::enabledNudge = false;
int Detector::specialImageCounter = 0;
std::mutex Detector::setupMutex;
double Detector::cacheStep = 0;
std::vector<double> Detector::millerTargetTable;

// MARK: initialisation and constructors

void Detector::setupCache()
{
    std::lock_guard<std::mutex> lg(setupMutex);

    if (millerTargetTable.size())
    {
        return;
    }

    double maxSqr = 1;
    double intervals = 100;
    cacheStep = maxSqr / intervals; // & lower down
    intervals = 1000;
    double distSq = 0;

    for (int i = 0; i < intervals; i++)
    {
        double dist = sqrt(distSq);
        double target = exp(-fabs(dist / maxSqr));
        millerTargetTable.push_back(target);
        distSq += cacheStep;
    }
}

double Detector::lookupCache(double distSq)
{
    double category = distSq / cacheStep;

    if (category > millerTargetTable.size())
    {
                std::lock_guard<std::mutex> lg(setupMutex);

                if (millerTargetTable.size() <= category)
                {
                        double maxDistSq = cacheStep * millerTargetTable.size();

                        for (int i = (int)millerTargetTable.size(); i <= category; i++)
                        {
                                double dist = sqrt(maxDistSq);
                                double target = exp(-fabs(dist));
                                millerTargetTable.push_back(target);
                                distSq += cacheStep;
                        }
                }
    }

    return millerTargetTable[category];
}

void Detector::initialiseZeros()
{
    gain = 1;
        _changed = true;
        _cycleNum = 0;

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
    slowRotated = new_vector(0, 1, 0);
    fastRotated = new_vector(1, 0, 0);
    memset(&rotationAngles.h, 0, sizeof(rotationAngles.h) * 3);
    memset(&nudgeTranslation.h, 0, sizeof(nudgeTranslation.h) * 3);
    memset(&nudgeRotation.h, 0, sizeof(nudgeRotation.h) * 3);
    memset(&quickMidPoint.h, 0, sizeof(quickMidPoint.h) * 3);
    memset(&quickAngles.h, 0, sizeof(quickAngles.h) * 3);
    memset(&poke.h, 0, sizeof(poke.h) * 3);
    memset(&interNudge.h, 0, sizeof(interNudge.h) * 3);
    memset(&originalNudgePosition.h, 0, sizeof(originalNudgePosition.h) * 3);
        memset(&pokeLateralAxis, 0, sizeof(pokeLateralAxis.h) * 3);
        memset(&pokeLongitudinalAxis, 0, sizeof(pokeLongitudinalAxis.h) * 3);
        memset(&pokePerpendicularAxis, 0, sizeof(pokePerpendicularAxis.h) * 3);

        memset(&requiredRotToOrigin, 0, sizeof(requiredRotToOrigin.h) * 3);

    for (int i = 0; i < 4; i++)
    {
        memset(&originalCorners[i].h, 0, sizeof(originalCorners[i].h) * 3);
    }

    nudgeTiltY = 0;
    nudgeTiltX = 0;
    nudgeStep = 0;

        addPixelOffsetX = 0;
        addPixelOffsetY = 0;

        parent = DetectorPtr();
    rotMat = MatrixPtr(new Matrix());
    changeOfBasisMat = MatrixPtr(new Matrix());
    fixedBasis = MatrixPtr(new Matrix());
    invWorkingBasisMat = MatrixPtr(new Matrix());

    originalNudgeMat = MatrixPtr(new Matrix());
    interRotation = MatrixPtr(new Matrix());

    mustUpdateMidPoint = true;

    if (mmPerPixel == 0)
    {
        mmPerPixel = FileParser::getKey("MM_PER_PIXEL", 0.11);
    }

    setupCache();
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

    _refinable = true;

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

void Detector::postInit()
{
    setUpdateMidPointForDetector();

    for (int i = 0; i < 4; i++)
    {
        int xPix = (i < 2) ? unarrangedTopLeftX : unarrangedBottomRightX;
        int yPix = (i % 2 == 0) ? unarrangedTopLeftY : unarrangedBottomRightY;

        spotCoordToAbsoluteVec(xPix, yPix, &originalCorners[i]);
    }

    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->postInit();
    }
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
        std::lock_guard<std::mutex> lg(getParent()->threadMutex);
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

void Detector::prepareInterNudges()
{
    rotateAxisRecursive(true);
    resetNudgeBasis();
}

// MARK: Housekeeping and additional initialisation

void Detector::updateUnarrangedMidPoint()
{
    unarrangedMidPointX = (double)(unarrangedBottomRightX + unarrangedTopLeftX) / 2;
    unarrangedMidPointY = (double)(unarrangedBottomRightY + unarrangedTopLeftY) / 2;
}

vec Detector::getArrangedTopLeft()
{
    vec pixXY_mmZ;
    spotCoordToAbsoluteVec(unarrangedTopLeftX, unarrangedTopLeftY, &pixXY_mmZ);
    pixXY_mmZ.l *= mmPerPixel;

    return pixXY_mmZ;
}

void Detector::setArrangedTopLeft(double newX, double newY, double newZ)
{
    assert(directionSanityCheck());

    /* This is where top left offset currently is from mid point */
    vec currentArrangedBottomRightOffset = new_vector(0, 0, 0);
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

void Detector::recalculateChangeOfBasis(vec *fastAxis, vec *slowAxis)
{
    changeOfBasisMat = calculateChangeOfBasis(fastAxis, slowAxis, MatrixPtr(new Matrix()));
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

    rotMat->setIdentity();
    rotMat->rotate(angles.h, angles.k, angles.l);

    chosenMat->multiply(*rotMat);
}

void Detector::resetNudgeBasis()
{
    vec midpoint = midPointOffsetFromParent();

    vec zAxis = new_vector(0, 0, 1);
    scale_vector_to_distance(&midpoint, 1);
    double angle = angleBetweenVectors(midpoint, zAxis);
    vec cross = cross_product_for_vectors(midpoint, zAxis);
    scale_vector_to_distance(&cross, 1);

    originalNudgeMat->setIdentity();

    if (cross.h == cross.h)
    {
        originalNudgeMat->rotateRoundUnitVector(cross, angle);
    }

        smartRatio = 0;
}

void Detector::rotateAxisRecursive(bool fix)
{
    /* Will mirror change of basis if in "fix" mode */
    MatrixPtr tempNewChange = MatrixPtr(new Matrix());

    changeOfBasisMat->copyComponents(fixedBasis);

    if (fix)
    {
        tempNewChange->copyComponents(changeOfBasisMat);
    }

    if (!isLUCA())
    {
        std::lock_guard<std::mutex> lg(getParent()->threadMutex);

        MatrixPtr mat = getParent()->getChangeOfBasis();
        changeOfBasisMat->multiply(*mat);
    }

    addToBasisChange(rotationAngles, changeOfBasisMat);
    addToBasisChange(rotationAngles, tempNewChange);

    vec tempNudgeRot = nudgeRotation;
    tempNudgeRot.l = 0;
    originalNudgeMat->multiplyVector(&tempNudgeRot);
    addToBasisChange(tempNudgeRot, changeOfBasisMat);
    addToBasisChange(tempNudgeRot, tempNewChange);

    if (enabledNudge)
    {
        setUpdateMidPointForDetector();
        midPointOffsetFromParent();

        changeOfBasisMat->multiply(*interRotation);
        tempNewChange->multiply(*interRotation);
    }

    if (fix)
    {
        midPointOffsetFromParent(true, fix);
        fixedBasis->copyComponents(tempNewChange);

        memset(&rotationAngles.h, 0, sizeof(rotationAngles.h) * 3);
        memset(&nudgeTranslation.h, 0, sizeof(nudgeTranslation.h) * 3);
        memset(&nudgeRotation.h, 0, sizeof(nudgeRotation.h) * 3);
        memset(&interNudge.h, 0, sizeof(interNudge.h) * 3);
                interRotation->setIdentity();

        resetNudgeBasis();
    }

    slowRotated = slowDirection;
    changeOfBasisMat->multiplyVector(&slowRotated);
    fastRotated = fastDirection;
    changeOfBasisMat->multiplyVector(&fastRotated);

    workingBasisMat = calculateChangeOfBasis(&fastRotated, &slowRotated, invWorkingBasisMat);
    invWorkingBasisMat = workingBasisMat->inverse3DMatrix();

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
    return new_vector(workingBasisMat->components[3], workingBasisMat->components[4], workingBasisMat->components[5]);
}

vec Detector::getRotatedFastDirection()
{
    return new_vector(workingBasisMat->components[0], workingBasisMat->components[1], workingBasisMat->components[2]);
}

void Detector::enableNudge()
{
    enabledNudge = true;
    std::vector<DetectorPtr> allDetectors;
    getMaster()->getAllSubDetectors(allDetectors);

    for (int i = 0; i < allDetectors.size(); i++)
    {
        allDetectors[i]->prepareInterNudges();
    }
}

vec Detector::midPointOffsetFromParent(bool useParent, bool resetNudge)
{
    if (!mustUpdateMidPoint && useParent && !resetNudge)
    {
        return quickMidPoint;
    }

    vec parentMidPoint = new_vector(0, 0, 0);

    if (!isLUCA())
    {
        parentMidPoint = getParent()->midPointOffsetFromParent(true);
    }

    vec parentMidPointCopy = copy_vector(parentMidPoint);

    vec midPointAdded = arrangedMidPoint;

    if (!isLUCA())
    {
        getParent()->getChangeOfBasis()->multiplyVector(&midPointAdded);
    }

    add_vector_to_vector(&parentMidPoint, midPointAdded);

    interRotation->setIdentity();
    addToBasisChange(interNudge, interRotation);
    interRotation->multiplyVector(&parentMidPoint);

    double newDistance = length_of_vector(parentMidPoint) + nudgeTranslation.l;
    scale_vector_to_distance(&parentMidPoint, newDistance);

    if (!useParent || resetNudge)
    {
        vec diff = vector_between_vectors(parentMidPointCopy, parentMidPoint);

        if (!useParent)
        {
            return diff;
        }

        if (resetNudge)
        {

            if (!isLUCA())
            {
                getParent()->getChangeOfBasis()->inverse3DMatrix()->multiplyVector(&diff);
                arrangedMidPoint = diff;
            }
            else
            {
                arrangedMidPoint = parentMidPoint;
            }
        }
    }

    quickMidPoint = parentMidPoint;
    mustUpdateMidPoint = false;

    return parentMidPoint;
}

void Detector::removeMidPointRelativeToParent()
{
    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->removeMidPointRelativeToParent();
    }

    std::lock_guard<std::mutex> lg(getParent()->threadMutex);
    vec parentPos = getParent()->midPointOffsetFromParent(false, false);

    if (getParent()->isLUCA())
    {
        parentPos.l = 0;
    }

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

    vec midpoint = midPointOffsetFromParent();
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
    vec shift = new_vector(aShift->first, aShift->second, 0);

    invWorkingBasisMat->multiplyVector(&shift);

    aShift->first = shift.h;
    aShift->second = shift.k;
}

bool Detector::isAncestorOf(DetectorPtr detector)
{
    if (shared_from_this() == detector)
    {
        return true;
    }

    if (ancestorMap.count(detector))
    {
        return ancestorMap[detector];
    }

        for (int i = 0; i < childrenCount(); i++)
        {
                if (getChild(i)->isAncestorOf(detector))
                {
                        ancestorMap[detector] = true;
                        return true;
                }
        }

        std::lock_guard<std::mutex> lg(threadMutex);
    ancestorMap[detector] = false;
    return false;
}

// MARK: Spot to physical position where X-rays hit the detector

void Detector::spotCoordToRelativeVec(double unarrangedX, double unarrangedY,
                                      vec *arrangedPos)
{
    double unArrangedOffsetX = unarrangedX - unarrangedMidPointX;
    double unArrangedOffsetY = unarrangedY - unarrangedMidPointY;

    vec rearrangedOffset = new_vector(unArrangedOffsetX, unArrangedOffsetY, 0);
    invWorkingBasisMat->multiplyVector(&rearrangedOffset);

    add_vector_to_vector(arrangedPos, rearrangedOffset);
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
                                      vec *arrangedPos, ImagePtr image)
{
    *arrangedPos = midPointOffsetFromParent();

        if (image)
        {
                image->augmentMidpoint(arrangedPos);
        }

    spotCoordToRelativeVec(unarrangedX, unarrangedY, arrangedPos);
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

DetectorPtr Detector::spotToAbsoluteVec(SpotPtr spot, vec *arrangedPos, ImagePtr image)
{
        double xSpot = spot->getRawX();
        double ySpot = spot->getRawY();

        if (spot->hasDetector())
        {
                spot->getDetector()->spotCoordToAbsoluteVec(spot->getRawX(), spot->getRawY(), arrangedPos, image);
                return spot->getDetector();
        }

        DetectorPtr detector = findDetectorPanelForSpotCoord(xSpot, ySpot);

        if (detector)
        {
                detector->spotCoordToAbsoluteVec(xSpot, ySpot, arrangedPos, image);
        }

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

void Detector::intersectionToSpotCoord(vec absoluteVec, double *xSpot, double *ySpot, ImagePtr image)
{
    // make the assumption that absoluteVec is indeed in the plane of the detector

    /* Need the cumulative midpoint for this panel */

    vec cumulative = midPointOffsetFromParent();

        if (image)
        {
                image->augmentMidpoint(&cumulative);
        }

    /* Need offset from centre of the panel */
    take_vector_away_from_vector(cumulative, &absoluteVec);

    /* how many slow directions vs how many fast directions...? */
    /* i.e. change of basis... */

    workingBasisMat->multiplyVector(&absoluteVec);

    /* Now to add back the unarranged mid point */
    absoluteVec.h += unarrangedMidPointX;
    absoluteVec.k += unarrangedMidPointY;

    *xSpot = absoluteVec.h;
    *ySpot = absoluteVec.k;
}

void Detector::intersectionWithRay(vec ray, vec *intersection, ImagePtr image)
{
    vec cumulative = midPointOffsetFromParent();

        if (image)
        {
                image->augmentMidpoint(&cumulative);
        }

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
    vec hkl = miller->getRay();
        ImagePtr image = miller->getImage();

    DetectorPtr probe;
    bool hasDetector = miller->hasDetector();

    if (!hasDetector)
    {
        probe = detectorForRayIntersection(hkl, intersection);

        if (probe)
        {
            miller->setDetector(probe);
        }
    }
    else
    {
        probe = miller->getDetector();
        probe->intersectionWithRay(hkl, intersection, image);
    }

    return probe;
}

DetectorPtr Detector::spotCoordForMiller(MillerPtr miller, double *xSpot, double *ySpot)
{
        ImagePtr image = miller->getImage();
        vec intersection;
    DetectorPtr probe = intersectionForMiller(miller, &intersection);

    if (probe)
    {
        probe->intersectionToSpotCoord(intersection, xSpot, ySpot, image);
    }

    return probe;
}

// MARK: Keeping track of added Millers

void Detector::addMillerCarefully(MillerPtr miller)
{
    millerMutex.lock();

        if (miller->reachesThreshold())
        {
                miller->centrePeak();
                millers.push_back(miller);
        }

        allMillers.push_back(miller);

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

vec Detector::getDifferenceFromOriginal(int corner)
{
        int xPix = (corner < 2) ? unarrangedTopLeftX : unarrangedBottomRightX;
        int yPix = (corner % 2 == 0) ? unarrangedTopLeftY : unarrangedBottomRightY;

        vec diff;
        spotCoordToAbsoluteVec(xPix, yPix, &diff);

        take_vector_away_from_vector(originalCorners[corner], &diff);

        return diff;
}

std::string Detector::writeGeometryFile(int fileCount, int indentCount, CSVPtr differenceCSV)
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
        output << indents(indentCount + 1) << "refine = " << (_refinable ? "true" : "false") << std::endl;
        output << indents(indentCount + 1) << "gain = " << gain << std::endl;

    if (!differenceCSV)
    {
        differenceCSV = CSVPtr(new CSV(3, "xDiff", "yDiff", "zDiff"));
    }

    for (int i = 0; i < 4; i++)
    {
                if (hasChildren())
                {
                        break;
                }

                vec diff = getDifferenceFromOriginal(i);
                differenceCSV->addEntry(3, diff.h, diff.k, diff.l);
    }

    for (int i = 0; i < childrenCount(); i++)
    {
        output << std::endl << getChild(i)->writeGeometryFile(fileCount, indentCount + 1, differenceCSV);
    }

    output << indents(indentCount) << "}" << std::endl;

    if (indentCount == 0)
    {
        differenceCSV->writeToFile("cgeom_diff_" + i_to_str(fileCount) + ".csv");
    }

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
        staticLogAndExit(logged);
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

void Detector::calculateMillerShifts(bool stdev)
{
        Miller::refreshMillerPositions(millers);

        xShifts.clear();
        yShifts.clear();

        for (int i = 0; i < millerCount(); i++)
        {
                MillerPtr myMiller = miller(i);

                if (!myMiller || !myMiller->reachesThreshold())
                        continue;

                float *shiftX = NULL;
                float *shiftY = NULL;

                if (stdev)
                {
                        // pixels
                        shiftX = myMiller->getXShiftPointer();
                        shiftY = myMiller->getYShiftPointer();
                }
                else if (!stdev)
                {
                        shiftX = myMiller->getRecipXShiftPtr();
                        shiftY = myMiller->getRecipYShiftPtr();
                }

                xShifts.push_back(shiftX);
                yShifts.push_back(shiftY);
        }
}

double Detector::millerStdevScoreWrapper(void *object)
{
    return static_cast<Detector *>(object)->millerScore(false, true);
}

double Detector::millerScoreWrapper(void *object)
{
    return static_cast<Detector *>(object)->millerScore();
}

double Detector::millerScore(bool ascii, bool stdev, int number)
{
        if (!millers.size())
        {
                return 0;
        }

        double geometryRlpSize = FileParser::getKey("INDEXING_RLP_SIZE", 0.0020);

        calculateMillerShifts(stdev);

    std::vector<double> distances;

    CSVPtr csv = CSVPtr(new CSV(2, "x", "y"));

    double totalScore = 0;
        double biggestContribution = 0;
        Coord bestPix;

    int count = 0;
    double maxSqr = 9;

    if (!stdev)
    {
        maxSqr = 5.0 * geometryRlpSize;
        maxSqr *= maxSqr;
    }

    for (int i = 0; i < xShifts.size(); i++)
    {
        double x = *xShifts[i];
        double y = *yShifts[i];

        if (x != x || y != y)
        {
            continue;
        }

                double contribution = 0;
                double x0 = 0;
        double y0 = 0;

        if (stdev)
        {
            x0 = x;
            y0 = y;
        }

        for (int j = 0; j < xShifts.size(); j++)
        {
            count++;

            double x1 = *xShifts[j];
            double y1 = *yShifts[j];

            double xDiff = x1 - x0;
            double yDiff = y1 - y0;

            double distSq = xDiff * xDiff + yDiff * yDiff;
            distSq /= maxSqr;

            if ((ascii || number >= 0) && i == 0)
            {
                csv->addEntry(0, x1, y1);
            }

            contribution += lookupCache(distSq);

        }

                totalScore -= contribution;

                if (stdev && contribution > biggestContribution)
                {
                        biggestContribution = contribution;
                        bestPix = std::make_pair(x0, y0);
                }

                if (!stdev || (stdev && number >= 0))
        {
            break;
        }
    }

        // assuming all Millers have same wavelength
        double wavelength = miller(0)->getImage()->getWavelength();
        double radius = 1 / wavelength;

        // to get the same thing in radians.
        requiredRotToOrigin = new_vector(bestPix.first * radius, bestPix.second * radius, 0);
        //workingBasisMat->multiplyVector(&requiredRotToOrigin);

        if (count)
    {
        totalScore /= count;
    }

    if (ascii)
    {
                bool approximate = FileParser::getKey("GEOMETRY_IS_APPROXIMATE", false);

        double edge = FileParser::getKey("METROLOGY_SEARCH_SIZE", 3) + 0.5;

        if (!stdev)
        {
            edge = 5 * geometryRlpSize;
                }

                if (approximate)
                {
                        edge *= 3;
                }

        if (number < 0)
        {
            logged << std::endl << "ASCII plot of coordinates for " << getTag() << std::endl;
            sendLog();
            csv->setMinMaxXY(-edge, -edge, +edge, +edge);
            csv->plotColumns(0, 1);
        }

        if (number >= 0)
        {
            std::map<std::string, std::string> plotMap;
            plotMap["filename"] = getTag() + "_miller_" + (stdev ? "pix_" : "recip_") + i_to_str(number);
            plotMap["height"] = "1000";
            plotMap["width"] = "900";
            plotMap["xHeader0"] = "x";
            plotMap["yHeader0"] = "y";

            plotMap["xMax0"] = f_to_str(edge);
            plotMap["xMin0"] = f_to_str(-edge);
            plotMap["yMax0"] = f_to_str(edge);
            plotMap["yMin0"] = f_to_str(-edge);

            if (stdev)
            {
                plotMap["xTitle0"] = "x shift in pix";
                plotMap["yTitle0"] = "y shift in pix";
            }
            else
            {
                plotMap["xTitle0"] = "x shift in inv Ang";
                plotMap["yTitle0"] = "y shift in inv Ang";
            }
            plotMap["style0"] = "scatter";
            plotMap["colour0"] = "blue";
                        csv->writeToFile("csv_" + plotMap["filename"] + ".csv");
            csv->plotPNG(plotMap);
        }
    }

    return totalScore;
}

void Detector::quickJumpToOrigin()
{
        millerScore();
        logged << getTag() << " - " << requiredRotToOrigin.h << ", " << requiredRotToOrigin.k << std::endl;
        sendLog();
        resetNudgeBasis();
        setInterNudgeX(this, requiredRotToOrigin.h);
        setInterNudgeY(this, requiredRotToOrigin.k);
        resetNudgeBasis();
}

void Detector::transferSingleMillers()
{
        for (int i = 0; i < allMillers.size(); i++)
        {
                MillerPtr myMiller = allMillers[i].lock();

                if (!myMiller || !myMiller->reachesThreshold())
                        continue;

                addMillerCarefully(myMiller);
        }
}

double Detector::transferMillers(void *object)
{
        Detector *me = static_cast<Detector *>(object);

        Miller::refreshMillerPositions(me->allMillers, true);
        Detector::getMaster()->clearMillers();

        std::vector<DetectorPtr> allChildren;
        getMaster()->getAllSubDetectors(allChildren, true);

        for (int i = 0; i < allChildren.size(); i++)
        {
                allChildren[i]->transferSingleMillers();
        }

        return 0;
}

double Detector::peakScore()
{
        if (!millers.size())
        {
                return 0;
        }

        double peakScore = 0;
        double total = 0;

        Miller::refreshMillerPositions(allMillers, true);

        for (int i = 0; i < allMillers.size(); i++)
        {
                MillerPtr myMiller = allMillers[i].lock();

                if (!myMiller)
                        continue;

                if (myMiller->getRawestIntensity() == myMiller->getRawestIntensity())
                {
                        total++;
                }

                if (myMiller->reachesThreshold())
                {
                        peakScore++;
                }
        }

        return -peakScore;
}

double Detector::peakScoreWrapper(void *object)
{
        return static_cast<Detector *>(object)->peakScore();
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

    drawImage->drawSpotsOnPNG(filename, true);
}

// MARK: others

vec Detector::getAverageMidpoint()
{
        vec aveMidpoint = new_vector(0, 0, 0);

        for (int i = 0; i < childrenCount(); i++)
        {
                vec myMidpoint = getChild(i)->midPointOffsetFromParent(false);
                getChangeOfBasis()->multiplyVector(&myMidpoint);
                add_vector_to_vector(&aveMidpoint, myMidpoint);
        }

        multiply_vector(&aveMidpoint, 1 / (double)childrenCount());

        return aveMidpoint;
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
        setUpdateMidPointForDetector();
        return;
    }

        vec aveMidpoint = getAverageMidpoint();

    setArrangedMidPoint(aveMidpoint);

    for (int i = 0; i < childrenCount(); i++)
    {
        vec childPoint = getChild(i)->getArrangedMidPoint();
        take_vector_away_from_vector(aveMidpoint, &childPoint);
        getChild(i)->setArrangedMidPoint(childPoint);
    }
}

void Detector::reportMillerScores(int refinementNum)
{
    if (isLUCA())
    {
        logged << "Writing pixel and reciprocal offset graph..." << std::endl;
        sendLog();
 }
        millerScore(true, false, refinementNum);
        millerScore(true, true, refinementNum);

    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->reportMillerScores(refinementNum);
    }
}

std::vector<DetectorPtr> Detector::getSubDetectorsOnLevel(int level)
{
        std::vector<DetectorPtr> lastDetectors, currentDetectors;
        lastDetectors.push_back(getMaster());
        int currentLevel = 0;

        while (currentLevel < level)
        {
                currentDetectors.clear();

                for (int i = 0; i < lastDetectors.size(); i++)
                {
                        for (int j = 0; j < lastDetectors[i]->childrenCount(); j++)
                        {
                                currentDetectors.push_back(lastDetectors[i]->getChild(j));
                        }
                }

                lastDetectors = currentDetectors;
                currentLevel++;
        }

        return lastDetectors;
}

void Detector::getAllSubDetectors(std::vector<DetectorPtr> &array, bool childrenOnly)
{
    for (int i = 0; i < childrenCount(); i++)
    {
        getChild(i)->getAllSubDetectors(array, childrenOnly);
    }

    if (childrenOnly)
    {
                if (!childrenCount() && !isSetRefinable())
                {
                        return;
                }

                bool anyRefinable = false;

                for (int i = 0; i < childrenCount(); i++)
                {
                        if (getChild(i)->isSetRefinable())
                        {
                                anyRefinable = true;
                        }
                }

                if ((!isSetRefinable() && !anyRefinable) ||
                        (isSetRefinable() && anyRefinable))
                {
                        return;
                }
    }

    array.push_back(shared_from_this());

}

void Detector::resetPoke()
{
    if (childrenCount() < 2)
    {
        return;
    }

    prepareInterNudges();

    poke.h = 0;
    poke.k = 0;
    poke.l = 0;

    vec firstChildMidpoint = getChild(0)->midPointOffsetFromParent();
    vec secondChildMidpoint = getChild(1)->midPointOffsetFromParent();
    vec travelToSecond = vector_between_vectors(firstChildMidpoint, secondChildMidpoint);
    multiply_vector(&travelToSecond, 0.5);

    vec betweenMidpoints = firstChildMidpoint;
    add_vector_to_vector(&betweenMidpoints, travelToSecond);

    pokePerpendicularAxis = betweenMidpoints;
    scale_vector_to_distance(&pokePerpendicularAxis, 1);
    vec horiz = new_vector(1, 0, -pokePerpendicularAxis.h / pokePerpendicularAxis.l);
    scale_vector_to_distance(&horiz, 1);
    vec vert = cross_product_for_vectors(pokePerpendicularAxis, horiz);

    double matComponents[9] = {horiz.h, vert.h, cross.h,
        horiz.k, vert.k, cross.k,
        horiz.l, vert.l, cross.l};

    MatrixPtr changeToXYZ = MatrixPtr(new Matrix(matComponents));

    scale_vector_to_distance(&betweenMidpoints, 1);
    changeToXYZ->multiplyVector(&travelToSecond);
    travelToSecond.l = 0;
    changeToXYZ->inverse3DMatrix()->multiplyVector(&travelToSecond);
    scale_vector_to_distance(&travelToSecond, 1);

    pokeLongitudinalAxis = travelToSecond;
    pokeLateralAxis = cross_product_for_vectors(pokePerpendicularAxis, pokeLongitudinalAxis);
}

void Detector::setPokeN(int pokeNum, double pokeValue)
{
    if (pokeNum == 0)
    {
        poke.h = pokeValue;
    }
    else if (pokeNum == 1)
    {
        poke.k = pokeValue;
    }
    else if (pokeNum == 2)
    {
        poke.l = pokeValue;
    }

    for (int i = 0; i < 2; i++)
    {
        int modifier = (i == 0) ? -1 : 1;

        vec latPoke = copy_vector(pokeLateralAxis);
        multiply_vector(&latPoke, poke.k);
        vec longPoke = copy_vector(pokeLongitudinalAxis);
        multiply_vector(&longPoke, poke.h);
        vec crossPoke = copy_vector(pokePerpendicularAxis);
        multiply_vector(&crossPoke, poke.l);
        vec allPokes = copy_vector(latPoke);
        add_vector_to_vector(&allPokes, longPoke);
        add_vector_to_vector(&allPokes, crossPoke);

        setInterNudgeX(&*getChild(i), modifier * allPokes.h);
        setInterNudgeY(&*getChild(i), modifier * allPokes.k);
        setInterNudgeZ(&*getChild(i), modifier * allPokes.l);
    }
}

double Detector::distanceFromSample()
{
    vec midpoint = midPointOffsetFromParent();
    return length_of_vector(midpoint);
}

void Detector::nudgeTiltAndStep(double *nudgeTiltX, double *nudgeTiltY, double *nudgeStep, double *interNudge)
{
    double expectedPixels = FileParser::getKey("EXPECTED_GEOMETRY_MOVEMENT", 0.02);
        *nudgeStep = expectedPixels;

    vec midpoint = midPointOffsetFromParent();
    double distance = length_of_vector(midpoint);
    MatrixPtr nudgeTranspose = originalNudgeMat->transpose();

    vec xWorkingAxis = new_vector(1, 0, 0);
    vec yWorkingAxis = new_vector(0, 1, 0);
    vec xAxisNudge = new_vector(1, 0, 0);
    vec yAxisNudge = new_vector(0, 1, 0);
    vec zAxisNudge = new_vector(0, 0, 1);
    vec zAxisReal = new_vector(0, 0, 1);

    // how many basis X and basis Y vectors do we need

    originalNudgeMat->multiplyVector(&xAxisNudge);
    originalNudgeMat->multiplyVector(&yAxisNudge);
    originalNudgeMat->multiplyVector(&zAxisNudge);

    changeOfBasisMat->multiplyVector(&zAxisReal);

    MatrixPtr zRotate = rotation_between_vectors(zAxisNudge, zAxisReal);
    zRotate->multiplyVector(&xWorkingAxis);
    zRotate->multiplyVector(&yWorkingAxis);

    // heading towards zero means we must send the detector BACK, so tilt is opposite
    // sign of step. The step is proportional to both tan(angle). For shorter distances,
    // the nudged step is smaller than for large distances.

    double xAngle = angleBetweenVectors(xWorkingAxis, xAxisNudge);
    *nudgeTiltX = expectedPixels * tan(xAngle) / distance * (xAngle < 0 ? -1 : 1);

    if (*nudgeTiltX != *nudgeTiltX)
    {
        *nudgeTiltX = 0;
    }

    double yAngle = angleBetweenVectors(yWorkingAxis, yAxisNudge);
    *nudgeTiltY = expectedPixels * tan(yAngle) / distance * (xAngle < 0 ? -1 : 1);

    if (*nudgeTiltY != *nudgeTiltY)
    {
        *nudgeTiltY = 0;
    }

    if (interNudge)
    {
        *interNudge = atan(expectedPixels / distance);
    }

    this->nudgeTiltX = *nudgeTiltX;
    this->nudgeTiltY = *nudgeTiltY;
    this->nudgeStep = *nudgeStep;
}

int Detector::spotCountFromImages(std::vector<ImagePtr> images, bool _delete)
{
        int total = 0;

        for (int i = 0; i < images.size(); i++)
        {
                bool modifiedImage = false;
                int imageTotal = 0;

                for (int j = 0; j < images[i]->spotCount(); j++)
                {
                        SpotPtr spot = images[i]->spot(j);
                        spot->estimatedVector();

                        if (spot->hasDetector())
                        {
                                if (isAncestorOf(spot->getDetector()))
                                {
                                        imageTotal++;

                                        if (_delete)
                                        {
                                                modifiedImage = true;
                                                images[i]->removeSpot(j);
                                                j--;
                                        }
                                }
                        }
                }
                // no point having one spot on an image
                if (imageTotal > 2) total += imageTotal;

                if (modifiedImage)
                {
                        images[i]->compileDistancesFromSpots();
                }
        }

        return total;
}

void Detector::flattenDetector()
{
        for (int i = 0; i < childrenCount(); i++)
        {
                getChild(i)->flattenDetector();
        }

        if (!isLUCA())
        {
                arrangedMidPoint.l = 0;
        }

        closest_major_axis(slowDirection);
        closest_major_axis(fastDirection);

        changeOfBasisMat->setIdentity();
        fixedBasis->setIdentity();

        if (hasChildren())
        {
                workingBasisMat->setIdentity();
                invWorkingBasisMat->setIdentity();
        }

        if (isLUCA())
        {
                rotateAxisRecursive(true);
        }

        mustUpdateMidPoint = true;
}

std::string toCrystFELAxis(vec axis)
{
        bool posH = (axis.h >= 0);
        bool posK = (axis.k >= 0);
        bool posL = (axis.l >= 0);

        std::ostringstream aString;
        aString << std::fixed << std::setprecision(10);
        aString << (posH ? "+" : "") << axis.h << "x " << (posK ? "+" : "") << axis.k << "y " << (posL ? "+" : "") << axis.l << "z" << std::endl;

        return aString.str();
}

void Detector::addCrystFELInfo(std::ostringstream &geometry, double clen)
{
        double pixSizeInMetres = 1000 / mmPerPixel;

        vec intersection;
        spotCoordToAbsoluteVec(unarrangedTopLeftX, unarrangedTopLeftY, &intersection);
        double coffset = intersection.l / pixSizeInMetres - clen;

        geometry << std::fixed << std::setprecision(10);
        geometry << getTag() << "/min_fs = " << unarrangedTopLeftX << std::endl;
        geometry << getTag() << "/min_ss = " << unarrangedTopLeftY << std::endl;
        geometry << getTag() << "/max_fs = " << unarrangedBottomRightX << std::endl;
        geometry << getTag() << "/max_ss = " << unarrangedBottomRightY << std::endl;
        geometry << getTag() << "/res = " << pixSizeInMetres << std::endl;
        geometry << getTag() << "/coffset = " << coffset << std::endl;

        vec slow = new_vector(workingBasisMat->components[1], workingBasisMat->components[5], workingBasisMat->components[9]);
        vec fast = new_vector(workingBasisMat->components[0], workingBasisMat->components[4], workingBasisMat->components[8]);

        geometry << getTag() << "/fs = " << toCrystFELAxis(fast);
        geometry << getTag() << "/ss = " << toCrystFELAxis(slow);

        geometry << getTag() << "/corner_x = " << intersection.h << std::endl;
        geometry << getTag() << "/corner_y = " << intersection.k << std::endl;
        geometry << std::endl;
}

std::string Detector::writeCrystFELFile()
{
    std::ostringstream geometry;

        geometry << "; CrystFEL geometry file spit out of cppxfel." << std::endl;

    double pixSizeInMetres = 1000 / mmPerPixel;
        vec lucaMidpoint = midPointOffsetFromParent();

        double detDist = lucaMidpoint.l / pixSizeInMetres;

        geometry << "res = " << pixSizeInMetres << std::endl;
        geometry << "clen = " << detDist << std::endl << std::endl;

        geometry << "; You may insert some things here probably." << std::endl << std::endl;

        std::vector<DetectorPtr> allDets;
        getAllSubDetectors(allDets);

        for (int i = 0; i < allDets.size(); i++)
        {
                if (allDets[i]->hasNoChildren())
                {
                        allDets[i]->addCrystFELInfo(geometry, detDist);
                }
        }

    return geometry.str();
}
