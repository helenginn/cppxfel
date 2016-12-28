//
//  Detector.h
//  cppxfel
//
//  Created by Helen Ginn on 26/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__Detector__
#define __cppxfel__Detector__

#include <stdio.h>
#include "parameters.h"
#include "LoggableObject.h"
#include "FileParser.h"

class Detector : public LoggableObject, public boost::enable_shared_from_this<Detector>
{
private:
    /* MARK: static simple types */
    static int detectorActive;
    static double mmPerPixel;
    static DetectorPtr masterPanel;
    
    /* MARK: Simple type class members */
    
    /* These map onto coordinates from Cheetah HDF5 */
    int unarrangedTopLeftX;
    int unarrangedTopLeftY;
    int unarrangedBottomRightX;
    int unarrangedBottomRightY;
    
    /* Easy-access: this unarranged pixel is considered the
     * centre of the panel. This is derived from the Cheetah
     * HDF5 coordinates */
    double unarrangedMidPointX;
    double unarrangedMidPointY;
    
    /* This is considered the top left of the panel in arranged space.*/
    /* This could be lifted straight out of CrystFEL */
    vec arrangedTopLeft;
    
    /* This is considered the centre of the panel in arranged space.
     * This should be derived from the arranged top left X/Y */
    vec arrangedMidPoint;
    
    /* Basis vectors can be lifted out of CrystFEL but these will be 
     * defined with regards to the top left corner. If there is a Z
     * component then this will mess up the mid point positioning.
     * Fix this later! */
    
    /* This is the basis vectors corresponding to mapped slow-speed axis */
    /* Think: s -> y or vertical */
    vec slowDirection;
    
    /* This is the basis vectors corresponding to mapped fast-speed axis */
    /* Think: f -> x or horizontal */
    vec fastDirection;
    
    /* This is stored for convenience */
    vec cross;
    
    /* These are the Tait-Bryan angles which are additionally applied */
    /* to the slow-speed and fast-speed axes. Radians. */
    vec rotationAngles;
    
    // MARK: Class type class members */
    
    std::string tag;
    std::vector<DetectorPtr> children;
    DetectorWeakPtr parent;
    MatrixPtr rotMat;
    MatrixPtr changeOfBasisMat;
    std::vector<MillerPtr> millers;
    std::mutex millerMutex;
    
    // MARK: private housekeeping
    
    void initialiseZeros();
    void updateUnarrangedMidPoint();
    bool directionSanityCheck();
    void updateCurrentRotation();
    
    // MARK: private calculations
    
    void spotCoordToRelativeVec(double unarrangedX, double unarrangedY,
                                vec *arrangedPos);
    vec midPointOffsetFromParent();
    void cumulativeMidPoint(vec *cumulative);
    
    

public:
    /* Initialise all variables to zero */
    Detector();
    
    /* For master panel: initialising using detector distance, beamX, beamY */
    /* Then set Detector::masterPanel to this pointer */
    Detector(double distance, double beamX, double beamY);
    
    /* For children panels, using initialisation variables from CrystFEL */
    /* Follow with parent->addChild to complete the process */
    Detector(DetectorPtr parent, Coord unarrangedTopLeft, Coord unarrangedBottomRight,
             vec slowDir, vec fastDir, vec arrangedTopLeft, bool lastIsMiddle = false);
    
    /* Are we using detectors? */
    
    static bool isActive()
    {
        if (detectorActive == -1)
        {
            detectorActive = FileParser::getKey("USE_NEW_DETECTOR_FORMAT", false);
        }
        
        return detectorActive;
    }
    
    // MARK: public initialisation getters/setters
    
    void setUnarrangedTopLeft(int newX, int newY)
    {
        unarrangedTopLeftX = newX;
        unarrangedTopLeftY = newY;
        
        updateUnarrangedMidPoint();
    }
    
    void setUnarrangedBottomRight(int newX, int newY)
    {
        unarrangedBottomRightX = newX;
        unarrangedBottomRightY = newY;
        
        updateUnarrangedMidPoint();
    }
    
    void setSlowDirection(double newX, double newY, double newZ)
    {
        slowDirection.h = newX;
        slowDirection.k = newY;
        slowDirection.l = newZ;
        
        updateCurrentRotation();
    }
    
    void setFastDirection(double newX, double newY, double newZ)
    {
        fastDirection.h = newX;
        fastDirection.k = newY;
        fastDirection.l = newZ;

        updateCurrentRotation();
    }
    
    void setArrangedTopLeft(double newX, double newY, double newZ);
    
    vec getArrangedMidPoint()
    {
        return arrangedMidPoint;
    }
    
    void prepareRotationAngles(double alpha, double beta, double gamma)
    {
        rotationAngles.h = alpha;
        rotationAngles.k = beta;
        rotationAngles.l = gamma;
    }
    
    /* Supply angles in radians */
    void setRotationAngles(double alpha, double beta, double gamma)
    {
        rotationAngles.h = alpha;
        rotationAngles.k = beta;
        rotationAngles.l = gamma;
        
        updateCurrentRotation();
    }
    
    void setTag(std::string newTag)
    {
        tag = newTag;
    }
    
    std::string getTag()
    {
        return tag;
    }
    
    // MARK: Higher level getters and setters
    
    void rotateDirection(vec *direction);
    vec getRotatedSlowDirection();
    vec getRotatedFastDirection();
    
    double halfSlow();
    double halfFast();
    
    MatrixPtr getChangeOfBasis()
    {
        return changeOfBasisMat;
    }
    
    MatrixPtr getRotation()
    {
        return rotMat;
    }
    
    // MARK: Public family functions
    
    void setParent(DetectorPtr myParent)
    {
        parent = myParent;
    }
    
    DetectorPtr getParent()
    {
        return parent.lock();
    }
    
    bool isLUCA()
    {
        return (!parent.lock());
    }
    
    size_t childrenCount()
    {
        return children.size();
    }
    
    bool hasNoChildren()
    {
        return (childrenCount() == 0);
    }
    
    bool hasChildren()
    {
        return (childrenCount() > 0);
    }
    
    bool isAncestorOf(DetectorPtr detector);
    
    DetectorPtr getChild(int i)
    {
        return children[i];
    }
    
    void addChild(DetectorPtr newD)
    {
        newD->setParent(shared_from_this());
        children.push_back(newD);
    }
    
    static DetectorPtr getMaster()
    {
        return masterPanel;
    }
    
    static void setMaster(DetectorPtr newMaster)
    {
        masterPanel = newMaster;
        newMaster->setTag("master");
    }

    // MARK: apply rotations bottom up
    
    void applyRotations();
    
    // MARK: Spot coord to absolute vec
    
    DetectorPtr spotToAbsoluteVec(SpotPtr spot, vec *arrangedPos);
    
    /* If you know the detector panel,
     * you can call this to get the arranged position */
    void spotCoordToAbsoluteVec(double unarrangedX, double unarrangedY,
                                vec *arrangedPos);
    
    /* If you don't know the detector panel, find it using this function */
    DetectorPtr findDetectorPanelForSpotCoord(double xSpot, double ySpot);
    
    /* Both at once - ask Master Panel */
    DetectorPtr findDetectorAndSpotCoordToAbsoluteVec(double unarrangedX, double unarrangedY,
                                                      vec *arrangedPos);
    
    /* Is spot within bounds of unarranged panels */
    bool spotCoordWithinBounds(double xSpot, double ySpot);
    
    // MARK: Absolute vec to spot coord
    
    /* If you know the detector panel, arranged position to unarranged */
    void intersectionToSpotCoord(vec absoluteVec, double *xSpot, double *ySpot);
    
    // MARK: ray trace to absolute vec
    
    DetectorPtr detectorForRayIntersection(vec ray, vec *intersection);
    void intersectionWithRay(vec ray, vec *intersection);
    DetectorPtr spotCoordForRayIntersection(vec ray, double *xSpot, double *ySpot);
    
    /* Keeping track of all Millers */
    
    void addMillerCarefully(MillerPtr miller);
    
    size_t millerCount()
    {
        return millers.size();
    }
    
    MillerPtr miller(int i)
    {
        return millers[i];
    }
    
    /* Write geometry file */
    
    std::string writeGeometryFile(int indentCount = 0);
    
    /* Refinement getter/setters */
    
    static void setArrangedMidPointX(void *object, double newX)
    {
        static_cast<Detector *>(object)->arrangedMidPoint.h = newX;
    }
    
    static double getArrangedMidPointX(void *object)
    {
        return static_cast<Detector *>(object)->arrangedMidPoint.h;
    }
    
    static void setArrangedMidPointY(void *object, double newY)
    {
        static_cast<Detector *>(object)->arrangedMidPoint.k = newY;
    }
    
    static double getArrangedMidPointY(void *object)
    {
        return static_cast<Detector *>(object)->arrangedMidPoint.k;
    }
    
    static void setArrangedMidPointZ(void *object, double newZ)
    {
        static_cast<Detector *>(object)->arrangedMidPoint.l = newZ;
    }
    
    static double getArrangedMidPointZ(void *object)
    {
        return static_cast<Detector *>(object)->arrangedMidPoint.l;
    }
    
    static void setAlpha(void *object, double newAlpha)
    {
        static_cast<Detector *>(object)->rotationAngles.h = newAlpha;
        static_cast<Detector *>(object)->updateCurrentRotation();
    }
    
    static double getAlpha(void *object)
    {
        return static_cast<Detector *>(object)->rotationAngles.h;
    }
    
    static void setBeta(void *object, double newBeta)
    {
        static_cast<Detector *>(object)->rotationAngles.k = newBeta;
        static_cast<Detector *>(object)->updateCurrentRotation();
    }
    
    static double getBeta(void *object)
    {
        return static_cast<Detector *>(object)->rotationAngles.k;
    }
    
    static void setGamma(void *object, double newGamma)
    {
        static_cast<Detector *>(object)->rotationAngles.l = newGamma;
        static_cast<Detector *>(object)->updateCurrentRotation();
    }
    
    static double getGamma(void *object)
    {
        return static_cast<Detector *>(object)->rotationAngles.l;
    }

};

#endif /* defined(__cppxfel__Detector__) */
