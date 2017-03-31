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
#include <mutex>

typedef enum
{
    DetectorTypeCSPAD,
    DetectorTypeMPCCD,
} DetectorType;

typedef std::map<DetectorPtr, bool> AncestorMap;

class Detector : public LoggableObject, public boost::enable_shared_from_this<Detector>
{
private:
    /* MARK: static simple types */
    IndexManagerPtr _manager;
    static int detectorActive;
    static double mmPerPixel;
    static DetectorPtr masterPanel;
    static bool noisy;
    static ImagePtr drawImage;
    static DetectorType detectorType;
    static int specialImageCounter;
    double gain;
    bool _refinable;
    std::mutex threadMutex;
    
    static bool enabledNudge;
    bool mustUpdateMidPoint;
    
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
    
    /* This is the latest mid point/angles NOT relative to the parent for
       retrieval if no updates have been made. */
    vec quickMidPoint;
    vec quickAngles;
    
    /* Basis vectors can be lifted out of CrystFEL but these will be 
     * defined with regards to the top left corner. If there is a Z
     * component then this will mess up the mid point positioning.
     * Fix this later! */
    
    /* Horizontal/vertical/closeness shift */
    vec nudgeTranslation;
    
    /* Nudge translation (x, y) and nudge rotation (z) which would be
       applied (equal but opposite) to the expected two child panels.
       In the spirit of ridiculous terminology, I am calling this a poke.
     */
    
    vec poke;
    
    /* The basis vectors of the poke - so now x will move panels laterally
       w.r.t. each other and y will move panels towards or away from each other.
     */
    
    vec pokeLateralAxis;
    vec pokeLongitudinalAxis;
    vec pokePerpendicularAxis;
    
    /* Internudge rotations */
    
    vec interNudge;
    
    /* Nudge (angles as per the intraPanel movement) */
    vec nudgeRotation;

    /* This is the basis vectors corresponding to mapped slow-speed axis */
    /* Think: s -> y or vertical */
    vec slowDirection;
    vec slowRotated;
    
    /* This is the basis vectors corresponding to mapped fast-speed axis */
    /* Think: f -> x or horizontal */
    vec fastDirection;
    vec fastRotated;
    
    /* This is stored for convenience */
    vec cross;
    
    /* These are the Tait-Bryan angles which are additionally applied */
    /* to the slow-speed and fast-speed axes. Radians. */
    vec rotationAngles;
    
    // MARK: Class type class members */
    
    std::string tag;
    std::vector<DetectorPtr> children;
    DetectorWeakPtr parent;
    
    /* Just for re-use it seems */
    MatrixPtr rotMat;
    
    /* To go from 'real' axes into 'panel' axes (but NOT map to pixels) */
    MatrixPtr changeOfBasisMat;
    
    /* To go from 'real' axes into 'panel' axes (AND map to pixels) */
    MatrixPtr workingBasisMat;
    /* From relative pixel movements into 'real' axes */
    MatrixPtr invWorkingBasisMat;
    
    /* Storing the bit which isn't re-calculated */
    MatrixPtr fixedBasis;
    
    /* From 'real' axes into 'eye placed at sample looking at midpoint
       of detector' axes */
    MatrixPtr originalNudgeMat;
     
    /* Matrix representing rotation around sample */
    MatrixPtr interRotation;
    
    /* midpoint of detector prior to any session of nudging */
    vec originalNudgePosition;
    
    /* Keep track of nudge step ratios */
    
    double nudgeStep;
    double nudgeTiltX;
    double nudgeTiltY;
    
    vec originalCorners[4];
    
    /* Ancestor map for quick access for previously asked "are you my ancestor"? */
    AncestorMap ancestorMap;
    
    // MARK: keeping track of millers, etc.
    std::vector<MillerWeakPtr> millers;
    std::mutex millerMutex;
    static double cacheStep;

    /* In the interests of caching distance target function */
    void setupCache();
    static std::mutex setupMutex;
    static std::vector<double> millerTargetTable;
    
    /* These keep track of the pointers to the Miller index shifts for easy access */
    std::vector<float *> xShifts;
    std::vector<float *> yShifts;
    std::vector<double *> intensities;
    
    // MARK: private housekeeping
    
    void initialiseZeros();
    void updateUnarrangedMidPoint();
    bool directionSanityCheck();
    
    // MARK: private calculations
    
    void spotCoordToRelativeVec(double unarrangedX, double unarrangedY,
                                vec *arrangedPos);
    void addToBasisChange(vec angles, MatrixPtr chosenMat = MatrixPtr());
    void fixBasisChange();

public:
    /* Initialise all variables to zero */
    Detector();
    
    /* For master panel: initialising using detector distance, beamX, beamY */
    /* Then set Detector::masterPanel to this pointer */
    Detector(double distance, double beamX, double beamY);
    
    /* For interpreting old panels.txt data */
    Detector(DetectorPtr parent, Coord arrangedTopLeft, Coord arrangedBottomRight,
             double angle, double offsetX, double offsetY, double gain);
    
    /* For children panels, just setting parent */
    /* Follow with parent->addChild to complete the process */
    Detector(DetectorPtr parent);
    
    /* using initialisation variables from CrystFEL */
    void initialise(Coord unarrangedTopLeft, Coord unarrangedBottomRight,
                    vec slowDir, vec fastDir, vec _arrangedTopLeft, bool lastIsMiddle = false, bool ghost = false);
    
    /* For detectors who have children and are not the master */
    Detector(DetectorPtr parent, vec arrangedMiddle, std::string tag);
    
    static void fullDescription();
    void description(bool propogate = false);
    
    void postInit();
    void resetNudgeBasis();
    void nudgeTiltAndStep(double *nudgeTiltX, double *nudgeTiltY, double *nudgeStep, double *interNudge = NULL);

	static double lookupCache(double distSq);

    static bool isActive()
    {
        if (detectorActive == -1)
        {
            detectorActive = FileParser::hasKey("DETECTOR_LIST");
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
    }
    
    void setFastDirection(double newX, double newY, double newZ)
    {
        fastDirection.h = newX;
        fastDirection.k = newY;
        fastDirection.l = newZ;
    }
    
    void setArrangedTopLeft(double newX, double newY, double newZ);
    vec getArrangedTopLeft();
    
    void setArrangedMidPoint(vec midpoint)
    {
        arrangedMidPoint = midpoint;
    }
    
    vec getArrangedMidPoint()
    {
        return arrangedMidPoint;
    }
    
    /* Scattered ray from crystal hits the centre of the panel
     * at a certain normal vector... */

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
    
    vec getRotatedSlowDirection();
    vec getRotatedFastDirection();
    
    double halfSlow();
    double halfFast();
    
    MatrixPtr calculateChangeOfBasis(vec *fastAxis, vec *slowAxis, MatrixPtr inv);
    void recalculateChangeOfBasis(vec *fastAxis, vec *slowAxis);
    
    MatrixPtr getChangeOfBasis()
    {
        return changeOfBasisMat;
    }
    
    void rotateAxisRecursive(bool fix = false);
    
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
        
        if (newMaster)
        {
            newMaster->setTag("master");
        }
    }

    /* Return resolution histogram */
    
    static CSVPtr resolutionHistogram();
    
    // MARK: apply rotations bottom up
    
    void updateCurrentRotation();

    // MARK: Fix midpoints of mid-level panels
    void fixMidpoints();
	vec getAverageMidpoint();
    
    void rearrangeCoord(std::pair<float, float> *aShift);
    
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
    DetectorPtr spotCoordForMiller(MillerPtr miller, double *xSpot, double *ySpot);
    DetectorPtr intersectionForMiller(MillerPtr miller, vec *intersection);
    
    /* Keeping track of all Millers */
    
    void addMillerCarefully(MillerPtr miller);

	void clearMillers()
	{
		millers.clear();
		std::vector<MillerWeakPtr>().swap(millers);

		for (int i = 0; i < childrenCount(); i++)
		{
			getChild(i)->clearMillers();
		}
	}

    size_t millerCount()
    {
        return millers.size();
    }
    
    MillerPtr miller(int i)
    {
        return millers[i].lock();
    }
    
    void removeMidPointRelativeToParent();

    /* Write geometry file */
    
    std::string writeGeometryFile(int fileCount, int indentCount = 0, CSVPtr differenceCSV = CSVPtr());
    std::string writeCrystFELFile();
    
    /* Refinement getter/setters */
    
    static void setArrangedMidPointX(void *object, double newX)
    {
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
        static_cast<Detector *>(object)->arrangedMidPoint.h = newX;
    }
    
    static double getArrangedMidPointX(void *object)
    {
        return static_cast<Detector *>(object)->arrangedMidPoint.h;
    }
    
    static void setArrangedMidPointY(void *object, double newY)
    {
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
        static_cast<Detector *>(object)->arrangedMidPoint.k = newY;
    }
    
    static double getArrangedMidPointY(void *object)
    {
        return static_cast<Detector *>(object)->arrangedMidPoint.k;
    }
    
    static void setArrangedMidPointZ(void *object, double newZ)
    {
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
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
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }
    
    static double getAlpha(void *object)
    {
        return static_cast<Detector *>(object)->rotationAngles.h;
    }
    
    static void setBeta(void *object, double newBeta)
    {
        static_cast<Detector *>(object)->rotationAngles.k = newBeta;
        static_cast<Detector *>(object)->updateCurrentRotation();
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }
    
    static double getBeta(void *object)
    {
        return static_cast<Detector *>(object)->rotationAngles.k;
    }
    
    static void setGamma(void *object, double newGamma)
    {
        static_cast<Detector *>(object)->rotationAngles.l = newGamma;
        static_cast<Detector *>(object)->updateCurrentRotation();
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }
    
    static double getGamma(void *object)
    {
        return static_cast<Detector *>(object)->rotationAngles.l;
    }
    
    static void setInterNudgeX(void *object, double horiz)
    {
        static_cast<Detector *>(object)->interNudge.h = horiz;
        static_cast<Detector *>(object)->updateCurrentRotation();
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }

    static double getInterNudgeX(void *object)
    {
        return static_cast<Detector *>(object)->interNudge.h;
    }
    
    static void setInterNudgeY(void *object, double vert)
    {
        static_cast<Detector *>(object)->interNudge.k = vert;
        static_cast<Detector *>(object)->updateCurrentRotation();
    }
    
    static double getInterNudgeY(void *object)
    {
        return static_cast<Detector *>(object)->interNudge.k;
    }
    
    static void setNudgeZ(void *object, double closeness)
    {
        static_cast<Detector *>(object)->nudgeTranslation.l = closeness;
        static_cast<Detector *>(object)->updateCurrentRotation();
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }
    
    static double getNudgeZ(void *object)
    {
        return static_cast<Detector *>(object)->nudgeTranslation.l;
    }

    static void setNudgeTiltX(void *object, double horizTilt)
    {
        static_cast<Detector *>(object)->nudgeRotation.h = horizTilt;
        static_cast<Detector *>(object)->updateCurrentRotation();
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }
    
    static void setSmartTiltX(void *object, double horizTilt)
    {
        Detector *det = static_cast<Detector *>(object);
        double angles = (det->nudgeRotation.k / det->nudgeTiltY
                         + det->nudgeRotation.h / det->nudgeTiltX);
        double distance = det->nudgeStep * angles;
        det->nudgeRotation.h = horizTilt;
        det->nudgeTranslation.l = distance;
        det->updateCurrentRotation();
    }
    
    static void setSmartTiltY(void *object, double vertTilt)
    {
        Detector *det = static_cast<Detector *>(object);
        det->nudgeRotation.k = vertTilt;
        double angles = (det->nudgeRotation.k / det->nudgeTiltY
                         + det->nudgeRotation.h / det->nudgeTiltX);
        double distance = det->nudgeStep * angles;
        det->nudgeTranslation.l = distance;
        det->updateCurrentRotation();
    }
    
    static double getNudgeTiltX(void *object)
    {
        return static_cast<Detector *>(object)->nudgeRotation.h;
    }
    
    static void setNudgeTiltY(void *object, double vertTilt)
    {
        static_cast<Detector *>(object)->nudgeRotation.k = vertTilt;
        static_cast<Detector *>(object)->updateCurrentRotation();
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }
    
    static double getNudgeTiltY(void *object)
    {
        return static_cast<Detector *>(object)->nudgeRotation.k;
    }
    
    static void setInterNudgeZ(void *object, double spin)
    {
        static_cast<Detector *>(object)->interNudge.l = spin;
        static_cast<Detector *>(object)->updateCurrentRotation();
        static_cast<Detector *>(object)->setUpdateMidPointForDetector();
    }
    
    static double getInterNudgeZ(void *object)
    {
        return static_cast<Detector *>(object)->nudgeRotation.l;
    }

    static double getPokeX(void *object)
    {
        return static_cast<Detector *>(object)->poke.h;
    }

    static double getPokeY(void *object)
    {
        return static_cast<Detector *>(object)->poke.k;
    }

    static double getPokeZ(void *object)
    {
        return static_cast<Detector *>(object)->poke.l;
    }

    void resetPoke();
    
    
    void setPokeN(int pokeNum, double pokeValue);
    
    static void setPokeX(void *object, double pokeValue)
    {
        static_cast<Detector *>(object)->setPokeN(0, pokeValue);
    }

    static void setPokeY(void *object, double pokeValue)
    {
        static_cast<Detector *>(object)->setPokeN(1, pokeValue);
    }

    static void setPokeZ(void *object, double pokeValue)
    {
        static_cast<Detector *>(object)->setPokeN(2, pokeValue);
    }

    static void setNoisy(bool _noisy)
    {
        noisy = _noisy;
    }
    
    static bool isNoisy()
    {
        return noisy;
    }
    
    static void setDrawImage(ImagePtr image)
    {
        if (!FileParser::getKey("DRAW_GEOMETRY_PNGS", true))
        {
            return;
        }
        
        drawImage = image;
    }
    
    /* MARK: Scoring functions */
    
    double millerScore(bool ascii = false, bool stdev = false, int number = -1);
    static double millerScoreWrapper(void *object);
    static double millerStdevScoreWrapper(void *object);
    
    /* Special image functions */
    
    static void drawSpecialImage(std::string filename = "");
    static ImagePtr getSpecialImage()
    {
        return drawImage;
    }
    
    /* Active index manager */
    
    void setIndexManager(IndexManagerPtr manager)
    {
        _manager = manager;
    }
    
    IndexManagerPtr getIndexManager()
    {
        return _manager;
    }
    
    double distanceFromSample();
    
    /* Resolution fun */
    
    double spotCoordToResolution(double unarrangedX, double unarrangedY, double wavelength);
    void resolutionOrZLimits(double *min, double *max, double wavelength, int type, CSVPtr csv = CSVPtr());
    void zLimits(double *min, double *max);
    void resolutionLimits(double *min, double *max, double wavelength);
    
    static void setDetectorType(DetectorType _type)
    {
        detectorType = _type;
    }
    
    static bool isMPCCD()
    {
        return (detectorType == DetectorTypeMPCCD);
    }

    static bool isCSPAD()
    {
        return (detectorType == DetectorTypeCSPAD);
    }
    
    vec midPointOffsetFromParent(bool useParent = true, bool resetNudge = false);
    void getAllSubDetectors(std::vector<DetectorPtr> &array, bool childrenOnly = false);
	std::vector<DetectorPtr> getSubDetectorsOnLevel(int level);
    
    void reportMillerScores(int refinementNum = 0);
    
    void prepareInterNudges();
    static void enableNudge();
    
    
    void setUpdateMidPointForDetector()
    {
        mustUpdateMidPoint = true;
        
        for (int i = 0; i < childrenCount(); i++)
        {
            getChild(i)->setUpdateMidPointForDetector();
        }
    }
    
    double getGain()
    {
        return gain;
    }
    
    void setGain(double _gain)
    {
        gain = _gain;
    }

    bool isRefinable(GeometryScoreType scoreType)
    {
        if (scoreType == GeometryScoreTypeInterMiller ||
            scoreType == GeometryScoreTypeIntraMiller)
        {
            for (int i = 0; i < childrenCount(); i++)
            {
                if (getChild(i)->isRefinable(GeometryScoreTypeBeamCentre))
                {
                    return false;
                }
            }

			return true;
        }
        
        if (scoreType == GeometryScoreTypeIntrapanelParent ||
            scoreType == GeometryScoreTypeBeamCentre)
        {
            return _refinable;
        }
        else if (scoreType == GeometryScoreTypeIntrapanel)
        {
            for (int i = 0; i < childrenCount(); i++)
            {
                if (getChild(i)->isRefinable(GeometryScoreTypeBeamCentre))
                {
                    return false;
                }
            }
            
            return _refinable;
        }
        else if (scoreType == GeometryScoreTypeInterpanel)
        {
            for (int i = 0; i < childrenCount(); i++)
            {
                if (getChild(i)->isRefinable(GeometryScoreTypeBeamCentre))
                {
                    return true;
                }
            }
            
            return false;
        }
        
        return false;
    }
    
    void setRefinable(bool refinable)
    {
        _refinable = refinable;
        logged << "Setting " << getTag() << " panel to " << (refinable ? "refineable" : "fixed") << std::endl;
        sendLog();
    }
};

#endif /* defined(__cppxfel__Detector__) */
