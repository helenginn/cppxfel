/*
 * Spot.h
 *
 *  Created on: 30 Dec 2014
 *      Author: helenginn
 */

#ifndef SPOT_H_
#define SPOT_H_

#include <vector>
#include "Image.h"
#include "parameters.h"
#include "Vector.h"
#include "LoggableObject.h"

class Spot : LoggableObject, public boost::enable_shared_from_this<Spot>
{
private:
	vector<vector<double> > probe;
	ImageWeakPtr parentImage;
    DetectorWeakPtr lastDetector;
    double angleDetectorPlane;
    bool setAngle;
    bool checked;
    bool _isBeamCentre;
    int successfulCommonLines;
    double correctedX; double correctedY;
    double x; double y;
    double intensity;
    bool rejected;
    int height;
    int length;
    int background;
    int backgroundPadding;
    static double maxResolution;
    static double minIntensity;
    static double minCorrelation;
    static bool checkRes;
    static int useNewDetectorFormat;
    
public:
	Spot(ImagePtr image);
	virtual ~Spot();

    static void spotsAndVectorsToResolution(double lowRes, double highRes, std::vector<SpotPtr> spots, std::vector<SpotVectorPtr> spotVectors, std::vector<SpotPtr> *lowResSpots, std::vector<SpotVectorPtr> *lowResSpotVectors);
    void makeProbe(int height, int background, int size, int backPadding = 0);
    void addToMask(int *mask, int width, int height);
	void setXY(double x, double y);
    double scatteringAngle(ImagePtr image = ImagePtr());
	bool isAcceptable(ImagePtr image);
	double angleFromSpotToCentre(double centreX, double centreY);
    double angleInPlaneOfDetector(double centreX = 0, double centreY = 0, vec upBeam = new_vector(0, 1, 0));
    double resolution();
    bool isOnSameLineAsSpot(SpotPtr otherSpot, double toleranceDegrees);
    static void writeDatFromSpots(std::string filename, std::vector<SpotPtr> spots);
    bool isSameAs(SpotPtr spot2);
    double closeToSecondSpot(SpotPtr spot2, double squareMinDistance);
    double integrate();
    void recentreInWindow(int windowPadding = 0);

    Coord getXY();
    double getX(bool update = false);
    double getY(bool update = false);
    Coord getRawXY();
    vec estimatedVector();
    void setUpdate();
    double focusOnNearbySpot(double maxShift, double trialX, double trialY, int round = 0);
    
    std::string spotLine();
    
    void setRejected(bool isRejected = true)
    {
        rejected = isRejected;
    }
    
    bool isRejected()
    {
        return rejected;
    }
    
    int successfulLineCount()
    {
        return successfulCommonLines;
    }
    
    void addSuccessfulCommonLineCount()
    {
        successfulCommonLines++;
    }
    
    void deleteSuccessfulCommonLineCount()
    {
        successfulCommonLines--;
    }
    
    bool isChecked()
    {
        return checked;
    }
    
    void setChecked(bool newCheck = true)
    {
        checked = newCheck;
    }
    
    ImagePtr getParentImage()
	{
		return parentImage.lock();
	}

	void setParentImage(ImagePtr parentImage)
	{
		this->parentImage = parentImage;
	}
    
    DetectorPtr getDetector()
    {
        return lastDetector.lock();
    }
    
    void setDetector(DetectorPtr newD)
    {
        lastDetector = newD;
    }
    
    double getIntensity()
    {
        return intensity;
    }
    
    void setToBeamCentre()
    {
        _isBeamCentre = true;
    }
    
    bool isBeamCentre()
    {
        return _isBeamCentre;
    }
};

#endif /* SPOT_H_ */
