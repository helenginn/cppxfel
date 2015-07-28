/*
 * Spot.h
 *
 *  Created on: 30 Dec 2014
 *      Author: helenginn
 */

#ifndef SPOT_H_
#define SPOT_H_

#include <vector>
#include "parameters.h"
#include "Panel.h"

class Image;

class Spot
{
private:
	vector<vector<double> > probe;
	Image *parentImage;
    double angleDetectorPlane;
    bool setAngle;
    bool checked;

public:
	Spot(Image *image);
	virtual ~Spot();
	double x; double y;

	double weight();
	double maximumLift(Image *image, int x, int y, bool ignoreCovers);
	double maximumLift(Image *image, int x, int y);
	void makeProbe(int height, int size);
	void setXY(int x, int y);
	double scatteringAngle(Image *image = NULL);
	bool isAcceptable(Image *image);
	static void sortSpots(vector<Spot *> *spots);
	static bool spotComparison(Spot *a, Spot *b);
    double angleInPlaneOfDetector();
    bool isOnSameLineAsSpot(SpotPtr spot2, double tolerance);
    double resolution();
    static void writeDatFromSpots(std::string filename, std::vector<SpotPtr> spots);
    
    Coord getXY();
    double getX();
    double getY();
    Coord getRawXY();

    bool isChecked()
    {
        return checked;
    }
    
    void setChecked(bool newCheck = true)
    {
        checked = newCheck;
    }
    
    Image*& getParentImage()
	{
		return parentImage;
	}

	void setParentImage(Image*& parentImage)
	{
		this->parentImage = parentImage;
	}
    

};

#endif /* SPOT_H_ */
