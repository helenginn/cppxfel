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

class Image;

class Spot
{
private:
	vector<vector<double> > probe;
	Image *parentImage;

public:
	Spot();
	Spot(Image *image);
	virtual ~Spot();
	int x; int y;

	double weight();
	double maximumLift(Image *image, int x, int y, bool ignoreCovers);
	double maximumLift(Image *image, int x, int y);
	void makeProbe(int height, int size);
	void setXY(int x, int y);
	double scatteringAngle(Image *image);
	bool isAcceptable(Image *image);
	static void sortSpots(vector<Spot *> *spots);
	static bool spotComparison(Spot *a, Spot *b);

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
