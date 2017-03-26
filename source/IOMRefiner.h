/*
 * IOMRefiner.h
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#ifndef IOMRefiner_H_
#define IOMRefiner_H_

#include "Image.h"
#include "parameters.h"
#include "IndexManager.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "csymlib.h"
#include "Spot.h"
#include "Logger.h"

class Miller;

using namespace CSym;

class IOMRefiner : public hasSymmetry, public LoggableObject,
public boost::enable_shared_from_this<IOMRefiner>
{
private:
	ImageWeakPtr image;
	vector<MillerPtr> millers;
	vector<MillerPtr> nearbyMillers;
    vector<MillerPtr> roughMillers;
	vector<Spot *>spots;
    MatrixPtr matrix;
    MatrixPtr lastRotatedMatrix;
    MtzPtr lastMtz;

    double minResolution;
    bool roughCalculation;
    bool needsReintegrating;
    CCP4SPG *spaceGroup;
    bool complexUnitCell;
	double hRot;
	double kRot;
    double lRot;
    double bestHRot;
    double bestKRot;
    double bestLRot;
	double testWavelength;
	double testSpotSize;
	double initialStep;
    double getTotalStandardDeviation();
	int searchSize;
	double maxResolution;
	void calculateNearbyMillers();
	double lastTotal;
	double lastStdev;
	double testBandwidth;
    double lastScore;
    bool recalculateMillerPositions;
    IndexingSolutionPtr indexingSolution;
    
    double orientationTolerance;
    double hkScore(bool stdev, bool finalise = false);
    double lScore(bool silent);
    static double hkScoreFinalise(void *object);
    static double hkScoreStdevWrapper(void *object);
    static double lScoreWrapper(void *object);
    double getReflectionWavelengthStdev();

public:
	IOMRefiner(ImagePtr newImage = ImagePtr(), MatrixPtr matrix = MatrixPtr());
    virtual ~IOMRefiner();

    void checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox = false, bool perfectCalculation = true, bool fake = false);
	void recalculateMillers(bool thorough = false);
	MtzPtr newMtz(int i, bool silent = false);
	void getWavelengthHistogram(vector<double> &wavelengths,
			vector<int> &frequencies, LogLevel level = LogLevelDetailed, int whichAxis = 0);

	void refineOrientationMatrix();
    void fakeSpots();
	
    bool millerWithinBandwidth(MillerPtr miller);
	int getTotalReflections();
    
    void dropMillers();

    bool isGoodSolution();

    double getWavelength();
    std::string refinementSummary();
    static std::string refinementSummaryHeader();
    
    MatrixPtr getMatrix()
    {
        return matrix;
    }
    
    void setMatrix(MatrixPtr matrix)
    {
        this->matrix = matrix;
    }
    
    double getLastScore()
    {
        return lastScore;
    }
    
	ImagePtr getImage()
	{
		return image.lock();
	}

	void setImage(ImagePtr image)
	{
		this->image = image;
	}

	static double getHRot(void *object)
	{
		return static_cast<IOMRefiner *>(object)->hRot;
	}

	static void setHRot(void *object, double rot)
	{
		static_cast<IOMRefiner *>(object)->hRot = rot;
	}

    static double getKRot(void *object)
    {
        return static_cast<IOMRefiner *>(object)->kRot;
    }
    
    static void setKRot(void *object, double rot)
    {
        static_cast<IOMRefiner *>(object)->kRot = rot;
    }
    
    static double getLRot(void *object)
    {
        return static_cast<IOMRefiner *>(object)->lRot;
    }
    
    static void setLRot(void *object, double rot)
    {
        static_cast<IOMRefiner *>(object)->lRot = rot;
    }
    
    double getBestLRot()
    {
        return bestLRot;
    }

	double getSearchSize()
	{
		return searchSize;
	}

    IndexingSolutionPtr getIndexingSolution()
    {
        return indexingSolution;
    }
    
    void setIndexingSolution(IndexingSolutionPtr solution)
    {
        indexingSolution = solution;
    }
};

#endif /* IOMRefiner_H_ */
