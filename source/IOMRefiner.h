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

class IOMRefiner : public boost::enable_shared_from_this<IOMRefiner>
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
	vector<double> unitCell;
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
	int expectedSpots;
	int searchSize;
    static bool lowIntensityPenalty;
	double maxResolution;
	MtzManager *reference;
	void calculateNearbyMillers();
	double lastTotal;
	double lastStdev;
	double testBandwidth;
    double lastScore;
	vector<vector<double> > solutions;
    bool recalculateMillerPositions;
    IndexingSolutionPtr indexingSolution;
    
    double orientationTolerance;
    std::ostringstream logged;
    void sendLog(LogLevel priority = LogLevelNormal);
    double hkScore(bool stdev, bool finalise = false);
    double lScore(bool silent);
    void recalculateMillers(bool thorough = false);
    static double hkScoreWrapper(void *object);
    static double hkScoreFinalise(void *object);
    static double hkScoreStdevWrapper(void *object);
    static double lScoreWrapper(void *object);
    double getReflectionWavelengthStdev();
    
public:
	IOMRefiner(ImagePtr newImage = ImagePtr(), MatrixPtr matrix = MatrixPtr());
    void setComplexMatrix();
    virtual ~IOMRefiner();

    void lockUnitCellDimensions();
    void calculateOnce();
	void checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox = false, bool perfectCalculation = true);
	MtzPtr newMtz(int i, bool silent = false);
	void getWavelengthHistogram(vector<double> &wavelengths,
			vector<int> &frequencies, LogLevel level = LogLevelDetailed, int whichAxis = 0);

	static void duplicateSpots(vector<ImagePtr>images);
	
	double minimizeParameter(double *meanStep, double *param, int whichAxis = 0);
	void minimizeTwoParameters(double *meanStep1, double *meanStep2,
			double *param1, double *param2);

	void refineOrientationMatrix();
	
    void showHistogram(bool silent);
    bool millerWithinBandwidth(MillerPtr miller);
	int getTotalReflections();
    int getTotalReflectionsWithinBandwidth();
    
    double getRot(int rotNum);
	void dropMillers();

    bool isGoodSolution();

    double getDetectorDistance();
    double getWavelength();
    std::string refinementSummary();
    static std::string refinementSummaryHeader();
    
    void getBestRots(double *rot1, double *rot2, double *rot3)
    {
        *rot1 = bestHRot;
        *rot2 = bestKRot;
        *rot3 = bestLRot;
    }
    
    MatrixPtr getMatrix()
    {
        return matrix;
    }
    
    void setMatrix(MatrixPtr matrix)
    {
        this->matrix = matrix;
    }
    
    void setMatrixCopy(MatrixPtr matrix)
    {
        this->matrix = matrix->copy();
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

	double getTestBandwidth() const
	{
		return testBandwidth;
	}

	void setTestBandwidth(double testBandwidth)
	{
		this->testBandwidth = testBandwidth;
	}

	double getTestSpotSize() const
	{
		return testSpotSize;
	}

	void setTestSpotSize(double testSpotSize)
	{
		this->testSpotSize = testSpotSize;
	}

	double getMaxResolution() const
	{
		return maxResolution;
	}

	void setMaxResolution(double maxResolution)
	{
		this->maxResolution = maxResolution;
	}

	int getSearchSize() const
	{
		return searchSize;
	}

	void setSearchSize(int searchSize)
	{
		this->searchSize = searchSize;
	}

	double getInitialStep() const
	{
		return initialStep;
	}

	void setInitialStep(double initialStep = 1)
	{
		this->initialStep = initialStep;
	}

	CCP4SPG*& getSpaceGroup()
	{
		return spaceGroup;
	}

	void setSpaceGroup(CCP4SPG*& spaceGroup)
	{
     //   if (this->spaceGroup != NULL)
     //       ccp4spg_free(&this->spaceGroup);
        
		this->spaceGroup = spaceGroup;
	}

	const vector<Spot *>& getSpots() const
	{
		return spots;
	}

	int getExpectedSpots() const
	{
		return expectedSpots;
	}

	void setExpectedSpots(int expectedSpots)
	{
		this->expectedSpots = expectedSpots;
	}

	const vector<vector<double> >& getSolutions() const
	{
		return solutions;
	}

	void setSolutions(const vector<vector<double> >& solutions)
	{
		this->solutions = solutions;
	}

	vector<double>& getUnitCell()
	{
		return unitCell;
	}

	void setUnitCell(vector<double>& unitCell)
	{
		this->unitCell = unitCell;
	}
    
    void setOrientationTolerance(double newTolerance)
    {
        orientationTolerance = newTolerance;
    }
    
    double getBestLRot()
    {
        return bestLRot;
    }
    
    bool isCalculatingRough()
    {
        return roughCalculation;
    }
    
    void setCalculatingRough(bool rough)
    {
        roughCalculation = rough;
    }
    
    MtzPtr getLastMtz()
    {
        return lastMtz;
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
