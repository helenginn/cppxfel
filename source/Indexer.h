/*
 * Indexer.h
 *
 *  Created on: 11 Nov 2014
 *      Author: helenginn
 */

#ifndef INDEXER_H_
#define INDEXER_H_

#include "Image.h"
#include "Miller.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "headers/csymlib.h"
#include "Spot.h"
#include "parameters.h"
#include "Logger.h"

class MtzManager;
class Image;
class Miller;

using namespace CSym;

typedef enum
{
    RefinementTypeOrientationMatrixEarly = 0,
	RefinementTypeDetectorWavelength = 1,
    RefinementTypeOrientationMatrixVeryEarly = 2,
    RefinementTypeOrientationMatrixLate = 3,
	RefinementTypeOrientationMatrixSpots = 4,
    RefinementTypeOrientationMatrixExactSpots = 5,
	RefinementTypeOrientationMatrixRough = 6,
	RefinementTypeOrientationMatrixMedian = 7,
    RefinementTypeOrientationMatrixTotalSignal = 8,
} RefinementType;

class Indexer
{
private:
	Image *image;
	std::vector<MillerPtr> millers;
	std::vector<MillerPtr> nearbyMillers;
	std::vector<Spot *>spots;
    MatrixPtr matrix;

    CCP4SPG *spaceGroup;
	std::vector<double> unitCell;
	double hRot;
	double kRot;
    double bestHRot;
    double bestKRot;
	int search;
	double testDistance;
	double testWavelength;
	double testSpotSize;
	double initialStep;
    double getTotalStandardDeviation();
	int expectedSpots;
	int searchSize;
	static double intensityThreshold;
    static bool absoluteIntensity;
	double maxResolution;
	MtzManager *reference;
	void calculateNearbyMillers(bool rough);
	RefinementType refinement;
	double lastTotal;
	double lastStdev;
	double testBandwidth;
    double lastScore;
	std::vector<std::vector<double> > solutions;
    
    std::ostringstream logged;
    void sendLog(LogLevel priority);

public:
	Indexer(Image *newImage = NULL, MatrixPtr matrix = MatrixPtr());
	virtual ~Indexer();

	void checkAllMillers(double maxResolution, double bandwidth, bool complexShoebox = false);
	MtzPtr newMtz(int i);
	void getWavelengthHistogram(std::vector<double> &wavelengths,
			std::vector<int> &frequencies, LogLevel level = LogLevelDetailed);
	double score();

	static bool millerReachesThreshold(MillerPtr miller);
	void findSpots();
	static void duplicateSpots(std::vector<Image *>images);
	void writeDatFromSpots(std::string filename);
	static void scatterSpots(std::vector<Image *> images);

	void matchMatrixToSpots();
	void matchMatrixToSpots(RefinementType refinement);
	void minimizeParameter(double *meanStep, double *param);
	void minimizeTwoParameters(double *meanStep1, double *meanStep2,
			double *param1, double *param2);

	void refineDetectorAndWavelength(MtzManager *reference = NULL);
	void refineOrientationMatrix();
	void refineOrientationMatrix(RefinementType refinementType);
	double refineRoundBeamAxis(double start, double end, double wedge, bool allSolutions);
	void refineRoundBeamAxis();

	int getTotalReflections();
	int getTotalReflections(double threshold);
    int getTotalReflectionsWithinBandwidth();
    double medianIntensity();
	int identicalSpotsAndMillers();
	double getTotalIntegratedSignal();

	void dropMillers();

    double getDetectorDistance();
    double getWavelength();
    
    void getBestRots(double *hRot, double *kRot)
    {
        *hRot = bestHRot;
        *kRot = bestKRot;
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
    
	Image*& getImage()
	{
		return image;
	}

	void setImage(Image*& image)
	{
		this->image = image;
	}

	int getSearch() const
	{
		return search;
	}

	void setSearch(int search)
	{
		this->search = search;
	}

	double getHRot() const
	{
		return hRot;
	}

	void setHRot(double rot)
	{
		hRot = rot;
	}

	double getKRot() const
	{
		return kRot;
	}

	void setKRot(double rot)
	{
		kRot = rot;
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
        if (this->spaceGroup != NULL)
            ccp4spg_free(&this->spaceGroup);
        
		this->spaceGroup = spaceGroup;
	}

	const std::vector<Spot *>& getSpots() const
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

	const std::vector<std::vector<double> >& getSolutions() const
	{
		return solutions;
	}

	void setSolutions(const std::vector<std::vector<double> >& solutions)
	{
		this->solutions = solutions;
	}

	double getIntensityThreshold() const
	{
		return intensityThreshold;
	}

	void setIntensityThreshold(double intensityThreshold)
	{
		this->intensityThreshold = intensityThreshold;
	}

	std::vector<double>& getUnitCell()
	{
		return unitCell;
	}

	void setUnitCell(std::vector<double>& unitCell)
	{
		this->unitCell = unitCell;
	}
};

#endif /* INDEXER_H_ */
