/*
 * Holder.h
 *
 *  Created on: 25 Oct 2014
 *      Author: helenginn
 */

#ifndef REFLECTION_H_
#define REFLECTION_H_

#include "csymlib.h"

class Reflection;

#include "parameters.h"
#include "headers/csymlib.h"

class Miller;

struct LiteMiller
{
    double intensity;
    double weight;
    bool friedel;
} ;

class Reflection
{
private:
    
    std::vector<LiteMiller> liteMillers;
	vector<MillerPtr> millers;
	double resolution;
    float negativeFriedelIntensity;
    float positiveFriedelIntensity;
    static CSym::CCP4SPG *ccp4_space_group;
    static MatrixPtr unitCellMatrix;
    
    static unsigned char spgNum;
    static bool hasSetup;
    static std::mutex setupMutex;
    static bool setupUnitCell;
    static MatrixPtr customAmbiguity;
    static std::vector<MatrixPtr> flipMatrices;
    static double rejectSigma;
    static bool shouldReject;
    unsigned char activeAmbiguity;
    vector<unsigned int> reflectionIds;
    MutexPtr millerMutex;
public:
    Reflection(float *unitCell = NULL, CSym::CCP4SPG *group = NULL);
    void setUnitCell(float *unitCell);
	virtual ~Reflection();

    MillerPtr acceptedMiller(int num);
	MillerPtr miller(int i);
    void printDescription();
	void addMiller(MillerPtr miller);
    void addMillerCarefully(MillerPtr miller);
    void addLiteMiller(MillerPtr miller);
    
	int millerCount();
	ReflectionPtr copy(bool copyMillers = false);
    
	static int indexForReflection(int h, int k, int l, CSym::CCP4SPG *lowspgroup, bool inverted = false);
    static int reflectionIdForCoordinates(int h, int k, int l);
    
    int checkOverlaps();
    int checkSpotOverlaps(std::vector<SpotPtr> *spots, bool actuallyDelete = true);
    void reflectionDescription();
	void calculateResolution(MtzManager *mtz);
	void clearMillers();
    void removeMiller(int index);

	double meanIntensityWithExclusion(std::string *filename, int start = 0, int end = 0);
	double meanIntensity(bool withCutoff = true, int start = 0, int end = 0);
	double meanSigma();
	double meanSigma(bool friedel);
	double meanPartiality(bool withCutoff = true);
	double mergedIntensity(WeightType weighting);
    double meanWeight(bool cutoff = true);
    double rMergeContribution(double *numerator, double *denominator);

	bool betweenResolutions(double lowAngstroms, double highAngstroms);
	bool anyAccepted();
	int acceptedCount();
	int rejectCount();

	void merge(WeightType weighting, double *intensity, double *sigma, bool calculateRejections);
    void medianMerge(double *intensity, double *sigma, int *rejected, signed char friedel);
    void liteMerge(double *intensity, double *sigma, int *rejected, signed char friedel = -1);
    void clearLiteMillers();
	double standardDeviation(WeightType weighting);
    
    void detailedDescription();
    double observedPartiality(MtzManager *reference, Miller *miller);

    double mergeSigma();
    
    void generateReflectionIds();
    
    static CSym::CCP4SPG *getSpaceGroup()
    {
        return ccp4_space_group;
    }
    
    void incrementAmbiguity();
    static int ambiguityCount();
    static MatrixPtr matrixForAmbiguity(int i);

    static void setSpaceGroup(int spaceGroupNum);
    void setUnitCellDouble(double *theUnitCell);

    void resetFlip();
    void setFlip(int i);
    void setFlipAsActiveAmbiguity();
    
    LiteMiller liteMiller(int i)
    {
        return liteMillers[i];
    }
    
    int getActiveAmbiguity()
    {
        return activeAmbiguity;
    }
    
    void setActiveAmbiguity(int newAmbiguity)
    {
        activeAmbiguity = newAmbiguity;
    }

	const vector<MillerPtr>& getMillers() const
	{
		return millers;
	}

	void setMillers(const vector<MillerPtr>& millers)
	{
		this->millers = millers;
	}
    
    static MatrixPtr getCustomAmbiguity()
    {
        return customAmbiguity;
    }

    long unsigned int getReflId()
    {
   /*     if (activeAmbiguity > ambiguityCount())
        {
            std::cout << "Active ambiguity error" << std::endl;
        }
        */
        return reflectionIds[activeAmbiguity];
    }
    
	double getResolution() const
	{
		return resolution;
	}

	void setResolution(double resolution)
	{
		this->resolution = resolution;
	}
    
    int liteMillerCount()
    {
        return (int)liteMillers.size();
    }
    
    static MatrixPtr getFlipMatrix(int i);
    
    static bool reflLessThan(ReflectionPtr refl1, ReflectionPtr refl2)
    {
        return (refl1->getReflId() < refl2->getReflId());
    };
};

#endif /* REFLECTION_H_ */
