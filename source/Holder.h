/*
 * Holder.h
 *
 *  Created on: 25 Oct 2014
 *      Author: helenginn
 */

#ifndef HOLDER_H_
#define HOLDER_H_

#include "headers/csymlib.h"

#include "Miller.h"
#include "MtzGrouper.h"
#include "parameters.h"
#include "headers/csymlib.h"
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <cctbx/uctbx.h>

class Miller;

using cctbx::sgtbx::space_group;
using cctbx::sgtbx::space_group_symbols;
using cctbx::sgtbx::space_group_type;
using cctbx::sgtbx::reciprocal_space::asu;

class Holder
{
private:
	std::vector<MillerPtr> millers;
	int refl_id;
	int inv_refl_id;
	double refIntensity;
	double refSigma;
	double resolution;
    static space_group spaceGroup;
    int spgNum;
    static space_group_type spgType;
    static asu asymmetricUnit;
    static bool hasSetup;
    static bool setupUnitCell;
    int activeAmbiguity;
    std::vector<int> reflectionIds;
    static cctbx::uctbx::unit_cell unitCell;
public:
    Holder(float *unitCell = NULL, CSym::CCP4SPG *group = NULL);
    void setUnitCell(float *unitCell);
	virtual ~Holder();

	MillerPtr miller(int i);
    void printDescription();
	void addMiller(MillerPtr miller);
	int millerCount();
	Holder *copy(bool copyMillers);
	Holder *copy();
    
	static int indexForReflection(int h, int k, int l, CSym::CCP4SPG *lowspgroup, bool inverted);

    int checkOverlaps();
	void holderDescription();
	void calculateResolution(MtzManager *mtz);
	void clearMillers();
    void removeMiller(int index);

	double meanIntensityWithExclusion(std::string *filename);
	double meanIntensity();
	double meanSigma();
	double meanSigma(bool friedel);
	double meanPartiality();
	double mergedIntensity(WeightType weighting);
    double meanWeight();
    double rMergeContribution(double *numerator, double *denominator);

	bool betweenResolutions(double lowAngstroms, double highAngstroms);
	bool anyAccepted();
	int acceptedCount();
	int rejectCount();

	void merge(WeightType weighting, double *intensity, double *sigma, bool calculateRejections);
	double standardDeviation(WeightType weighting);
    
    void detailedDescription();
    double observedPartiality(MtzManager *reference, Miller *miller);

    double mergeSigma();
    
    int reflectionIdForMiller(cctbx::miller::index<> cctbxMiller);
    void generateReflectionIds();
    void setAdditionalWeight(double weight);
    
    asu *getAsymmetricUnit()
    {
        return &asymmetricUnit;
    }
    
    space_group *getSpaceGroup()
    {
        return &spaceGroup;
    }
    
    void incrementAmbiguity();
    int ambiguityCount();
    MatrixPtr matrixForAmbiguity(int i);

    void setSpaceGroup(int spaceGroupNum);
    void setSpaceGroup(CSym::CCP4SPG *ccp4spg, cctbx::sgtbx::space_group_type newSpgType, asu newAsymmetricUnit);
    
    void resetFlip();
    void setFlip(int i);
    void setFlipAsActiveAmbiguity();
    
    int getActiveAmbiguity()
    {
        return activeAmbiguity;
    }
    
    void setActiveAmbiguity(int newAmbiguity)
    {
        activeAmbiguity = newAmbiguity;
    }

	const std::vector<MillerPtr>& getMillers() const
	{
		return millers;
	}

	void setMillers(const std::vector<MillerPtr>& millers)
	{
		this->millers = millers;
	}

	double getRefIntensity() const
	{
		return refIntensity;
	}

	void setRefIntensity(double refIntensity)
	{
		this->refIntensity = refIntensity;
	}

    int getReflId()
    {
        if (activeAmbiguity > ambiguityCount())
        {
            std::cout << "Active ambiguity error" << std::endl;
        }
        
        return reflectionIds[activeAmbiguity];
    }
    
	double getRefSigma() const
	{
		return refSigma;
	}

	void setRefSigma(double refSigma)
	{
		this->refSigma = refSigma;
	}

	double getResolution() const
	{
		return resolution;
	}

	void setResolution(double resolution)
	{
		this->resolution = resolution;
	}
};

#endif /* HOLDER_H_ */
