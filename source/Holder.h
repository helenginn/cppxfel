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

class Miller;


class Holder
{
private:
	std::vector<MillerPtr> millers;
	int refl_id;
	int inv_refl_id;
	double refIntensity;
	double refSigma;
	double resolution;
    bool inverse;
    
public:
	Holder();
	virtual ~Holder();

	MillerPtr miller(int i);
    void printDescription();
	void addMiller(MillerPtr miller);
	int millerCount();
	Holder *copy(bool copyMillers);
	Holder *copy();
	void flip();
	void flip(MtzManager *mtz);
	void flipMiller(MillerPtr miller, int spg_num);

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

	bool betweenResolutions(double lowAngstroms, double highAngstroms);
	bool anyAccepted();
	int acceptedCount();
	int rejectCount();

	void merge(WeightType weighting, double *intensity, double *sigma, bool calculateRejections);
	double standardDeviation(WeightType weighting);
    
    double observedPartiality(MtzManager *reference, Miller *miller);
    
	int getInvReflId() const
	{
		return inv_refl_id;
	}

	void setInvReflId(int invReflId)
	{
		inv_refl_id = invReflId;
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
    
    void setInverse(bool newInverse)
    {
        inverse = newInverse;
    }
    
    bool isInverse()
    {
        return inverse;
    }

	int getReflId() const
	{
        return refl_id;
	}

	void setReflId(int reflId)
	{
		refl_id = reflId;
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
