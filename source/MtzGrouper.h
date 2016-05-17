/*
 * MtzGrouper.h
 *
 *  Created on: 28 Aug 2014
 *      Author: helenginn
 */

#ifndef MTZGROUPER_H_
#define MTZGROUPER_H_

#include <vector>

#include "MtzManager.h"
#include "parameters.h"
#include "Logger.h"

class MtzGrouper
{
private:
    bool usingNewRefinement;
    std::ostringstream logged;
	double correlationThreshold;
	ScalingType scalingType;
	bool excludeWorst;
	WeightType weighting;
	double acceptableResolution;
	bool cutResolution;
	double expectedResolution;
    bool isMtzAccepted(MtzPtr mtz);
    bool exclusionByCCHalf;

    static void checkCCHalf(vector<MtzPtr> *managers, int offset, int *total);
	void merge(MtzPtr *mergeMtz, MtzPtr *unmergedMtz, bool firstHalf,
			bool all, std::string *unmergedName = NULL);
public:
	MtzGrouper();
	virtual ~MtzGrouper();

	vector<MtzPtr> mtzManagers;

    void merge(MtzPtr *mergeMtz, MtzPtr *unmergedMtz = MtzPtr(), int cycle = -1, bool anom = false);

	void mergeAnomalous(MtzPtr *mergeMtz, MtzPtr *unmergedMtz,
			bool firstHalf, bool all, std::string filename = "anomalous_diff");
	void differenceBetweenMtzs(MtzPtr *mergeMtz, MtzPtr *positive, MtzPtr *negative);

	static void mergeWrapper(void *object, MtzPtr *mergeMtz,
			MtzPtr *unmergedMtz, bool firstHalf, bool all, bool anom);
	int groupMillers(MtzPtr *mergeMtz, MtzPtr *unmergedMtz, int start,
			int end);
	void mergeMillers(MtzPtr *mergeMtz, bool reject, int mtzCount);
	void unflipMtzs();

	int groupMillersWithAnomalous(MtzPtr *positive, MtzPtr *negative, int start,
			int end);
	void writeAnomalousMtz(MtzPtr *positive, MtzPtr *negative,
			std::string filename);
    
    void sendLog(LogLevel priority = LogLevelNormal);

	double getCorrelationThreshold() const
	{
		return correlationThreshold;
	}

	void setCorrelationThreshold(double correlationThreshold)
	{
		this->correlationThreshold = correlationThreshold;
	}

	const vector<MtzPtr>& getMtzManagers() const
	{
		return mtzManagers;
	}

	void setMtzManagers(const vector<MtzPtr>& mtzManagers);

	bool isExcludeWorst() const
	{
		return excludeWorst;
	}

	void setExcludeWorst(bool excludeWorst)
	{
		this->excludeWorst = excludeWorst;
	}

	WeightType getWeighting() const
	{
		return weighting;
	}

	void setWeighting(WeightType weighting)
	{
		this->weighting = weighting;
	}

	ScalingType getScalingType() const
	{
		return scalingType;
	}

	void setScalingType(ScalingType scalingType)
	{
		this->scalingType = scalingType;
	}

	double getAcceptableResolution() const
	{
		return acceptableResolution;
	}

	bool isCutResolution() const
	{
		return cutResolution;
	}

	void setCutResolution(bool cutResolution)
	{
		this->cutResolution = cutResolution;
	}

	double getExpectedResolution() const
	{
		return expectedResolution;
	}

	void setExpectedResolution(double expectedResolution)
	{
		this->expectedResolution = expectedResolution;
	}
};

#endif /* MTZGROUPER_H_ */
