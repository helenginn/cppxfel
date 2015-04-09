#ifndef scalingmanager
#define scalingmanager

#include <vector>
#include <string>
#include <iostream>
#include "MtzManager.h"
#include <liblbfgs.h>
#include <memory>
#include "definitions.h"

struct Scale_factor;

#include "ReflectionManager.h"

using namespace std;

struct Scale_factor
{
	MtzManager *image;
	double G;
	vector<double> params;
	vector<double> multipliers;
	vector<int> *reflections;
    vector<ReflectionManager *>refManagers;
};

class ScalingManager
{
private:
	int orderFinalReflection(void);
	void evaluate_gradient_at_l(int l);
    CCP4SPG *group;

public:
	vector<ReflectionManager> *refs;
	vector<MtzPtr> mtzs;
	vector<Scale_factor> Gs;
	vector<double> gradients;
    std::pair<int, std::vector<ReflectionManager *> > reflectionsPerImage;

	int mtz_num;

	void findReflectionInVector(int refl_id, ReflectionManager **ref, int *l);
	void processReflections(void);
	void evaluate_gradient(lbfgsfloatval_t **gradients);
	double evaluate_psi();
	double set_Gs(lbfgsfloatval_t **new_Gs);
	void outputSimpleScaling(lbfgsfloatval_t **scales);
    static void rMergeThreaded(ScalingManager *self, int offset, double *numerator, double *denominator, int imageNumber);
	double rMerge(int imageNumber = -1);
    double rSplit(void);
	static double rMergeWrapper(void *object);
	double gradientForL(int l, double currentR);
	double gradientForParam(int l, int paramNum, double currentR);
	double multiplierForParam(int l, int paramNum, double currentR,
			double resolution);
	double multiplierForL(int l, double currentR);
    void mergedHolders(MtzManager *templateMtz, vector<Holder *> &holders, bool half, bool all);

	void minimize(void);
//	void minimizeLimited(void);
	void gridSearch(void);
	void referenceScaling(MtzManager *reference);
	double residualsForImage(MtzManager *reference, MtzManager *image,
			double scale);

	ScalingManager(vector<MtzPtr> mtzs);
	ScalingManager(char **filenames, int filenum, int partiality);
	~ScalingManager(void);
    
    void setGroup(CCP4SPG *newGroup)
    {
        group = newGroup;
    }
};

#endif

