//#ifndef lbfgs_scaling
//#define lbfgs_scaling

#include <vector>
#include <string>
#include <iostream>
#include "ScalingManager.h"
#include <liblbfgs.h>
#include "parameters.h"

class Lbfgs_Scaling
{

public:
	Lbfgs_Scaling(char **filenames, int filenum);
	static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	~Lbfgs_Scaling(void);
	void run(void);
	Lbfgs_Scaling(vector<MtzPtr>mtzs);
};

//#endif
