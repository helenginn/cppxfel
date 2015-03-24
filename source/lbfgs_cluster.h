#ifndef lbfgs_cluster
#define lbfgs_cluster

#include <vector>
#include <string>
#include <iostream>
#include <liblbfgs.h>
#include "StatisticsManager.h"

using namespace std;


class Lbfgs_Cluster
{
private:
	int mtz_num;
public:
	//Lbfgs_Cluster(char **filenames, int filenum);
	static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	static lbfgsfloatval_t evaluate_angle_only(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
//	~Lbfgs_Cluster(void);

	void run(MtzManager **mergedMtz);
	void initialise_cluster_lbfgs(char **files, int filenum, MtzManager **mergedMtz);
	void initialise_cluster_lbfgs(vector<MtzPtr>mtzs, MtzManager **mergedMtz);
	double findDividingLine(std::vector<std::pair<double, double> > dots);

	int getMtzNum() const {
		return mtz_num;
	}

	void setMtzNum(int mtzNum) {
		mtz_num = mtzNum;
	}
};

#endif //lbfgs_cluster
