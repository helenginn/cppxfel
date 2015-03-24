#include "lbfgs_cluster.h"
#include <vector>
#include <string>
#include <iostream>
#include "stdio.h"
#include <cmath>

#include "FileParser.h"
#include "Logger.h"

StatisticsManager *statsManager;

lbfgsfloatval_t Lbfgs_Cluster::evaluate(void *instance,
		const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
		const lbfgsfloatval_t step)
{
	lbfgsfloatval_t fx = 0;

	for (int i = 0; i <= n - 4; i += 2)
	{
		for (int j = i + 2; j <= n - 2; j += 2)
		{
			lbfgsfloatval_t correl = (*statsManager).cc_array[i / 2][j / 2];
			lbfgsfloatval_t inv_correl = 1
					- (*statsManager).inv_cc_array[i / 2][j / 2];
			// lbfgsfloatval_t hits = (*statsManager).hits[i/2][j/2];
			//		lbfgsfloatval_t inv_hits = (*statsManager).inv_hits[i/2][j/2];

			lbfgsfloatval_t dot_product = x[i] * x[j] + x[i + 1] * x[j + 1];

			if (correl != -1)
			{
				double addition = pow(correl - dot_product, 2);
				fx += addition;
				//		cout << "Correl 0/0 = " << addition << endl;
			}
			if (inv_correl != 2)
			{ // 1 - (-1)
				double addition = pow(inv_correl - dot_product, 2);
				fx += addition;
				//		cout << "Correl 0/1 = " << addition << endl;
			}
		}
	}

//	for (int i = 0; i < n; i += 2)
//		cout << "x[" << i << "] = " << x[i] << ", x[" << i + 1 << "] = "
//				<< x[i + 1] << endl;

//	cout << endl;

//	printf("f(x) = %.3f\n", fx);

	for (int i = 0; i <= n - 2; i += 2)
	{
		double sigma_i = 0;
		double sigma_j = 0;
		double ai = x[i];
		double aj = x[i + 1];

		for (int j = 0; j <= n - 2; j += 2)
		{
			if (i == j)
				continue;

			double correl = (*statsManager).cc_array[i / 2][j / 2];
			double inv_correl = 1 - (*statsManager).inv_cc_array[i / 2][j / 2];

			// lbfgsfloatval_t hits = (*statsManager).hits[i/2][j/2];
			//	lbfgsfloatval_t inv_hits = (*statsManager).inv_hits[i/2][j/2];

			double bi = x[j];
			double bj = x[j + 1];

			//	double magnitude_term = sqrt(ai * ai + aj * aj) * sqrt(bi * bi + bj * bj);
			double dot_product = (ai * bi + aj * bj);

			double left_hand_side = -2 * (correl - dot_product);
			double inv_left_hand_side = -2 * (inv_correl - dot_product);

			double right_hand_side = bi;

			double addition_i = left_hand_side * right_hand_side;
			double inv_addition_i = inv_left_hand_side * right_hand_side;

			//		cout << "For " << i << ": left_hand_side = " << left_hand_side << ", right_hand_side_left = " << right_hand_side_left << ", right_hand_side_right = " << right_hand_side_right << endl;

			// and for aj

			right_hand_side = bj;

			double addition_j = left_hand_side * right_hand_side;
			double inv_addition_j = inv_left_hand_side * right_hand_side;

			//		cout << "For " << i + 1 << ": left_hand_side = " << left_hand_side << ", right_hand_side_left = " << right_hand_side_left << ", right_hand_side_right = " << right_hand_side_right << endl;

			if ((*statsManager).cc_array[i / 2][j / 2] != -1)
			{
				sigma_i += addition_i;
				sigma_j += addition_j;

				//		cout << "Correl 0/0 gradient addition = " << addition_i << ", " << addition_j << endl;
			}
			if ((*statsManager).inv_cc_array[i / 2][j / 2] != -1)
			{
				sigma_i += inv_addition_i;
				sigma_j += inv_addition_j;

				//		cout << "Correl 0/1 gradient addition = " << inv_addition_i << ", " << inv_addition_j << endl;
			}
		}

		g[i] = sigma_i;
		g[i + 1] = sigma_j;

//		cout << "g[" << i << "] = " << g[i] << ", g[" << i + 1 << "] = "
//				<< g[i + 1] << endl;
	}

//	cout << endl;

	return fx;
}

lbfgsfloatval_t Lbfgs_Cluster::evaluate_angle_only(void *instance,
		const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
		lbfgsfloatval_t step)
{
	lbfgsfloatval_t fx = 0;

	for (int i = 0; i <= n - 4; i += 2)
	{
		for (int j = i + 2; j <= n - 2; j += 2)
		{
			lbfgsfloatval_t correl = (*statsManager).cc_array[i / 2][j / 2];
			lbfgsfloatval_t inv_correl = 1
					- (*statsManager).inv_cc_array[i / 2][j / 2];

			lbfgsfloatval_t dot_product = x[i] * x[j] + x[i + 1] * x[j + 1];
			lbfgsfloatval_t magnitude_term = sqrt(
					x[i] * x[i] + x[i + 1] * x[i + 1])
					* sqrt(x[j] * x[j] + x[j + 1] * x[j + 1]);

			if (correl != -1)
			{
				double addition = pow(correl - dot_product / magnitude_term, 2);
				fx += addition;
				//		cout << "Correl 0/0 = " << addition << endl;
			}
			if (inv_correl != 2)
			{
				double addition = pow(inv_correl - dot_product / magnitude_term,
						2);
				fx += addition;
				//		cout << "Correl 0/1 = " << addition << endl;
			}
		}
	}

	for (int i = 0; i < n; i += 2)
		cout << "x[" << i << "] = " << x[i] << ", x[" << i + 1 << "] = "
				<< x[i + 1] << endl;

	cout << endl;

	printf("f(x) = %.3f\n", fx);

	for (int i = 0; i <= n - 2; i += 2)
	{
		double sigma_i = 0;
		double sigma_j = 0;
		double ai = x[i];
		double aj = x[i + 1];

		for (int j = 0; j <= n - 2; j += 2)
		{
			if (i == j)
				continue;

			double correl = (*statsManager).cc_array[i / 2][j / 2];
			double inv_correl = 1 - (*statsManager).inv_cc_array[i / 2][j / 2];

			double bi = x[j];
			double bj = x[j + 1];

			double magnitude_term = sqrt(ai * ai + aj * aj)
					* sqrt(bi * bi + bj * bj);
			double dot_product = (ai * bi + aj * bj);

			double left_hand_side = -2 * (correl - dot_product / magnitude_term)
					/ pow(magnitude_term, 2);
			double inv_left_hand_side = -2
					* (inv_correl - dot_product / magnitude_term)
					/ pow(magnitude_term, 2);

			double right_hand_side_left = magnitude_term * bi;
			double right_hand_side_right = dot_product * ai
					* sqrt(bi * bi + bj * bj)
					* std::pow((double) (ai * ai + aj * aj), (double) -0.5);

			double addition_i = left_hand_side
					* (right_hand_side_left - right_hand_side_right);
			double inv_addition_i = inv_left_hand_side
					* (right_hand_side_left - right_hand_side_right);

			//		cout << "For " << i << ": left_hand_side = " << left_hand_side << ", right_hand_side_left = " << right_hand_side_left << ", right_hand_side_right = " << right_hand_side_right << endl;


			// and for aj

			right_hand_side_left = magnitude_term * bj;
			right_hand_side_right = dot_product * aj * sqrt(bi * bi + bj * bj)
					* std::pow((double) (ai * ai + aj * aj), (double) -0.5);

			double addition_j = left_hand_side
					* (right_hand_side_left - right_hand_side_right);
			double inv_addition_j = inv_left_hand_side
					* (right_hand_side_left - right_hand_side_right);

			//		cout << "For " << i + 1 << ": left_hand_side = " << left_hand_side << ", right_hand_side_left = " << right_hand_side_left << ", right_hand_side_right = " << right_hand_side_right << endl;

			if ((*statsManager).cc_array[i / 2][j / 2] != -1)
			{
				sigma_i += addition_i;
				sigma_j += addition_j;

				//		cout << "Correl 0/0 gradient addition = " << addition_i << ", " << addition_j << endl;
			}
			if ((*statsManager).inv_cc_array[i / 2][j / 2] != -1)
			{
				sigma_i += inv_addition_i;
				sigma_j += inv_addition_j;

				//		cout << "Correl 0/1 gradient addition = " << inv_addition_i << ", " << inv_addition_j << endl;
			}
		}

		g[i] = sigma_i;
		g[i + 1] = sigma_j;

		cout << "g[" << i << "] = " << g[i] << ", g[" << i + 1 << "] = "
				<< g[i + 1] << endl;
	}

//	cout << endl;

	return fx;
}

static int progress(void *instance, const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step, int n, int k, int ls)
{
    ostringstream logged;
    
    logged << "Iteration: " << k << std::endl;
    logged << "fx = " << fx << std::endl;
    
    Logger::mainLogger->addStream(&logged);
    
    return 0;
}

void Lbfgs_Cluster::initialise_cluster_lbfgs(char **files, int filenum,
		MtzManager **mergedMtz)
{
	statsManager = new StatisticsManager();

	statsManager->loadFiles(files, filenum, 1);
	statsManager->generate_cc_grid();

	run(mergedMtz);

	delete statsManager;
}

void Lbfgs_Cluster::initialise_cluster_lbfgs(vector<MtzPtr> mtzs,
		MtzManager **mergedMtz)
{
    mtz_num = (int)mtzs.size();
    
    statsManager = new StatisticsManager();
	statsManager->setMtzs(mtzs);

	statsManager->generate_cc_grid();

	run(mergedMtz);

	delete statsManager;
}

void Lbfgs_Cluster::run(MtzManager **mergedMtz)
{
	int n = mtz_num * 2;

    ostringstream logged;
    
	if (mtz_num <= 10)
	{
		for (int i = 0; i < mtz_num; i++)
		{
			for (int j = 0; j < mtz_num; j++)
            {
                logged << statsManager->cc_array[i][j] << "\t";
            }
            
            logged << std::endl;
        }
		
        logged << std::endl;

		for (int i = 0; i < mtz_num; i++)
		{
			for (int j = 0; j < mtz_num; j++)
            {
                logged << statsManager->inv_cc_array[i][j] << "\t";
            }
            
            logged << std::endl;
		}
        
        logged << std::endl;
    }
    
    Logger::mainLogger->addStream(&logged);
    logged.str("");

	int i, ret = 0;
	lbfgsfloatval_t fx;
	lbfgsfloatval_t *x = lbfgs_malloc(n);
	lbfgs_parameter_t param;

	if (x == NULL)
	{
		printf("ERROR: Failed to allocate a memory block for variables.\n");
		return;
	}

	/* Initialize the variables. */
	for (i = 0; i < n; i++)
	{
		//
		x[0] = 0;
		x[1] = 1;
		x[2] = 1;
		x[3] = 0;
		if (i > 3)
			x[i] = (double) rand() / (double) RAND_MAX;
	}

	lbfgs_parameter_init(&param);
//    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
	param.max_iterations = 200;

	/*
	 Start the L-BFGS optimization; this will invoke the callback functions
	 evaluate() and progress() when necessary.
	 */
	ret = lbfgs(n, x, &fx, evaluate, progress, NULL, &param);

	switch (ret)
    {
        case LBFGSERR_INVALID_N:
            logged << "LBFGSERR_INVALID_N";
            break;
        case LBFGSERR_INVALID_N_SSE:
            logged << "LBFGSERR_INVALID_N_SSE";
            break;
        case LBFGSERR_INVALID_X_SSE:
            logged << "LBFGSERR_INVALID_X_SSE";
            break;
        case LBFGSERR_ROUNDING_ERROR:
            logged << "LBFGSERR_ROUNDING_ERROR";
            break;
        case LBFGSERR_MINIMUMSTEP:
            logged << "LBFGSERR_MINIMUMSTEP";
            break;
        case LBFGSERR_MAXIMUMSTEP:
            logged << "LBFGSERR_MAXIMUMSTEP";
            break;
        case LBFGSERR_MAXIMUMLINESEARCH:
            logged << "LBFGSERR_MAXIMUMLINESEARCH";
            break;
        case LBFGSERR_MAXIMUMITERATION:
            logged << "LBFGSERR_MAXIMUMITERATION";
            break;
        default:
		break;
	}
    
    logged << std::endl;
    
	/* Report the result. */
    
    logged << "L-BFGS optimization terminated with status code = " << ret << std::endl;

    for (int i = 0; i < n; i += 2)
	{
        logged << statsManager->mtzs[i / 2]->getFilename() << "\t" << x[i] << "\t" << x[i + 1] << std::endl;
	}

    Logger::mainLogger->addStream(&logged);
    logged.str("");

	vector<MtzPtr> all;
	int idx = 0;
	int inv = 0;

	logged << "Flipping images" << std::endl;

	std::vector<std::pair<double, double> > dots =
			std::vector<std::pair<double, double> >();

	for (int i = 0; i < n; i += 2)
	{
		dots.push_back(std::make_pair(x[i], x[i + 1]));
	}

	double gradient = this->findDividingLine(dots);

	logged << "Best gradient: y = " << gradient << "x" << std::endl;

	double angle = atan(gradient);

	double wedgeRemoval = FileParser::getKey("REMOVE_WEDGE", 0);

	double wedgeAngle = wedgeRemoval * M_PI / 2;

	double leftAngle = angle + wedgeAngle / 2;
	double rightAngle = angle - wedgeAngle / 2;

	double leftGradient = tan(leftAngle);
	double rightGradient = tan(rightAngle);

	logged << "Upper gradient = " << leftGradient << ", lower gradient = " << rightGradient << std::endl;

	for (int i = 0; i < n; i += 2)
	{
		if (x[i + 1] > x[i] * leftGradient)
		{
			statsManager->mtzs[i / 2]->setInverse(true);
			idx++;
			MtzPtr ptr = MtzPtr();
			all.push_back(statsManager->mtzs[i / 2]);
		}
		else if (x[i + 1] < x[i] * rightGradient)
		{
			statsManager->mtzs[i / 2]->setInverse(false);
			inv++;
			all.push_back(statsManager->mtzs[i / 2]);
		}
		else
		{
			statsManager->mtzs[i / 2]->setRejected(true);
		}
	}

	int excluded = statsManager->mtz_num - (int)all.size();

	logged << "N: Excluded images: " << excluded << std::endl;
	logged << "N: idx count: " << idx << std::endl;
	logged << "N: inv count: " << inv << std::endl;

    Logger::mainLogger->addStream(&logged);
    logged.str("");
    
	for (int i=0; i < mtz_num; i++)
	{
        statsManager->mtzs[i]->applyUnrefinedPartiality();
//		statsManager->mtzs[i]->setPartialityToUnity();
	}

	if (mergedMtz != NULL)
	{
		logged << "Merging from cluster algorithm" << std::endl;

		MtzGrouper *idxGrouper = new MtzGrouper();
	//	idxGrouper->setExpectedResolution(20);
		idxGrouper->setWeighting(WeightTypeAverage);
		idxGrouper->setMtzManagers(all);
		idxGrouper->merge(mergedMtz);
		delete idxGrouper;
	}
    
	lbfgs_free(x);
}

double Lbfgs_Cluster::findDividingLine(std::vector<std::pair<double, double> > dots)
{
    ostringstream logged;
    
	double angle = 45 * M_PI / 180;
	int count = 0;
	int minimum = INT_MAX;
	double stepAngle = 10 * M_PI / 180;

	while (count < 50 && minimum > 2)
	{
		logged << "Count: " << count << ", angle: " << angle << ", minimum " << minimum << std::endl;

		int j = 0;
        int scores[3] = {0, 0, 0};
		double trials[3] = {0, 0, 0};

		for (double testAngle = angle - stepAngle; testAngle <= angle + stepAngle; testAngle +=
				stepAngle)
		{
			double gradient = tan(testAngle);

			int idxCount = 0;
			int invCount = 0;

			for (int i = 0; i < dots.size(); i++)
			{
				bool side = (dots[i].second > dots[i].first * gradient);
				idxCount += (side == true);
				invCount += (side == false);
			}

			scores[j] = abs(idxCount - invCount);
			trials[j] = testAngle;
			j++;

		}

		int bestTrial = 1;
		int bestScore = INT_MAX;

		for (int i=0; i < 3; i++)
		{
			if (scores[i] < bestScore)
			{
				bestScore = scores[i];
				bestTrial = i;
			}
		}

		if (bestTrial == 1)
			stepAngle /= 2;

		minimum = bestScore;
		angle = trials[bestTrial];

		count++;
	}

	logged << "Final minimum: " << minimum << std::endl;

    Logger::mainLogger->addStream(&logged);
    
	double gradient = tan(angle);

	return gradient;
}
