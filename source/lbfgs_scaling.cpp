#include "lbfgs_scaling.h"
#include <vector>
#include <string>
#include <iostream>
#include "stdio.h"
#include "definitions.h"

ScalingManager *scaling;


lbfgsfloatval_t Lbfgs_Scaling::evaluate(void *instance,
		const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
		const lbfgsfloatval_t step)
{
	(*scaling).set_Gs((lbfgsfloatval_t**) &x);

	//	scaling->evaluate_gradient(&g);

	lbfgsfloatval_t fx = 0;

	fx = scaling->rMerge();

	for (int i = 0; i < n; i++)
	{
		g[i] = scaling->gradientForL(i, fx);
	}

//	for (int i=0; i < n; i++)
//		std::cout << g[i] << std::endl;

	return fx;
}

static int progress(void *instance, const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step, int n, int k, int ls)
{
	/*	printf("Iteration %d:\n", k);
	 printf("  fx = %f, x[0] = %f, x[1] = %f, x[2] = %f\n", fx, x[0], x[1],
	 x[2]);
	 printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
	 printf("\n");*/
	return 0;
}

Lbfgs_Scaling::Lbfgs_Scaling(char **filenames, int filenum)
{
	scaling = new ScalingManager(filenames, filenum, 1);

	(*scaling).processReflections();
}

Lbfgs_Scaling::Lbfgs_Scaling(vector<MtzPtr> mtzs)
{
	scaling = new ScalingManager(mtzs);

	scaling->processReflections();
}

void Lbfgs_Scaling::run(void)
{
	int n = (*scaling).mtz_num;

	int ret = 0;
	lbfgsfloatval_t fx = 0;
	lbfgsfloatval_t *x = lbfgs_malloc(n);
	lbfgs_parameter_t param;

	if (x == NULL)
	{
		printf("ERROR: Failed to allocate a memory block for variables.\n");
		return;
	}

	std::cout << "n = " << n << std::endl;

    for (int i = 0; i < n; i++)
    {
        x[i] = 1;
    }
    
	fx = scaling->rMerge();
	std::cout << "R meas after gradient scaling: " << fx << std::endl;

	lbfgs_parameter_init(&param);
	param.max_iterations = 50;
    
    
    

	ret = lbfgs(n, x, &fx, evaluate, progress, NULL, &param);

	switch (ret)
	{
	case LBFGS_ALREADY_MINIMIZED:
		printf("LBFGS_ALREADY_MINIMIZED\n");
		break;
	case LBFGSERR_LOGICERROR:
		printf("LBFGSERR_LOGICERROR\n");
		break;
	case LBFGSERR_OUTOFMEMORY:
		printf("LBFGSERR_OUTOFMEMORY %i\n", ret);
		break;
	case LBFGSERR_CANCELED:
		printf("LBFGSERR_CANCELED %i\n", ret);
		break;
	case LBFGSERR_INVALID_N:
		printf("LBFGSERR_INVALID_N %i\n", ret);
		break;
	case LBFGSERR_INVALID_N_SSE:
		printf("LBFGSERR_INVALID_N_SSE %i\n", ret);
		break;
	case LBFGSERR_INVALID_X_SSE:
		printf("LBFGSERR_INVALID_X_SSE %i\n", ret);
		break;
	case LBFGSERR_INVALID_EPSILON:
		printf("LBFGSERR_INVALID_EPSILON %i\n", ret);
		break;
	case LBFGSERR_INVALID_TESTPERIOD:
		printf("LBFGSERR_INVALID_TESTPERIOD %i\n", ret);
		break;
	case LBFGSERR_ROUNDING_ERROR:
		printf("LBFGSERR_ROUNDING_ERROR %i\n", ret);
		break;
	case LBFGSERR_MINIMUMSTEP:
		printf("LBFGSERR_MINIMUMSTEP %i\n", ret);
		break;
	case LBFGSERR_MAXIMUMSTEP:
		printf("LBFGSERR_MAXIMUMSTEP %i\n", ret);
		break;
	case LBFGSERR_MAXIMUMLINESEARCH:
		printf("LBFGSERR_MAXIMUMLINESEARCH %i\n", ret);
		break;
	case LBFGSERR_MAXIMUMITERATION:
		printf("LBFGSERR_MAXIMUMITERATION %i\n", ret);
		break;
	default:
		break;
	}

	/* Report the result. */
	printf("L-BFGS optimization terminated with status code = %d\n", ret);
	//  printf(" fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);

	for (int i = 0; i < n; i++)
	{
        std::string filename = scaling->mtzs[i]->getFilename();

	//	std::cout << filename << " " << scaling->Gs[i].G << std::endl;

		if (scaling->Gs[i].G > 0.5)
			scaling->mtzs[i]->applyScaleFactor(scaling->Gs[i].G);
		else
			scaling->mtzs[i]->setRefCorrelation(0);
        
        x[i] = 1;
	}
    
    fx = scaling->rMerge();
    scaling->set_Gs(&x);
    
    std::cout << "R meas after correcting scales: " << fx << std::endl;

    double rsplit = scaling->rSplit();
    std::cout << "Estimated R split: " << rsplit << std::endl;


	lbfgs_free(x);
}

Lbfgs_Scaling::~Lbfgs_Scaling(void)
{
	delete scaling;
}
