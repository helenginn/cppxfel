/*
 * MtzGrouper.cpp
 *
 *  Created on: 28 Aug 2014
 *      Author: helenginn
 */

#include "MtzGrouper.h"

#include <cmath>
#include "StatisticsManager.h"
#include "misc.h"
#include "headers/csymlib.h"
#include <boost/thread/thread.hpp>
#include "Miller.h"
#include "lbfgs_scaling.h"
#include "headers/ccp4_spg.h"
#include "headers/ccp4_general.h"
#include "headers/ccp4_parser.h"

#include "FileParser.h"
#include "GraphDrawer.h"

MtzGrouper::MtzGrouper()
{
	correlationThreshold = 0;
	excludeWorst = true;
	weighting = WeightTypeAverage;
	scalingType = ScalingTypeAverage;
	acceptableResolution = 1;
	cutResolution = false;
    expectedResolution = FileParser::getKey("MAX_RESOLUTION_ALL", 1.6);
}

MtzGrouper::~MtzGrouper()
{

}

bool MtzGrouper::isMtzAccepted(MtzPtr mtz)
{
    int minimumReflectionCutoff = FileParser::getKey(
                                                     "MINIMUM_REFLECTION_CUTOFF",
                                                     MINIMUM_REFLECTION_CUTOFF);

    double refCorrelation = mtz->getRefCorrelation();
    
    if ((refCorrelation < correlationThreshold && excludeWorst
         && !mtz->isFreePass())
        || mtz->isRejected())
        return false;
    
    if (refCorrelation < 0 || refCorrelation == 1)
        return false;
    
    if (mtz->accepted() < minimumReflectionCutoff)
        return false;
    
    return true;
}

void MtzGrouper::setMtzManagers(const std::vector<MtzPtr>& mtzManagers)
{
	this->mtzManagers = mtzManagers;

	double averageCorrelation = 0;

	for (int i = 0; i < mtzManagers.size(); i++)
	{
		averageCorrelation += mtzManagers[i]->getRefCorrelation();
	}

	averageCorrelation /= mtzManagers.size();

	this->setCorrelationThreshold(averageCorrelation - 0.05);
}

void MtzGrouper::merge(MtzManager **mergeMtz, MtzManager **unmergedMtz,
		int cycle, bool anom)
{
	logged << "N: ==== Merge cycle " << cycle << " ====" << std::endl;
	logged << "N: Scaling type: ";

	switch (scalingType)
	{
	case ScalingTypeAverage:
		logged << "Average merge" << std::endl;
		break;
	case ScalingTypeReference:
		logged << "Against reference (gradient)" << std::endl;
		break;
	case ScalingTypeReferenceLeastSquares:
		logged << "Against reference (least squares)" << std::endl;
		break;
	case ScalingTypeMinimizeRMerge:
		logged << "Minimize R merge" << std::endl;
		break;
	case ScalingTypeBFactor:
		logged << "B factor and scale" << std::endl;
		break;
	case ScalingTypeResolutionShells:
		logged << "Gradients per resolution shell" << std::endl;
		break;
	default:
		logged << "Unsupported option" << std::endl;
		break;
	}

	double averageCorrelation = 0;
	double averageAboveCutoff = 0;
	int aboveCutoffNum = 0;
	double rotationCorrection = 0;

	for (int i = 0; i < mtzManagers.size(); i++)
	{
		double correl = mtzManagers[i]->getRefCorrelation();
		averageCorrelation += correl;

		double hRot = mtzManagers[i]->getHRot();
		double kRot = mtzManagers[i]->getKRot();

		double correction = sqrt(hRot * hRot + kRot * kRot);
		rotationCorrection += correction;

		if (correl > correlationThreshold)
		{
			averageAboveCutoff += correl;
			aboveCutoffNum++;
		}

        double a, b, c;
        mtzManagers[i]->getMatrix()->orientationMatrixUnitCell(&a, &b, &c);
        
        double rSplit = 0;
        if (MtzManager::getReferenceManager() != NULL)
            rSplit = mtzManagers[i]->rSplit(0, expectedResolution);
        
		logged << mtzManagers[i]->getFilename() << "\t" << correl << "\t"
				<< mtzManagers[i]->accepted() << "\t"
				<< mtzManagers[i]->getMosaicity() << "\t"
				<< mtzManagers[i]->getWavelength() << "\t"
				<< mtzManagers[i]->getBandwidth() << "\t" << hRot << "\t"
				<< kRot << "\t" << mtzManagers[i]->getSpotSize() << "\t"
				<< mtzManagers[i]->getExponent() << "\t" << rSplit << "\t" << a << std::endl;
	}

	averageCorrelation /= mtzManagers.size();
	rotationCorrection /= mtzManagers.size();
	averageAboveCutoff /= aboveCutoffNum;

	logged << "N: Average correlation per image: " << averageCorrelation
			<< std::endl;
	logged << "N: Average correlation for those above threshold: "
			<< averageAboveCutoff << std::endl;
	logged << "N: Average rotation correction in degrees: "
    << rotationCorrection << std::endl;
    
    sendLog();
    
    if (MtzManager::getReferenceManager() != NULL)
    {
        double refScale = 1000 / MtzManager::getReferenceManager()->averageIntensity();
        MtzManager::getReferenceManager()->applyScaleFactor(refScale);
    }
    
	double scale = 1;

	for (int i = 0; i < mtzManagers.size(); i++)
	{
		if (scalingType == ScalingTypeAverage)
		{
			scale = 1000 / mtzManagers[i]->averageIntensity();
		}
		else if (scalingType == ScalingTypeReference
                 || scalingType == ScalingTypeReferenceLeastSquares
                 || scalingType == ScalingTypeMinimizeRMerge)
		{
			bool leastSquares =
					(scalingType == ScalingTypeReferenceLeastSquares);

			scale = mtzManagers[i]->gradientAgainstManager(
					*MtzManager::getReferenceManager(), leastSquares);
		}
		else if (scalingType == ScalingTypeBFactor)
		{
			double newScale = 1;
			double bFactor = 0;

			mtzManagers[i]->bFactorAndScale(&newScale, &bFactor);
			mtzManagers[i]->applyScaleFactor(newScale, bFactor);
		}
		else if (scalingType == ScalingTypeResolutionShells)
		{
			mtzManagers[i]->applyScaleFactorsForBins();
		}

		mtzManagers[i]->applyScaleFactor(scale);
	}
/*
#ifdef MAC
	GraphDrawer drawer = GraphDrawer(MtzManager::getReferenceManager());

	drawer.resolutionStatsPlot(mtzManagers, "resolution");
	drawer.resolutionStatsPlot(mtzManagers, "intensity_ref", GraphMap(), true, false);
	drawer.resolutionStatsPlot(mtzManagers, "intensity_img", GraphMap(), true, true);
#endif*/

	if (scalingType == ScalingTypeMinimizeRMerge)
	{
        std::vector<MtzPtr> acceptedMtzs;
        
        for (int i = 0; i < mtzManagers.size(); i++)
        {
            if (isMtzAccepted(mtzManagers[i]))
                acceptedMtzs.push_back(mtzManagers[i]);
        }
        
		Lbfgs_Scaling scaling = Lbfgs_Scaling(acceptedMtzs);
		scaling.run();
	}

	logged << "Altered scales." << std::endl;

	MtzManager *idxMerge = NULL;
	MtzManager **unmerged = NULL;
	MtzManager *invMerge = NULL;

    string idxName = string("idxMerge.mtz");
    string invName = string("invMerge.mtz");
    string unmergedName = string("unmerged.mtz");
    
    if (cycle >= 0)
    {
        unmergedName = string("unmerged") + i_to_str(cycle) + string(".mtz");
        idxName = string("idxMerge") + i_to_str(cycle) + string(".mtz");
        invName = string("invMerge") + i_to_str(cycle) + string(".mtz");
    }
    
	if (anom == false)
	{
        time_t startcputime;
        time(&startcputime);
        
        merge(mergeMtz, unmergedMtz, false, true, &unmergedName);
		
        time_t endcputime;
        time(&endcputime);
        
        time_t difference = endcputime - startcputime;
        double seconds = difference;
        
        int finalSeconds = (int) seconds % 60;
        int minutes = seconds / 60;
        
        std::cout << "N: Clock time " << minutes << " minutes, " << finalSeconds << " seconds to merge full data set" << std::endl;
        
        merge(&idxMerge, unmerged, true, false);
		merge(&invMerge, unmerged, false, false);
	}
	else
	{
		mergeAnomalous(mergeMtz, unmergedMtz, false, true);

		mergeAnomalous(&idxMerge, unmerged, true, false);
		mergeAnomalous(&invMerge, unmerged, false, false);
	}

	idxMerge->writeToFile(idxName, true);
	invMerge->writeToFile(invName, true);
    
	logged << "N: === R split ===" << std::endl;
    sendLog();
	idxMerge->rSplitWithManager(invMerge, false, false, 0, expectedResolution);
	logged << "N: === CC half ===" << std::endl;
    sendLog();
	idxMerge->correlationWithManager(invMerge, false, false, 0,
			expectedResolution);
    
    sendLog();

	delete idxMerge;
	delete invMerge;
}

void MtzGrouper::mergeWrapper(void *object, MtzManager **mergeMtz,
		MtzManager **unmergedMtz, bool firstHalf, bool all, bool anom)
{
	if (!anom)
	{
		static_cast<MtzGrouper *>(object)->merge(mergeMtz, unmergedMtz,
				firstHalf, all);
	}
	else
	{
		static_cast<MtzGrouper *>(object)->mergeAnomalous(mergeMtz, unmergedMtz,
				firstHalf, all);
	}
}

void MtzGrouper::merge(MtzManager **mergeMtz, MtzManager **unmergedMtz,
                       bool firstHalf, bool all, std::string *unmergedName)
{
	*mergeMtz = new MtzManager();
	(*mergeMtz)->setFilename("Merged");
	(*mergeMtz)->copySymmetryInformationFromManager(mtzManagers[0]);

	int start = firstHalf ? 0 : (int)mtzManagers.size() / 2;
	int end = firstHalf ? (int)mtzManagers.size() / 2 : (int)mtzManagers.size();

	if (all)
	{
		start = 0;
		end = (int)mtzManagers.size();
	}

	int mtzCount = groupMillers(mergeMtz, unmergedMtz, start, end);

	if (unmergedMtz != NULL)
	{
		(*unmergedMtz) = (*mergeMtz)->copy();
	}
    
    if (unmergedName != NULL)
    {
        (*mergeMtz)->writeToFile(*unmergedName);
    }

	logged << "Merging miller indices" << std::endl;

	bool outlier_rejection = FileParser::getKey("OUTLIER_REJECTION",
	REJECTING_MILLERS);

	mergeMillers(mergeMtz, all && outlier_rejection, mtzCount);

	unflipMtzs();
}

void MtzGrouper::mergeAnomalous(MtzManager **mergeMtz, MtzManager **unmergedMtz,
		bool firstHalf, bool all, std::string filename)
{
	MtzManager *positive = new MtzManager();
	MtzManager *negative = new MtzManager();
	positive->setFilename("positive_pair");
	negative->setFilename("negative_pair");
	positive->copySymmetryInformationFromManager(mtzManagers[0]);
	negative->copySymmetryInformationFromManager(mtzManagers[0]);

	*mergeMtz = new MtzManager();
	(*mergeMtz)->setFilename(filename);
	(*mergeMtz)->copySymmetryInformationFromManager(mtzManagers[0]);

	int start = firstHalf ? 0 : (int)mtzManagers.size() / 2;
	int end = firstHalf ? (int)mtzManagers.size() / 2 : (int)mtzManagers.size();

	if (all)
	{
		start = 0;
		end = (int)mtzManagers.size();
	}

	int mtzCount = groupMillersWithAnomalous(&positive, &negative, start, end);

	bool outlier_rejection = FileParser::getKey("OUTLIER_REJECTION", REJECTING_MILLERS);

	mergeMillers(&positive, outlier_rejection && all, mtzCount);
	mergeMillers(&negative, outlier_rejection && all, mtzCount);

	std::cout << "Merging miller indices" << std::endl;

	differenceBetweenMtzs(mergeMtz, &positive, &negative);

	if (all)
	{
		writeAnomalousMtz(&positive, &negative, "anomalous_diff.mtz");
	}

	delete negative;
	delete positive;

	(*mergeMtz)->description();

	unflipMtzs();
}

int MtzGrouper::groupMillers(MtzManager **mergeMtz, MtzManager **unmergedMtz,
		int start, int end)
{
	int flipCount = 0;
	int mtzCount = 0;

	for (int i = start; i < end; i++)
	{
		if (!isMtzAccepted(mtzManagers[i]))
            continue;
        
		if (mtzManagers[i]->isInverse())
		{
			mtzManagers[i]->setFlipped(true);
			flipCount++;
		}

		mtzCount++;

		double cutoffRes = 1 / acceptableResolution;

		if (cutResolution)
		{
			std::cout << "Cutoff res not supported anymore!" << std::endl;
		}

		for (int j = 0; j < mtzManagers[i]->holderCount(); j++)
		{
			if (mtzManagers[i]->holder(j)->getResolution() > cutoffRes)
				continue;

			int refl_id = mtzManagers[i]->holder(j)->getReflId();

			Holder *holder = NULL;
			(*mergeMtz)->findHolderWithId(refl_id, &holder);

			if (holder == NULL)
			{
				Holder *newHolder = mtzManagers[i]->holder(j)->copy(false);
				(*mergeMtz)->addHolder(newHolder);
			}
			else
			{
				for (int k = 0; k < mtzManagers[i]->holder(j)->millerCount();
						k++)
				{
                    MillerPtr newMiller = mtzManagers[i]->holder(j)->miller(k);
					holder->addMiller(newMiller);
					newMiller->setImageCorrelation(
							mtzManagers[i]->getRefCorrelation());
				}
			}
		}
	}

	int last_refl_id = 0;

	for (int i = 0; i < (*mergeMtz)->holderCount(); i++)
	{
		Holder *holder = (*mergeMtz)->holder(i);

		if (holder->getReflId() == last_refl_id)
		{
			std::cout << "Same" << std::endl;
		}
		if (holder->getReflId() < last_refl_id)
		{
			std::cout << "Less" << std::endl;
		}

		last_refl_id = holder->getReflId();
	}

	std::cout << "N: MTZs used in merge: " << mtzCount << std::endl;
	std::cout << "N: Total flipped: " << flipCount << std::endl;
	std::cout << "N: Holders used: " << (*mergeMtz)->holderCount() << std::endl;

	return mtzCount;
}

int MtzGrouper::groupMillersWithAnomalous(MtzManager **positive,
		MtzManager **negative, int start, int end)
{
	int flipCount = 0;
	int mtzCount = 0;

	for (int i = start; i < end; i++)
	{
		double refCorrelation = mtzManagers[i]->getRefCorrelation();

		if (refCorrelation < correlationThreshold && excludeWorst
				&& !mtzManagers[i]->isFreePass())
			continue;

		if (refCorrelation < 0 || refCorrelation == 1)
			continue;

		if (mtzManagers[i]->isInverse())
		{
			mtzManagers[i]->setFlipped(true);
			flipCount++;
		}

		mtzCount++;

		double cutoffRes = 1 / acceptableResolution;

		if (cutResolution)
		{
			std::cout << "Cutoff res not supported anymore!" << std::endl;
		}

		for (int j = 0; j < mtzManagers[i]->holderCount(); j++)
		{
			for (int k = 0; k < mtzManagers[i]->holder(j)->millerCount(); k++)
			{
				bool friedel =
						mtzManagers[i]->holder(j)->miller(k)->positiveFriedel();
				MtzManager *friedelMtz = (friedel ? *positive : *negative);

				if (mtzManagers[i]->holder(j)->getResolution() > cutoffRes)
					continue;

				Holder *holder = NULL;
				friedelMtz->findHolderWithId(
						mtzManagers[i]->holder(j)->getReflId(), &holder);

				// there is a bug here which would mis-sort friedel pairs if
				// there were opposing ones of the same symmetry in the same image

				if (holder == NULL)
				{
					Holder *newHolder = mtzManagers[i]->holder(j)->copy(true);
					newHolder->clearMillers();
					MillerPtr newMiller = mtzManagers[i]->holder(j)->miller(k);
					newHolder->addMiller(newMiller);
					friedelMtz->addHolder(newHolder);
					friedelMtz->sortLastHolder();
				}
				else
				{
					{
						MillerPtr newMiller = mtzManagers[i]->holder(j)->miller(
								k);
						holder->addMiller(newMiller);
						newMiller->setImageCorrelation(
								mtzManagers[i]->getRefCorrelation());
					}
				}
			}
		}
	}

	std::cout << "N: MTZs used in merge: " << mtzCount << std::endl;
	std::cout << "N: Total flipped: " << flipCount << std::endl;

	return mtzCount;
}

void MtzGrouper::mergeMillers(MtzManager **mergeMtz, bool reject, int mtzCount)
{
	int holderCount = 0;
	int millerCount = 0;
	int rejectCount = 0;
	double aveStdErr = 0;
	int aveStdErrCount = 0;
    
    bool recalculateSigma = FileParser::getKey("RECALCULATE_SIGMA", false);

	for (int i = 0; i < (*mergeMtz)->holderCount(); i++)
	{
		Holder *holder = (*mergeMtz)->holder(i);

		if (!holder->anyAccepted())
		{
			(*mergeMtz)->removeHolder(i);
			i--;
			continue;
		}

		double totalIntensity = 0;
		double totalStdev = 0;

		holder->merge(weighting, &totalIntensity, &totalStdev, reject);

		double error = totalStdev / totalIntensity;

		if (error == error && totalStdev != 100)
		{
			aveStdErr += abs(error);
			aveStdErrCount++;

		}

        double totalSigma = holder->meanSigma() / holder->meanPartiality();
        
        if (recalculateSigma)
        {
            int n = holder->acceptedCount();
            
            totalSigma = totalStdev / sqrt(n);
        }

		millerCount += holder->acceptedCount();
		rejectCount += holder->rejectCount();
		holderCount++;

		int h, k, l = 0;
		MillerPtr firstMiller = (*mergeMtz)->holder(i)->miller(0);
		ccp4spg_put_in_asu((*mergeMtz)->getLowGroup(), firstMiller->h,
				firstMiller->k, firstMiller->l, &h, &k, &l);

		MillerPtr newMiller = MillerPtr(new Miller((*mergeMtz), h, k, l));

		newMiller->setData(totalIntensity, totalSigma, 1, 0);
		newMiller->setParent((*mergeMtz)->holder(i));
		(*mergeMtz)->holder(i)->calculateResolution(*mergeMtz);

		(*mergeMtz)->holder(i)->clearMillers();
		(*mergeMtz)->holder(i)->addMiller(newMiller);
	}

	aveStdErr /= aveStdErrCount;

	double multiplicity = (double) millerCount / (double) holderCount;
	double aveRejection = (double) rejectCount / (double) mtzCount;

	std::cout << "N: Total MTZs: " << mtzManagers.size() << std::endl;
	std::cout << "N: Multiplicity before merge: " << multiplicity << std::endl;
	std::cout << "N: Rejects per image: " << aveRejection << std::endl;
	std::cout << "N: Average error per reflection: " << aveStdErr << std::endl;

	(*mergeMtz)->insertionSortHolders();
}

void MtzGrouper::writeAnomalousMtz(MtzManager **positive, MtzManager **negative,
		std::string filename)
{
    double doubleCell[6];
	float cell[6], wavelength, fdata[9];

	/* variables for symmetry */
	CCP4SPG *mtzspg = (*positive)->getLowGroup();
	float rsm[192][4][4];
	char ltypex[2];

	/* variables for MTZ data structure */
	MTZ *mtzout;
	MTZXTAL *xtal;
	MTZSET *set;
	MTZCOL *colout[9];

	/*  Removed: General CCP4 initializations e.g. HKLOUT on command line */

	(*positive)->getUnitCell(&doubleCell[0], &doubleCell[1], &doubleCell[2], &doubleCell[3], &doubleCell[4],
			&doubleCell[5]);
    
    for (int i = 0; i < 6; i++)
    {
        cell[i] = (float)doubleCell[i];
    }

	wavelength = (*positive)->getWavelength();

	mtzout = MtzMalloc(0, 0);
	ccp4_lwtitl(mtzout, "Anomalous dataset ", 0);
	mtzout->refs_in_memory = 0;
	mtzout->fileout = MtzOpenForWrite(filename.c_str());

// then add symm headers...
	for (int i = 0; i < mtzspg->nsymop; ++i)
		CCP4::rotandtrn_to_mat4(rsm[i], mtzspg->symop[i]);
	strncpy(ltypex, mtzspg->symbol_old, 1);
	ccp4_lwsymm(mtzout, mtzspg->nsymop, mtzspg->nsymop_prim, rsm, ltypex,
			mtzspg->spg_ccp4_num, mtzspg->symbol_old, mtzspg->point_group);

// then add xtals, datasets, cols
	xtal = MtzAddXtal(mtzout, "XFEL crystal", "XFEL project", cell);
	set = MtzAddDataset(mtzout, xtal, "Dataset", wavelength);
	colout[0] = MtzAddColumn(mtzout, set, "H", "H");
	colout[1] = MtzAddColumn(mtzout, set, "K", "H");
	colout[2] = MtzAddColumn(mtzout, set, "L", "H");
	colout[3] = MtzAddColumn(mtzout, set, "I(+)", "K");
	colout[4] = MtzAddColumn(mtzout, set, "SIGI(+)", "M");
	colout[5] = MtzAddColumn(mtzout, set, "I(-)", "K");
	colout[6] = MtzAddColumn(mtzout, set, "SIGI(-)", "M");
	colout[7] = MtzAddColumn(mtzout, set, "IMEAN", "J");
	colout[8] = MtzAddColumn(mtzout, set, "SIGIMEAN", "Q");

	int num = 0;

	int hits = 0;
	vector<Holder *> posHolders;
	vector<Holder *> negHolders;

	(*positive)->findCommonReflections(*negative, posHolders, negHolders, &hits,
			0);

	for (int i = 0; i < posHolders.size(); i++)
	{
		double intensityPlus = posHolders[i]->meanIntensity();
		double sigmaPlus = posHolders[i]->meanSigma();

		double intensityMinus = negHolders[i]->meanIntensity();
		double sigmaMinus = negHolders[i]->meanSigma();

		double intensity = (intensityPlus + intensityMinus) / 2;
		double sigma = (sigmaPlus + sigmaMinus) / 2;

		if (intensity != intensity)
		{
			continue;
		}

		num++;

		int h = negHolders[i]->miller(0)->h;
		int k = negHolders[i]->miller(0)->k;
		int l = negHolders[i]->miller(0)->l;
		int _h, _k, _l;
		ccp4spg_put_in_asu(mtzspg, h, k, l, &_h, &_k, &_l);

		fdata[0] = _h;
		fdata[1] = _k;
		fdata[2] = _l;
		fdata[3] = intensityPlus;
		fdata[4] = sigmaPlus;
		fdata[5] = intensityMinus;
		fdata[6] = sigmaMinus;
		fdata[7] = intensity;
		fdata[8] = sigma;
		ccp4_lwrefl(mtzout, fdata, colout, 9, num);
	}

// print header information, just for info
//	ccp4_lhprt(mtzout, 1);
	MtzPut(mtzout, " ");
	MtzFree(mtzout);
}

void MtzGrouper::differenceBetweenMtzs(MtzManager **mergeMtz,
		MtzManager **positive, MtzManager **negative)
{
	int hits = 0;
	vector<Holder *> posHolders;
	vector<Holder *> negHolders;

	(*positive)->findCommonReflections(*negative, posHolders, negHolders, &hits,
			0);

	for (int i = 0; i < posHolders.size(); i++)
	{
		Holder *newHolder = posHolders[i]->copy(true);
		double posInt = posHolders[i]->meanIntensity();
		double negInt = negHolders[i]->meanIntensity();

		double posSigma = posHolders[i]->meanSigma();
		double negSigma = negHolders[i]->meanSigma();

		newHolder->miller(0)->setRawIntensity(posInt - negInt);
		newHolder->miller(0)->setPartiality(1);
		newHolder->miller(0)->setSigma(posSigma + negSigma);

		(*mergeMtz)->addHolder(newHolder);
	}

	(*mergeMtz)->insertionSortHolders();
}

void MtzGrouper::unflipMtzs()
{
	for (int i = 0; i < mtzManagers.size(); i++)
	{
		mtzManagers[i]->setFlipped(false);
	}
}

void MtzGrouper::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

