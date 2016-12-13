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
#include "csymlib.h"
#include <boost/thread/thread.hpp>
#include "Miller.h"
#include "Reflection.h"
#include "lbfgs_scaling.h"
#include "ccp4_general.h"
#include "ccp4_parser.h"
#include <fstream>

#include "FileReader.h"
#include "FileParser.h"
#include "GraphDrawer.h"
#include "FreeMillerLibrary.h"
#include "CSV.h"

MtzGrouper::MtzGrouper()
{
	correlationThreshold = 0;
	excludeWorst = true;
	weighting = WeightTypePartialitySigma;
	scalingType = ScalingTypeAverage;
	acceptableResolution = 1;
	cutResolution = false;
    expectedResolution = FileParser::getKey("MAX_RESOLUTION_ALL", 1.6);
    
    exclusionByCCHalf = FileParser::getKey("EXCLUSION_BY_CC_HALF", false);
        
}

MtzGrouper::~MtzGrouper()
{

}

bool MtzGrouper::isMtzAccepted(MtzPtr mtz)
{
    if (!excludeWorst)
        return true;
    
    int minimumReflectionCutoff = FileParser::getKey(
                                                     "MINIMUM_REFLECTION_CUTOFF",
                                                     MINIMUM_REFLECTION_CUTOFF);

    double refPartCorrelThreshold = FileParser::getKey(
                                                     "PARTIALITY_CORRELATION_THRESHOLD",
                                                     0.0);

    if (refPartCorrelThreshold > 0)
    {
        if (mtz->getRefPartCorrel() < refPartCorrelThreshold)
        {
            Logger::mainLogger->addString("Rejecting " + mtz->getFilename() + " due to low partiality correlation", LogLevelDetailed);
            return false;
        }
    }
    
    if (mtz->accepted() < minimumReflectionCutoff)
    {
        Logger::mainLogger->addString("Rejecting " + mtz->getFilename() + " due to not reaching minimum number of reflections", LogLevelDetailed);
        return false;
    }
    
    double refCorrelation = mtz->getRefCorrelation();
    
    if (refCorrelation < 0 || refCorrelation == 1)
    {
        Logger::mainLogger->addString("Rejecting " + mtz->getFilename() + " due to suspicious correlation with reference", LogLevelDetailed);
        
        return false;
    }
    
    double minimumRSplit = FileParser::getKey("R_FACTOR_THRESHOLD", 0.0);
    
    bool needsRSplit = mtz->getReferenceManager() != NULL;
    double rSplit = 0;
    
    if (needsRSplit)
        rSplit = mtz->rSplit(0, 0);
    
    if (refCorrelation < correlationThreshold)
    {
        Logger::mainLogger->addString("Rejecting " + mtz->getFilename() + " due to poor correlation with reference", LogLevelDetailed);
        
        return false;
    }
    
    if (mtz->isRejected())
    {
        Logger::mainLogger->addString("Rejecting " + mtz->getFilename() + " due to true rejection flag", LogLevelDetailed);
        
        return false;
    }
    
    if (needsRSplit && minimumRSplit > 0 && rSplit > minimumRSplit)
    {
        Logger::mainLogger->addString("Rejecting " + mtz->getFilename() + " due to R split being too high", LogLevelDetailed);
        
        return false;
    }
    
    return true;
}

void MtzGrouper::setMtzManagers(const vector<MtzPtr>& mtzManagers)
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

void MtzGrouper::merge(MtzPtr *mergeMtz, MtzPtr *unmergedMtz,
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
    
    sendLog();

	double averageCorrelation = 0;
	double averageAboveCutoff = 0;
	int aboveCutoffNum = 0;
	double rotationCorrection = 0;

    logged << "Filename\tCorrel\tRsplit\tPartcorrel\tRefcount\tMosaicity\tWavelength\tBandwidth\t";
    
    logged << "hRot\tkRot\t";
    
    logged << "rlpSize\texp\tcellA\tcellB\tcellC\tscale" << std::endl;

	for (int i = 0; i < mtzManagers.size(); i++)
    {
        MtzPtr mtz = mtzManagers[i];

		double correl = mtz->getRefCorrelation();
        if (correl != -1)
            averageCorrelation += correl;

		double hRot = mtz->getHRot();
		double kRot = mtz->getKRot();
        
		double correction = sqrt(hRot * hRot + kRot * kRot);
        
		rotationCorrection += correction;

		if (correl > correlationThreshold)
		{
			averageAboveCutoff += correl;
			aboveCutoffNum++;
		}

        double a, b, c;
        mtz->getMatrix()->orientationMatrixUnitCell(&a, &b, &c);
        
        double rSplit = 0;
        if (MtzManager::getReferenceManager() != NULL)
            rSplit = mtz->rSplit(0, expectedResolution, true, true);
        
        double partCorrel = mtz->getRefPartCorrel();
        
        double *cellDims = new double[3];
        mtz->getMatrix()->unitCellLengths(&cellDims);
        
        logged << mtzManagers[i]->getFilename() << "\t" << correl << "\t" << rSplit << "\t" << partCorrel << "\t"
				<< mtzManagers[i]->accepted() << "\t"
				<< mtzManagers[i]->getMosaicity() << "\t"
				<< mtzManagers[i]->getWavelength() << "\t"
				<< mtzManagers[i]->getBandwidth() << "\t";
        logged << hRot << "\t" << kRot << "\t";
        
        logged << mtzManagers[i]->getSpotSize() << "\t"
        << mtzManagers[i]->getExponent() << "\t" << cellDims[0] << "\t" << cellDims[1] << "\t" << cellDims[2] << "\t";
        
        logged << mtzManagers[i]->getScale() << std::endl;
        
        delete [] cellDims;
	}
    
    std::string tabbedParams = logged.str();
    std::replace(tabbedParams.begin(), tabbedParams.end(), '\t', ',');
    
    std::ofstream paramLog;
    std::string paramLogName = "params_cycle_" + i_to_str(cycle) + ".csv";
    
    std::string fullPath = FileReader::addOutputDirectory(paramLogName);
    
    paramLog.open(fullPath);
    paramLog << tabbedParams << std::endl;
    paramLog.close();
    
    logged << "Written parameter table to " << paramLogName << "." << std::endl;

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
    
	for (int i = 0; i < mtzManagers.size(); i++)
	{
        double scale = 1;

        if (scalingType == ScalingTypeAverage)
		{
			scale = 1000 / mtzManagers[i]->averageIntensity();
		}
		else if (scalingType == ScalingTypeReference
                 || scalingType == ScalingTypeMinimizeRMerge)
		{
            MtzManager *reference = MtzManager::getReferenceManager();
            
			scale = mtzManagers[i]->gradientAgainstManager(
					reference, false);
       //     std::cout << mtzManagers[i]->getFilename() << " " << scale << std::endl;
		}
		else if (scalingType == ScalingTypeBFactor)
		{
			double newScale = 1;
			double bFactor = 0;

			mtzManagers[i]->bFactorAndScale(&newScale, &bFactor);
		}
		else if (scalingType == ScalingTypeResolutionShells)
		{
			mtzManagers[i]->applyScaleFactorsForBins();
		}

		mtzManagers[i]->applyScaleFactor(scale);
	}

	if (scalingType == ScalingTypeMinimizeRMerge)
	{
        vector<MtzPtr> acceptedMtzs;
        
        for (int i = 0; i < mtzManagers.size(); i++)
        {
            if (isMtzAccepted(mtzManagers[i]))
                acceptedMtzs.push_back(mtzManagers[i]);
        }
        
		Lbfgs_Scaling scaling = Lbfgs_Scaling(acceptedMtzs);
		scaling.run();
	}

	logged << "Altered scales." << std::endl;

	MtzPtr idxMerge;
	MtzPtr unmerged;
	MtzPtr invMerge;

    std::string idxName = std::string("half1Merge.mtz");
    std::string invName = std::string("half2Merge.mtz");
    std::string unmergedName = std::string("unmerged.mtz");
    
    std::string anomString = anom ? "anom_" : "";
    std::string anomalousName = std::string("anomalous_diff.mtz");
    
    if (cycle >= 0)
    {
        anomalousName = std::string("anomalous_diff_") + i_to_str(cycle) + std::string(".mtz");;
        unmergedName = anomString + std::string("unmerged") + i_to_str(cycle) + std::string(".mtz");
        idxName = anomString + std::string("half1Merge") + i_to_str(cycle) + std::string(".mtz");
        invName = anomString + std::string("half2Merge") + i_to_str(cycle) + std::string(".mtz");
    }
    
	if (anom == false)
	{
        time_t startcputime;
        time(&startcputime);
        
        logged << "**** MERGING ALL DATA ****" << std::endl;
        sendLog();
        merge(mergeMtz, unmergedMtz, false, true, &unmergedName);
		
        time_t endcputime;
        time(&endcputime);
        
        time_t difference = endcputime - startcputime;
        double seconds = difference;
        
        int finalSeconds = (int) seconds % 60;
        int minutes = seconds / 60;
        
        logged << "N: Clock time " << minutes << " minutes, " << finalSeconds << " seconds to merge full data set" << std::endl;
        
        logged << "**** MERGING HALF DATA (1) ****" << std::endl;
        sendLog();
        merge(&idxMerge, &unmerged, true, false);

        logged << "**** MERGING HALF DATA (2) ****" << std::endl;
        sendLog();
        merge(&invMerge, &unmerged, false, false);
	}
	else
	{
		mergeAnomalous(mergeMtz, unmergedMtz, false, true, anomalousName);

		mergeAnomalous(&idxMerge, &unmerged, true, false);
		mergeAnomalous(&invMerge, &unmerged, false, false);
	}

	idxMerge->writeToFile(idxName, true);
	invMerge->writeToFile(invName, true);
    
	logged << "N: === R split ===" << std::endl;
    sendLog();
	idxMerge->rSplitWithManager(&*invMerge, false, false, 0, expectedResolution, 20, NULL, true);
	logged << "N: === CC half ===" << std::endl;
    sendLog();
	idxMerge->correlationWithManager(&*invMerge, false, false, 0,
			expectedResolution, 20, NULL, true);
    
    if (FreeMillerLibrary::active())
    {
        logged << "N: === R split (free only) ===" << std::endl;
        sendLog();
        idxMerge->rSplitWithManager(&*invMerge, false, false, 0, expectedResolution, 20, NULL, true, true);
        logged << "N: === CC half (free only) ===" << std::endl;
        sendLog();
        idxMerge->correlationWithManager(&*invMerge, false, false, 0,
                                         expectedResolution, 20, NULL, false, true);
    }
        
    sendLog();
}

void MtzGrouper::mergeWrapper(void *object, MtzPtr *mergeMtz,
		MtzPtr *unmergedMtz, bool firstHalf, bool all, bool anom)
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

void MtzGrouper::merge(MtzPtr *mergeMtz, MtzPtr *unmergedMtz,
                       bool firstHalf, bool all, std::string *unmergedName)
{
	*mergeMtz = MtzPtr(new MtzManager());
    
    std::string filename = std::string("merged_") + (all ? std::string("all_") : std::string("half_")) + (firstHalf ? std::string("0") : std::string("1"));
    
	(*mergeMtz)->setFilename(filename);
	(*mergeMtz)->copySymmetryInformationFromManager(mtzManagers[0]);

	int start = firstHalf ? 0 : (int)mtzManagers.size() / 2;
	int end = firstHalf ? (int)mtzManagers.size() / 2 : (int)mtzManagers.size();

	if (all)
	{
		start = 0;
		end = (int)mtzManagers.size();
	}

	int mtzCount = groupMillers(mergeMtz, unmergedMtz, start, end);
 //   std::cout << "N: Accepted " << total << " due to increase in CC half" << std::endl;
    
	if (!unmergedMtz)
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

void MtzGrouper::mergeAnomalous(MtzPtr *mergeMtz, MtzPtr *unmergedMtz,
		bool firstHalf, bool all, std::string filename)
{
	MtzPtr positive = MtzPtr(new MtzManager());
	MtzPtr negative = MtzPtr(new MtzManager());
	positive->setFilename("positive_pair");
	negative->setFilename("negative_pair");
	positive->copySymmetryInformationFromManager(mtzManagers[0]);
	negative->copySymmetryInformationFromManager(mtzManagers[0]);

	*mergeMtz = MtzPtr(new MtzManager());
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
		writeAnomalousMtz(&positive, &negative, filename);
	}

	(*mergeMtz)->description();

	unflipMtzs();
}

int MtzGrouper::groupMillers(MtzPtr *mergeMtz, MtzPtr *unmergedMtz,
		int start, int end)
{
	int mtzCount = 0;
    bool fastMerge = FileParser::getKey("FAST_MERGE", false);
    
    std::map<int, int> flipCounts;
    
    for (int i = 0; i < mtzManagers[0]->ambiguityCount(); i++)
    {
        flipCounts[i] = 0;
    }
    
	for (int i = start; i < end; i++)
	{
        if (!isMtzAccepted(mtzManagers[i]))
        {
            continue;
        }
        
        mtzManagers[i]->flipToActiveAmbiguity();
        int ambiguity = mtzManagers[i]->getActiveAmbiguity();
        flipCounts[ambiguity]++;

		mtzCount++;

		double cutoffRes = 1 / acceptableResolution;

		if (cutResolution)
		{
			std::cout << "Cutoff res not supported anymore!" << std::endl;
		}

		for (int j = 0; j < mtzManagers[i]->reflectionCount(); j++)
		{
            if (fastMerge && mtzManagers[i]->reflection(j)->acceptedCount() == 0 && mtzManagers[i]->reflection(j)->rejectCount() == 0)
                continue;
            
            if (mtzManagers[i]->reflection(j)->getResolution() > cutoffRes)
				continue;
            
			long unsigned int refl_id = mtzManagers[i]->reflection(j)->getReflId();

			ReflectionPtr reflection;
			(*mergeMtz)->findReflectionWithId(refl_id, &reflection);

			if (!reflection)
			{
				ReflectionPtr newReflection = mtzManagers[i]->reflection(j)->copy(false);
				(*mergeMtz)->addReflection(newReflection);
			}
			else
			{
				for (int k = 0; k < mtzManagers[i]->reflection(j)->millerCount();
						k++)
				{
                    MillerPtr newMiller = mtzManagers[i]->reflection(j)->miller(k);
                    
                    if (fastMerge && !newMiller->accepted() && !newMiller->isRejected())
                    {
                        continue;
                    }
					reflection->addMiller(newMiller);
                }
			}
		}
	}

	std::cout << "N: MTZs used in merge: " << mtzCount << std::endl;
	std::cout << "N: Flip ratios: ";
    
    for (std::map<int, int>::iterator it = flipCounts.begin(); it != flipCounts.end(); it++)
    {
        std::cout << flipCounts[it->first] << " ";
    }
    
    std::cout << std::endl;
    
	std::cout << "N: Reflections used: " << (*mergeMtz)->reflectionCount() << std::endl;

	return mtzCount;
}

int MtzGrouper::groupMillersWithAnomalous(MtzPtr *positive,
		MtzPtr *negative, int start, int end)
{
	int flipCount = 0;
	int mtzCount = 0;

    logged << "Grouping Millers for anomalous merge." << std::endl;
    sendLog();
    
	for (int i = start; i < end; i++)
	{
        if (!isMtzAccepted(mtzManagers[i]))
			continue;

        mtzManagers[i]->flipToActiveAmbiguity();

		mtzCount++;

		double cutoffRes = 1 / acceptableResolution;

		if (cutResolution)
		{
			std::cout << "Cutoff res not supported anymore!" << std::endl;
		}

		for (int j = 0; j < mtzManagers[i]->reflectionCount(); j++)
		{
			for (int k = 0; k < mtzManagers[i]->reflection(j)->millerCount(); k++)
			{
                bool friedel = false;
				mtzManagers[i]->reflection(j)->miller(k)->positiveFriedel(&friedel);
                
				MtzPtr friedelMtz = (friedel ? *positive : *negative);

				if (mtzManagers[i]->reflection(j)->getResolution() > cutoffRes)
					continue;

				ReflectionPtr reflection;
				friedelMtz->findReflectionWithId(
						mtzManagers[i]->reflection(j)->getReflId(), &reflection);

				if (!reflection)
				{
					ReflectionPtr newReflection = mtzManagers[i]->reflection(j)->copy(true);
					newReflection->clearMillers();
					MillerPtr newMiller = mtzManagers[i]->reflection(j)->miller(k);
					newReflection->addMiller(newMiller);
					friedelMtz->addReflection(newReflection);
					friedelMtz->sortLastReflection();
				}
				else
				{
					{
						MillerPtr newMiller = mtzManagers[i]->reflection(j)->miller(
								k);
						reflection->addMiller(newMiller);
					}
				}
			}
		}
	}

	std::cout << "N: MTZs used in merge: " << mtzCount << std::endl;
	std::cout << "N: Total flipped: " << flipCount << std::endl;

	return mtzCount;
}

void MtzGrouper::mergeMillers(MtzPtr *mergeMtz, bool reject, int mtzCount)
{
	int reflectionCount = 0;
	int millerCount = 0;
	int rejectCount = 0;
	double aveStdErr = 0;
	int aveStdErrCount = 0;
    
    bool recalculateSigma = FileParser::getKey("RECALCULATE_SIGMA", false);
    bool minimumMultiplicity = FileParser::getKey("MINIMUM_MULTIPLICITY", 0);
    
    for (int i = 0; i < (*mergeMtz)->reflectionCount(); i++)
	{
		ReflectionPtr reflection = (*mergeMtz)->reflection(i);
        double addedRejects = reflection->rejectCount();
        
        for (int j = 0; j < (*mergeMtz)->reflection(i)->millerCount(); j++)
        {
            MillerPtr miller = (*mergeMtz)->reflection(i)->miller(j);
            
            if (!miller->isRejected())
            {
                continue;
            }
            
       //     logged << miller->getH() << ", " << miller->getK() << ", " << miller->getL() << ", from " << miller->getMtzParent()->getFilename() << std::endl;
       //     sendLog();

        }
        
        int accepted = reflection->acceptedCount();
        
        if (accepted <= minimumMultiplicity || accepted == 0)
        {
            (*mergeMtz)->removeReflection(i);
            rejectCount += addedRejects;
            i--;
            continue;
        }
        
        sendLog();
        double totalIntensity = 0;
		double totalStdev = 0;

		reflection->merge(weighting, &totalIntensity, &totalStdev, reject);

        rejectCount += reflection->rejectCount();
		double error = totalStdev / totalIntensity;

		if (error == error && totalStdev != 100)
		{
			aveStdErr += fabs(error);
			aveStdErrCount++;

		}

        double totalSigma = 0;
        
        if (!recalculateSigma)
        {
            totalSigma = reflection->meanSigma() / reflection->meanPartiality();
        }
        else
        {
            totalSigma = reflection->mergeSigma();
            
            if (totalSigma == 0)
            {
                totalSigma = sqrt(fabs(totalIntensity));
            }
            else if (std::isnan(totalSigma))
            {
                totalSigma = sqrt(fabs(totalIntensity)); // anything reasonable
            }
            
            if (std::isnan(totalIntensity)) // give up
            {
                continue;
            }
        }
        
		millerCount += reflection->acceptedCount();
        
        reflectionCount++;

		int h, k, l = 0;
		MillerPtr firstMiller = (*mergeMtz)->reflection(i)->miller(0);
		ccp4spg_put_in_asu((*mergeMtz)->getLowGroup(), firstMiller->getH(),
				firstMiller->getK(), firstMiller->getL(), &h, &k, &l);

		MillerPtr newMiller = MillerPtr(new Miller(&*(*mergeMtz), h, k, l));

        if (totalSigma != totalSigma || totalIntensity != totalIntensity)
        {
            continue;
        }
        
		newMiller->setData(totalIntensity, totalSigma, 1, 0);
		newMiller->setParent((*mergeMtz)->reflection(i));
		(*mergeMtz)->reflection(i)->calculateResolution(&**mergeMtz);

		(*mergeMtz)->reflection(i)->clearMillers();
		(*mergeMtz)->reflection(i)->addMiller(newMiller);
	}

	aveStdErr /= aveStdErrCount;

	double multiplicity = (double) millerCount / (double) reflectionCount;
    double aveRejection = (double) rejectCount;// / (double) mtzCount;

	std::cout << "N: Total MTZs: " << mtzManagers.size() << std::endl;
	std::cout << "N: Multiplicity before merge: " << multiplicity << std::endl;
	std::cout << "N: Rejects per image: " << aveRejection << std::endl;
	std::cout << "N: Average error per reflection: " << aveStdErr << std::endl;

	(*mergeMtz)->insertionSortReflections();
}

void MtzGrouper::writeAnomalousMtz(MtzPtr *positive, MtzPtr *negative,
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
    
    std::string fullPath = FileReader::addOutputDirectory(filename);

	mtzout = MtzMalloc(0, 0);
	ccp4_lwtitl(mtzout, "Anomalous dataset ", 0);
	mtzout->refs_in_memory = 0;
	mtzout->fileout = MtzOpenForWrite(fullPath.c_str());

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
	vector<ReflectionPtr> posReflections;
	vector<ReflectionPtr> negReflections;

	(*positive)->findCommonReflections(&**negative, posReflections, negReflections, &hits);

	for (int i = 0; i < posReflections.size(); i++)
	{
		double intensityPlus = posReflections[i]->meanIntensity();
		double sigmaPlus = posReflections[i]->meanSigma();

		double intensityMinus = negReflections[i]->meanIntensity();
		double sigmaMinus = negReflections[i]->meanSigma();

        double intensity = (intensityPlus + intensityMinus) / 2;
		double sigma = (sigmaPlus + sigmaMinus) / 2;
        
        if (intensityPlus != intensityPlus || intensityMinus != intensityMinus)
            continue;
        
        if (sigmaPlus < 0.0001)
            sigmaPlus = nan(" ");

        if (sigmaMinus < 0.0001)
            sigmaMinus = nan(" ");

        if (intensityPlus != intensityPlus)
        {
            intensity = intensityMinus;
            sigma = sigmaMinus;
        }
        
        if (intensityMinus != intensityMinus)
        {
            intensity = intensityPlus;
            sigma = sigmaPlus;
        }

		num++;

		int h = negReflections[i]->miller(0)->getH();
		int k = negReflections[i]->miller(0)->getK();
		int l = negReflections[i]->miller(0)->getL();
		int _h, _k, _l;
		ccp4spg_put_in_asu(mtzspg, h, k, l, &_h, &_k, &_l);

		fdata[0] = _h;
		fdata[1] = _k;
		fdata[2] = _l;
		fdata[3] = intensityMinus;
		fdata[4] = sigmaMinus;
		fdata[5] = intensityPlus;
		fdata[6] = sigmaPlus;
		fdata[7] = intensity;
		fdata[8] = sigma;
		ccp4_lwrefl(mtzout, fdata, colout, 9, num);
	}

// print header information, just for info
//	ccp4_lhprt(mtzout, 1);
	MtzPut(mtzout, " ");
	MtzFree(mtzout);
}

void MtzGrouper::differenceBetweenMtzs(MtzPtr *mergeMtz,
		MtzPtr *positive, MtzPtr *negative)
{
	int hits = 0;
	vector<ReflectionPtr> posReflections;
	vector<ReflectionPtr> negReflections;

	(*positive)->findCommonReflections(&**negative, posReflections, negReflections, &hits);

	for (int i = 0; i < posReflections.size(); i++)
	{
		ReflectionPtr newReflection = posReflections[i]->copy(true);
		double posInt = posReflections[i]->meanIntensity();
		double negInt = negReflections[i]->meanIntensity();

		double posSigma = posReflections[i]->meanSigma();
		double negSigma = negReflections[i]->meanSigma();

		newReflection->miller(0)->setRawIntensity(posInt - negInt);
		newReflection->miller(0)->setPartiality(1);
		newReflection->miller(0)->setSigma(posSigma + negSigma);

		(*mergeMtz)->addReflection(newReflection);
	}

	(*mergeMtz)->insertionSortReflections();
}

void MtzGrouper::unflipMtzs()
{
	for (int i = 0; i < mtzManagers.size(); i++)
	{
        mtzManagers[i]->resetFlip();
	}
}

void MtzGrouper::sendLog(LogLevel priority)
{
    Logger::mainLogger->addStream(&logged, priority);
    logged.str("");
    logged.clear();
}

