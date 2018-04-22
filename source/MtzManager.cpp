#include "MtzManager.h"

#include <string>
#include <iostream>
#include "cmtzlib.h"
#include "csymlib.h"
#include "ccp4_spg.h"
#include "ccp4_general.h"
#include "Vector.h"
#include "FileReader.h"
#include <cerrno>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/variant.hpp>
#include "StatisticsManager.h"
#include "Reflection.h"
#include "Miller.h"
#include "Detector.h"
#include "Image.h"
#include "Hdf5Image.h"
#include "Hdf5ManagerProcessing.h"
#include "RefinementStepSearch.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "definitions.h"
#include "FileParser.h"
#include "CSV.h"

using namespace CMtz;

using namespace CSym;

vector<double> MtzManager::superGaussianTable;
bool MtzManager::setupSuperGaussian = false;
std::mutex MtzManager::tableMutex;
double MtzManager::superGaussianScale = 0;

MtzManager *MtzManager::referenceManager = NULL;
MtzPtr MtzManager::differenceManager = MtzPtr();

std::string MtzManager::describeScoreType()
{
	switch (scoreType)
	{
		case ScoreTypeCorrelation:
			return std::string("correl");
		case ScoreTypeMinimizeRSplit:
			return std::string("rfactor");
		case ScoreTypeReward:
			return std::string("reward");
		case ScoreTypeRSplitIntensity:
			return std::string("r+int");
        case ScoreTypeReferenceFree:
            return std::string("refree");
        case ScoreTypePartFree:
            return std::string("partfree");
        case ScoreTypePartRefScale:
            return std::string("partrefscale");
        case ScoreTypeGeometric:
            return std::string("geometric");
		default:
			return std::string("unknown");
	}
}

void MtzManager::makeSuperGaussianLookupTable(double exponent)
{
	if (setupSuperGaussian)
		return;

	tableMutex.lock();

	if (setupSuperGaussian)
	{
		tableMutex.unlock();
		return;
	}

	superGaussianScale = pow((M_PI / 2), (2 / exponent - 1));

	superGaussianTable.clear();

	const double min = 0;
	const double max = MAX_SUPER_GAUSSIAN;
	const double step = SUPER_GAUSSIAN_STEP;
	int count = 0;

	for (double x = min; x < max; x += step)
	{
		double value = super_gaussian(x, 0, superGaussianScale, exponent);

		//      std::cout << x << "\t" << value << std::endl;

		superGaussianTable.push_back(value);

		count++;
	}

	setupSuperGaussian = true;
	tableMutex.unlock();
}

bool MtzManager::reflection_comparison(ReflectionPtr i, ReflectionPtr j)
{
	return (i->getReflId() < j->getReflId());
}

void MtzManager::sortReflections()
{
	std::sort(reflections.begin(), reflections.end(), reflection_comparison);
}

void MtzManager::loadParametersMap()
{
	optimisingWavelength = FileParser::getKey("OPTIMISING_WAVELENGTH",
											  !OPTIMISED_WAVELENGTH);
	optimisingBandwidth = FileParser::getKey("OPTIMISING_BANDWIDTH",
											 !OPTIMISED_BANDWIDTH);
	optimisingMosaicity = FileParser::getKey("OPTIMISING_MOSAICITY",
											 !OPTIMISED_MOSAICITY);
	optimisingOrientation = FileParser::getKey("OPTIMISING_ORIENTATION",
											   !OPTIMISED_ROT);
	optimisingRlpSize = FileParser::getKey("OPTIMISING_RLP_SIZE",
										   !OPTIMISED_SPOT_SIZE);
	optimisingExponent = FileParser::getKey("OPTIMISING_EXPONENT",
											!OPTIMISED_EXPONENT);

	resetDefaultParameters();

	stepSizeWavelength = FileParser::getKey("STEP_SIZE_WAVELENGTH",
											MEAN_STEP);
	stepSizeBandwidth = FileParser::getKey("STEP_SIZE_BANDWIDTH",
										   BANDWIDTH_STEP);
	stepSizeMosaicity = FileParser::getKey("STEP_SIZE_MOSAICITY",
										   MOS_STEP);
	stepSizeRlpSize = FileParser::getKey("STEP_SIZE_RLP_SIZE", SPOT_STEP);
	stepSizeOrientation = FileParser::getKey("STEP_SIZE_ORIENTATION",
											 ROT_STEP);
	stepSizeOrientABC = FileParser::getKey("STEP_SIZE_ORIENTATION_ABC",
										   ROT_STEP * 3);
	stepSizeExponent = FileParser::getKey("STEP_SIZE_EXPONENT",
										  EXPONENT_STEP);

	toleranceWavelength = FileParser::getKey("TOLERANCE_WAVELENGTH",
											 MEAN_TOLERANCE);
	toleranceBandwidth = FileParser::getKey("TOLERANCE_BANDWIDTH",
											BANDWIDTH_TOLERANCE);
	toleranceMosaicity = FileParser::getKey("TOLERANCE_MOSAICITY",
											MOS_TOLERANCE);
	toleranceRlpSize = FileParser::getKey("TOLERANCE_RLP_SIZE",
										  SPOT_SIZE_TOLERANCE);
	toleranceOrientation = FileParser::getKey("TOLERANCE_ORIENTATION",
											  ROT_TOLERANCE);
	toleranceExponent = FileParser::getKey("TOLERANCE_EXPONENT",
										   EXPO_TOLERANCE);

	bandwidth = FileParser::getKey("INITIAL_BANDWIDTH", INITIAL_BANDWIDTH);
	mosaicity = FileParser::getKey("INITIAL_MOSAICITY", INITIAL_MOSAICITY);
	spotSize = FileParser::getKey("INITIAL_RLP_SIZE", INITIAL_SPOT_SIZE);
	exponent = FileParser::getKey("INITIAL_EXPONENT", INITIAL_EXPONENT);

	hRot = 0;
	kRot = 0;

	int defaultScoreInt = FileParser::getKey("DEFAULT_TARGET_FUNCTION",
											 (int) DEFAULT_SCORE_TYPE);
	defaultScoreType = (ScoreType) defaultScoreInt;

	setScale(1);
	setActiveAmbiguity(0);
    maxResolutionAll = FileParser::getKey("MAX_REFINED_RESOLUTION",
                                          MAX_OPTIMISATION_RESOLUTION);
    
    

    
	minResolutionAll = FileParser::getKey("MIN_REFINED_RESOLUTION", 0.);
}

MtzManager::MtzManager()
{
	refineOrientations = FileParser::getKey("REFINE_ORIENTATIONS", true);

	testBandwidth = FileParser::getKey("OVER_PRED_BANDWIDTH",
									   OVER_PRED_BANDWIDTH) / 2;
	testSpotSize = FileParser::getKey("OVER_PRED_RLP_SIZE",
									  OVER_PRED_SPOT_SIZE);

	lastTotal = 0;
	lastStdev = 0;
	fullyLoaded = false;
	_bin = 0;

	initialStep = FileParser::getKey("INITIAL_ORIENTATION_STEP", INITIAL_ORIENTATION_STEP);

	searchSize = FileParser::getKey("METROLOGY_SEARCH_SIZE",
									METROLOGY_SEARCH_SIZE);
	orientationTolerance = FileParser::getKey("INDEXING_ORIENTATION_TOLERANCE",
											  INDEXING_ORIENTATION_TOLERANCE);


	lastReference = NULL;
	reflections.resize(0);
	bandwidth = INITIAL_BANDWIDTH;
	hRot = 0;
	kRot = 0;
	lRot = 0;
	mosaicity = INITIAL_MOSAICITY;
	spotSize = INITIAL_SPOT_SIZE;
	wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.);
	refCorrelation = 0;
	exponent = INITIAL_EXPONENT;
	scoreType = ScoreTypeCorrelation;
	trust = TrustLevelBad;

    //bool checkRes = checkMaxResolution();
    //logged << "checkRes: " << checkRes << std::endl;
    //sendLog();
    
    //maxResolutionAll = getImagePtr()->Image::getFoundMaxResolution();
    //bool checkpointer = !getImagePtr();
    
    //logged << "Firstly, setting maximum resolution to found resolution" << maxResolutionAll << std::endl;
    //logged << "HELLOS! Setting maximum resolution to resolution" << maxResolutionAll << std::endl;
    //sendLog();

    maxResolutionAll = FileParser::getKey("MAX_REFINED_RESOLUTION", MAX_OPTIMISATION_RESOLUTION);
                                              //MAX_OPTIMISATION_RESOLUTION);

    logged << "HELLOS2! Setting maximum resolution to resolution" << maxResolutionAll << std::endl;
    sendLog();
    
	minResolutionAll = FileParser::getKey("MIN_REFINED_RESOLUTION", 0.);
	defaultScoreType = DEFAULT_SCORE_TYPE;
	rejected = false;
	scale = 1;
	externalScale = -1;
	previousReference = NULL;
	previousAmbiguity = -1;
	activeAmbiguity = 0;
	bFactor = 0;
	setInitialValues = false;
	refPartCorrel = 0;
	dropped = false;
	lastRSplit = 0;
	timeDelay = 0;
	_xPos = 0;
	_yPos = 0;

	loadParametersMap();

	matrix = MatrixPtr();
}

void MtzManager::clearReflections()
{
	reflections.clear();
	vector<ReflectionPtr>().swap(reflections);
}

void MtzManager::removeReflection(int i)
{
	reflections.erase(reflections.begin() + i);
}

void MtzManager::addReflection(ReflectionPtr reflection)
{
	size_t lowestId = 0;
	ReflectionPtr refl = findReflectionWithId(reflection, &lowestId);

	reflections.insert(reflections.begin() + lowestId, reflection);
}

void MtzManager::addMiller(MillerPtr miller)
{
	miller->setMtzParent(this);

	int reflection_index = Reflection::indexForReflection(miller->getH(), miller->getK(), miller->getL(), getSpaceGroup());
	ReflectionPtr prevReflection;

	this->findReflectionWithId(reflection_index, &prevReflection);

	miller->setImage(getImagePtr());

	if (prevReflection != NULL)
	{
		prevReflection->addMiller(miller); // TODO
		miller->setParent(prevReflection);
	}
	else
	{
		ReflectionPtr newReflection = ReflectionPtr(new Reflection());
		newReflection->setUnitCell(getUnitCell());
		newReflection->setSpaceGroup(getSpaceGroupNum());
		newReflection->addMiller(miller);
		miller->setParent(newReflection);
		newReflection->calculateResolution(this);
		addReflection(newReflection);
	}
}

int MtzManager::millerCount()
{
	int sum = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		sum += reflection(i)->millerCount();
	}

	return sum;
}

std::vector<MillerPtr> MtzManager::strongMillers()
{
	std::vector<MillerPtr> strongs;

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);

			if (miller->reachesThreshold())
			{
				strongs.push_back(miller);
			}
		}
	}

	return strongs;
}

void MtzManager::setMatrix(MatrixPtr newMat)
{
	matrix = newMat;
	rotatedMatrix = matrix->copy();
}

void MtzManager::setDefaultMatrix()
{
	std::vector<double> unitCell = getUnitCell();
	MatrixPtr mat = Matrix::matrixFromUnitCell(unitCell);
	matrix = mat->copy();
}

void MtzManager::copySymmetryInformationFromManager(MtzPtr toCopy)
{
	std::vector<double> otherUnitCell =  toCopy->getUnitCell();
	this->setUnitCell(otherUnitCell);

	this->setSpaceGroup(toCopy->getSpaceGroup());
}

void MtzManager::dropReflections()
{
	writeToFile("tmp-" + getFilename());

	reflections.clear();
	std::vector<ReflectionPtr>().swap(reflections);

	dropped = true;
}

void MtzManager::getWavelengthFromHDF5()
{
	if (!dropped)
	{
		bool useHdf5Wavelength = FileParser::getKey("USE_HDF5_WAVELENGTH", true);

		if (!useHdf5Wavelength)
		{
			wavelength = 0;
			return;
		}

		if (!Hdf5ManagerCheetah::cheetahManagerCount())
		{
			return;
		}

		unsigned long length = getFilename().length();
		std::string noImg = getFilename().substr(4, length - 4); // [img-tag-XXXXX_0 to tag-XXXXX_0]
		unsigned long underscoreIndex = noImg.rfind("_");
		std::string noCryst = noImg.substr(0, underscoreIndex);

		Hdf5ManagerCheetahPtr manager = Hdf5ManagerCheetah::hdf5ManagerForImage(noCryst);
		double aWavelength = 0;
		double *wavePtr = &aWavelength;

		if (!manager)
		{
			return;
		}

		manager->wavelengthForImage(noCryst, (void **)&wavePtr);

		if (aWavelength > 0)
		{
			this->wavelength = aWavelength;
		}

		logged << "Searching for wavelength in " << noCryst << " and found wavelength " << aWavelength << std::endl;
		sendLog();
	}
}

void MtzManager::loadReflections()
{
	if (reflectionCount())
	{
		return;
	}

	if (getFilename().length() == 0)
	{
		std::cerr
		<< "Cannot load reflections as no filename has been specified."
		<< std::endl;
		return;
	}

	if (!FileReader::exists(getFilename()))
	{
		logged << "Cannot find MTZ file for " << getFilename() << std::endl;
		sendLog();
		return;
	}

	bool recalculateWavelengths = FileParser::getKey("RECALCULATE_WAVELENGTHS", false);

	std::string fullFilename = (dropped ? "tmp-" : "") + getFilename();

	if (dropped)
	{
		fullFilename = FileReader::addOutputDirectory(fullFilename);
	}

	MTZ *mtz = MtzGet(fullFilename.c_str(), 0);

	int fromMtzNum = MtzSpacegroupNumber(mtz);

	int spgnum = FileParser::getKey("SPACE_GROUP", fromMtzNum);

	bool throwaway = FileParser::getKey("THROWAWAY_UNACCEPTED", false);

	setSpaceGroupNum(spgnum);

	if (getSpaceGroup() == NULL)
	{
		setRejected(true);
		return;
	}

	if (mtz == NULL)
		return;

	float *refldata = (float *) CCP4::ccp4_utils_malloc(
														(mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));
	memset(refldata, '\0',
		   (mtz->ncol_read + 1) * mtz->nref_filein * sizeof(float));

	float *adata = (float *) CCP4::ccp4_utils_malloc(
													 (mtz->ncol_read) * sizeof(float));

	MtzRrefl(mtz->filein, mtz->ncol_read * mtz->nref_filein, refldata);

	MTZCOL *col_f = MtzColLookup(mtz, "I");

	if (col_f == NULL)
	{
		col_f = MtzColLookup(mtz, "IMEAN");
	}

	if (col_f == NULL)
	{
		col_f = MtzColLookup(mtz, "F");

		if (col_f == NULL)
		{
			col_f = MtzColLookup(mtz, "FP");

			if (col_f == NULL)
				col_f = MtzColLookup(mtz, "FC");

			if (col_f != NULL)
			{
				std::cout << "Warning: using observed amplitude instead of intensity" << std::endl;
			}
			else
			{
				std::cout << "Warning: could not find any amplitude columns" << std::endl;
			}
		}
		else
		{
			std::cout << "Warning: using calculated amplitude instead of intensity" << std::endl;
		}
	}

	if (col_f == NULL)
	{
		std::cout << "Warning: intensity column not labelled I or IMEAN or amplitude" << std::endl;
		exit(1);
	}

	MTZCOL *col_csigf = MtzColLookup(mtz, "CSIGI");
	MTZCOL *col_sigf = MtzColLookup(mtz, "SIGI");

	if (col_sigf == NULL)
	{
		col_sigf = MtzColLookup(mtz, "SIGIMEAN");
	}

	if (col_sigf == NULL)
	{
		col_sigf = MtzColLookup(mtz, "SIGF");

		if (col_sigf != NULL)
		{
			std::cout << "Warning: using SIGF column" << std::endl;
		}
	}

	if (col_sigf == NULL)
	{
		std::cout << "Warning: sigma column not labelled SIGI or SIGIMEAN" << std::endl;
	}

	MTZCOL *col_wave = MtzColLookup(mtz, "WAVE");
	MTZCOL *col_partials = MtzColLookup(mtz, "PART");
	MTZCOL *col_phase = MtzColLookup(mtz, "PHIC");
	MTZCOL *col_shiftx = MtzColLookup(mtz, "SHIFTX");
	MTZCOL *col_shifty = MtzColLookup(mtz, "SHIFTY");
	MTZCOL *col_reject = MtzColLookup(mtz, "REJECT");
	MTZCOL *col_h = MtzColLookup(mtz, "H");
	MTZCOL *col_k = MtzColLookup(mtz, "K");
	MTZCOL *col_l = MtzColLookup(mtz, "L");

	MTZXTAL **xtals = MtzXtals(mtz);
	std::vector<float> cell;
	cell.resize(6);
	ccp4_lrcell(xtals[0], &cell[0]);

	bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);

	if (!fixUnitCell)
	{
		setUnitCell(cell);
	}
	else
	{
		vector<double> givenUnitCell = FileParser::getKey("UNIT_CELL", vector<double>());

		if (givenUnitCell.size() == 6)
		{
			if (!matrix)
			{
				matrix = MatrixPtr(new Matrix());
			}

			setUnitCell(givenUnitCell);
			lockUnitCellDimensions();
			getMatrix()->changeOrientationMatrixDimensions(givenUnitCell);
		}
		else
		{
			setUnitCell(cell);
		}
	}

	std::vector<double> unitCell = getUnitCell();

	for (int i = 0; i < mtz->nref_filein * mtz->ncol_read; i += mtz->ncol_read)
	{
		memcpy(adata, &refldata[i], mtz->ncol_read * sizeof(float));

		int h = adata[col_h->source - 1];
		int k = adata[col_k->source - 1];
		int l = adata[col_l->source - 1];

		float intensity = adata[col_f->source - 1];

		float countingSigma = 1;
		float sigma = 1;
		if (col_sigf != NULL)
		{
			sigma = adata[col_sigf->source - 1];
		}

		if (col_csigf != NULL)
		{
			countingSigma = adata[col_csigf->source - 1];
		}
		else if (col_csigf == NULL && col_sigf != NULL)
		{
			countingSigma = adata[col_sigf->source - 1];
		}

		float wavelength = 1;
		float partiality = 1;
		float phase = 0;
		float shiftX = 0;
		float shiftY = 0;
		int rejectFlags = 0;

		if (col_reject != NULL)
		{
			rejectFlags = adata[col_reject->source - 1];
		}

		if (col_wave != NULL && col_partials != NULL)
		{
			wavelength = adata[col_wave->source - 1];
			partiality = adata[col_partials->source - 1];
		}

		if (partiality < 0.05 && throwaway) {
			continue;
		}

		if (col_shiftx != NULL && col_shifty != NULL)
		{
			shiftX = adata[col_shiftx->source - 1];
			shiftY = adata[col_shifty->source - 1];
		}

		if (col_phase != NULL)
		{
			phase = adata[col_phase->source - 1];
		}

		if (col_partials == NULL && intensity != intensity)
			continue;

		MillerPtr miller = MillerPtr(new Miller(this, h, k, l));
		miller->setData(intensity, sigma, partiality, wavelength);
		miller->setCountingSigma(countingSigma);

		/*   if (dropped)
		 {
		 miller->setSigma(sigma);
		 }
		 */
		miller->setPhase(phase);
		miller->setCorrectedX(shiftX);
		miller->setCorrectedY(shiftY);
		miller->setRejected(rejectFlags);
		miller->matrix = this->matrix;
		miller->setScale(scale);
		addMiller(miller);

		if (miller->isSpecial())
		{
			logged << "Adding chosen Miller from " << getFilename() << " with intensity " << intensity << std::endl;
			sendLog();
		}
	}

	free(refldata);
	free(adata);

	std::ostringstream log;

	log << "Loaded " << mtz->nref_filein << " reflections (" << accepted()
	<< " accepted) for " << getFilename() << std::endl;

	bool lowMem = FileParser::getKey("LOW_MEMORY_MODE", false);

	Logger::mainLogger->addStream(&log, (lowMem ? LogLevelDetailed : LogLevelNormal));

	setSigmaToUnity();

	if (recalculateWavelengths && matrix)
		this->recalculateWavelengths();

	MtzFree(mtz);

	getWavelengthFromHDF5();

	fullyLoaded = true;

}

void MtzManager::setReference(MtzManager *reference)
{
	if (reference != NULL)
		Logger::mainLogger->addString("Setting reference to " + reference->getFilename());
	MtzManager::referenceManager = reference;
}


ReflectionPtr MtzManager::findReflectionWithId(ReflectionPtr exampleRefl, size_t *lowestId)
{
	if (reflectionCount() == 0)
	{
		if (lowestId) *lowestId = 0;
		return ReflectionPtr();
	}

	std::vector<ReflectionPtr>::iterator low;
	low = std::lower_bound(reflections.begin(), reflections.end(), exampleRefl, Reflection::reflLessThan);
	if (lowestId) *lowestId = (low - reflections.begin());

	if (low < reflections.end() && (*low)->getReflId() == exampleRefl->getReflId())
	{
		return *low;
	}
	else
	{
		return ReflectionPtr();
	}
}

int MtzManager::findReflectionWithId(long unsigned int refl_id, ReflectionPtr *reflection, bool insertionPoint)
{
	if (reflectionCount() == 0)
	{
		*reflection = ReflectionPtr();
		return -1;
	}

	int lower = 0;
	int higher = reflectionCount() - 1;
	int new_bound = (higher + lower) / 2;

	if ((refl_id < this->reflection(lower)->getReflId())
		|| (refl_id > this->reflection(higher)->getReflId()))
	{
		if (insertionPoint)
		{
			if (refl_id < this->reflection(lower)->getReflId())
				return 0;

			if (refl_id > this->reflection(higher)->getReflId())
				return reflectionCount();
		}

		*reflection = ReflectionPtr();

		return -1;
	}

	while (this->reflection(new_bound)->getReflId() != refl_id)
	{
		if (new_bound == higher || new_bound == lower)
		{
			if (insertionPoint)
			{
				int lowest = 0;

				int start = lower - 2 >= 0 ? lower - 2 : 0;

				for (int i = start; i < higher + 2 && i < reflectionCount(); i++)
				{
					if (this->reflection(i)->getReflId() < refl_id)
						lowest = i;
				}

				return lowest + 1;
			}

			if (this->reflection(higher)->getReflId() == refl_id)
			{
				(*reflection) = this->reflection(higher);
				return -1;
			}

			if (this->reflection(lower)->getReflId() == refl_id)
			{
				(*reflection) = this->reflection(lower);
				return -1;
			}

			*reflection = ReflectionPtr();
			return -1;
		}

		if (this->reflection(new_bound)->getReflId() > refl_id)
		{
			higher = new_bound;
		}
		else if (this->reflection(new_bound)->getReflId() < refl_id)
		{
			lower = new_bound;
		}

		new_bound = (higher + lower) / 2;
	}

	(*reflection) = this->reflection(new_bound);

	return -1;
}

void MtzManager::findCommonReflections(MtzManager *other,
									   vector<ReflectionPtr> &reflectionVector1, vector<ReflectionPtr> &reflectionVector2,
									   int *num, bool acceptableOnly, bool preserve)
{
	for (int i = 0; i < reflectionCount(); i++)
	{
		ReflectionPtr myRef = reflection(i);
		if (!myRef->acceptedCount() && acceptableOnly)
		{
			continue;
		}

		ReflectionPtr otherReflection = other->findReflectionWithId(myRef);

		if (otherReflection && otherReflection->millerCount() > 0)
		{
			reflectionVector1.push_back(myRef);
			reflectionVector2.push_back(otherReflection);
		}
	}

	if (num != NULL)
	{
		*num = (int)reflectionVector1.size();
	}
}

MtzPtr MtzManager::getDifferenceManager()
{
	if (differenceManager && differenceManager->isFullyLoaded())
	{
		return differenceManager;
	}

	tableMutex.lock();

	if (differenceManager)
	{
		tableMutex.unlock();
		return differenceManager;
	}

	differenceManager = MtzPtr(new MtzManager());
	differenceManager->setFilename(FileParser::getKey("DIFFERENCE_MTZ", std::string("")));
	differenceManager->loadReflections();

	if (differenceManager->reflectionCount() <= 0)
	{
		std::ostringstream logged;
		logged << "Difference MTZ required, please specify with DIFFERENCE_MTZ." << std::endl;
		staticLogAndExit(logged);
	}

	tableMutex.unlock();

	return differenceManager;
}

std::vector<double> MtzManager::getDifferencesWith(MtzPtr other, std::vector<double> &refls)
{
	std::vector<ReflectionPtr> myRefls, yourRefls;
	findCommonReflections(&*other, myRefls, yourRefls, NULL, true);

	std::vector<double> millerDiffs;

	for (int i = 0; i < myRefls.size(); i++)
	{
		double myInt = myRefls[i]->meanIntensity();
		double yourInt = yourRefls[i]->meanIntensity();

		double diff = yourInt - myInt;

		if (diff != diff)
		{
			continue;
		}

		MtzPtr diffMtz = getDifferenceManager();
		ReflectionPtr refRefl = ReflectionPtr();
		diffMtz->findReflectionWithId(myRefls[i]->getReflId(), &refRefl);

		if (!refRefl)
		{
			continue;
		}

		millerDiffs.push_back(diff);
		refls.push_back(refRefl->meanIntensity());
	}

	return millerDiffs;
}

void MtzManager::bFactorAndScale(double *scale, double *bFactor, double exponent)
{
	applyScaleFactor(1, 0, 0, true);
	applyBFactor(0);

	MtzManager *reference = MtzManager::getReferenceManager();

	vector<ReflectionPtr> refReflections, imageReflections;
	this->findCommonReflections(reference, imageReflections, refReflections, NULL);

	vector<boost::tuple<double, double, double> > pointsToFit;

	for (int i = 0; i < imageReflections.size(); i++)
	{
		ReflectionPtr imgRef = imageReflections[i];
		ReflectionPtr refRef = refReflections[i];

		if (!imgRef->anyAccepted())
			continue;

		double refMean = refRef->meanIntensity();
		double imgMean = imgRef->meanIntensity();
		double imgWeight = imgRef->meanPartiality();
		imgWeight *= log(1 / imgRef->getResolution());

		double ratio = refMean / imgMean;
		ratio = 1 / ratio;

		if (ratio != ratio)
			continue;

		double resolution = 1 / imgRef->getResolution();

		double four_d_squared = 4 * pow(resolution, 2);

		double right_exp = 1 / four_d_squared;

		double intensityRatio = (imgMean / refMean);
		double logIntensityRatio = log(intensityRatio);

		if (logIntensityRatio != logIntensityRatio)
			continue;

		boost::tuple<double, double, double> point = boost::make_tuple(right_exp,
																	   logIntensityRatio, imgWeight);

		pointsToFit.push_back(point);
	}

	double gradient = 0;
	double intercept = 0;
    
 
    if (getFilename() == "img-lyso_58_0.mtz" || getFilename() == "img-lyso_83_0.mtz" || getFilename() == "img-lyso_91_0.mtz" || getFilename() == "img-lyso_1_0.mtz")
    {
        {
            CSVPtr csv = CSVPtr(new CSV(3, "x", "y", "weight"));
            for (int i=0; i < pointsToFit.size(); i++)
            {
                double x = boost::get<0>(pointsToFit[i]);
                double y = boost::get<1>(pointsToFit[i]);
                double weight = boost::get<2>(pointsToFit[i]);
                std::cout << x << "," << y << "," << weight << std::endl;
                csv->addEntry(3, (double)x, (double)y, (double)weight);
            }
            csv->writeToFile("pointsToFit" + getFilename() +".csv");
        }
    }
 


	regression_line(pointsToFit, intercept, gradient);

	double k = 1 / exp(intercept);

	double b = gradient / -2;//pow(fabs(gradient), 1 / (double) exponent);

	if (b != b)
		b = 0;

	logged << "Last scale for " << getFilename() << ": " << k << ", bFactor: " << b << std::endl;
	sendLog(LogLevelDetailed);

	*scale = k;
	*bFactor = b;

	applyScaleFactor(k);
	applyBFactor(b);
}

void MtzManager::scaleToMtz(MtzManager *otherManager, bool info, double lowRes, double highRes)
{
	vector<ReflectionPtr> reflections1;
	vector<ReflectionPtr> reflections2;

	int count = 0;
	int num = 0;

	this->findCommonReflections(otherManager, reflections1, reflections2, &num, true, true);

	if (num <= 1)
		return;

	double x_squared = 0;
	double x_y = 0;

	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < reflections1[i]->millerCount(); j++)
		{
			MillerPtr miller = reflections1[i]->miller(j);

			if (miller->isFree())
				continue;

			if (!miller->accepted())
				continue;

			if (!reflections1[i]->betweenResolutions(lowRes, highRes))
				continue;

			double int1 = miller->intensity();
			double int2 = reflections2[i]->meanIntensity();
			double weight = miller->getPartiality();

			if ((int1 != int1) || (int2 != int2) || (weight != weight))
				continue;

			x_y += int1 * int2 * weight;
			x_squared += int2 * int2 * weight;

			count++;
		}
	}

	double grad = (x_squared / x_y);

	if (grad < 0)
		grad = -1;

	applyScaleFactor(grad);
}





void MtzManager::applyBFactor(double bFactor)
{
	if (bFactor == 0)
	{
		return;
	}

	this->bFactor = bFactor;

	for (int i = 0; i < reflections.size(); i++)
	{
		for (int j = 0; j < reflections[i]->millerCount(); j++)
		{
			reflection(i)->miller(j)->setBFactor(bFactor);
		}
	}
}

void MtzManager::applyScaleFactor(double scaleFactor,
								  double lowRes, double highRes, bool absolute)
{
	if (scaleFactor == scaleFactor)
	{
		if (absolute)
			scale = scaleFactor;
		else
			scale *= scaleFactor;
	}

	logged << "Applying scale factor to " << getFilename() << " - " << scaleFactor << " (now " << this->scale << ")" << std::endl;

	for (int i = 0; i < reflections.size(); i++)
	{
		for (int j = 0; j < reflections[i]->millerCount(); j++)
		{
			if (reflection(i)->betweenResolutions(lowRes, highRes))
			{
				if (absolute)
					reflection(i)->miller(j)->setScale(scale);
				else
					reflection(i)->miller(j)->applyScaleFactor(scaleFactor);
			}
		}
	}

	sendLog(LogLevelDebug);

}

double MtzManager::averageIntensity(void)
{
	double total_intensity = 0;
	double total = 0;

	for (int i = 0; i < reflections.size(); i++)
	{
		double intensity = reflections[i]->meanIntensity();
		if (intensity == intensity)
		{
			double weight = reflections[i]->meanPartiality();
			total_intensity += intensity * weight;
			total += weight;
		}
	}

	total_intensity /= total;
	return total_intensity;
}

void MtzManager::applyPolarisation(void)
{
	for (int i = 0; i < reflections.size(); i++)
	{
		for (int j = 0; j < reflections[i]->millerCount(); j++)
		{
			reflections[i]->miller(j)->applyPolarisation(1.459);
		}
	}
}


void MtzManager::writeToFile(std::string newFilename, bool announce, bool plusAmbiguity,
							 bool withScale)
{
	int columns = 11;

	float cell[6], wavelength;
	float *fdata = new float[columns];
	std::vector<double> unitCell = getUnitCell();

	/* variables for symmetry */
	CCP4SPG *mtzspg = getSpaceGroup();
	float rsm[192][4][4];
	char ltypex[2];

	/* variables for MTZ data structure */
	MTZ *mtzout;
	MTZXTAL *xtal;
	MTZSET *set;
	MTZCOL *colout[11];

	/*  Removed: General CCP4 initializations e.g. HKLOUT on command line */

	cell[0] = unitCell[0];
	cell[1] = unitCell[1];
	cell[2] = unitCell[2];
	cell[3] = unitCell[3];
	cell[4] = unitCell[4];
	cell[5] = unitCell[5];
	wavelength = this->getWavelength();

	if (newFilename == "")
	{
		newFilename = getFilename();
	}

	std::string outputFile = FileReader::addOutputDirectory(newFilename);

	mtzout = MtzMalloc(0, 0);
	ccp4_lwtitl(mtzout, "Written from Helen's XFEL tasks ", 0);
	mtzout->refs_in_memory = 0;
	mtzout->fileout = MtzOpenForWrite(outputFile.c_str());

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
	colout[3] = MtzAddColumn(mtzout, set, "I", "J");
	colout[4] = MtzAddColumn(mtzout, set, "SIGI", "Q");
	colout[5] = MtzAddColumn(mtzout, set, "PART", "R");
	colout[6] = MtzAddColumn(mtzout, set, "WAVE", "R");
	colout[7] = MtzAddColumn(mtzout, set, "SHIFTX", "R");
	colout[8] = MtzAddColumn(mtzout, set, "SHIFTY", "R");
	colout[9] = MtzAddColumn(mtzout, set, "CSIGI", "R");
	colout[10] = MtzAddColumn(mtzout, set, "REJECT", "R");

	int num = 0;

	if (plusAmbiguity)
	{
		flipToActiveAmbiguity();
	}

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);

			double intensity = reflections[i]->miller(j)->getRawestIntensity();

			if (withScale)
			{
				intensity = reflections[i]->miller(j)->intensity();
			}
			double sigma = reflections[i]->miller(j)->getRawSigma();
			double partiality = reflections[i]->miller(j)->getPartiality();
			double bFactor = reflections[i]->miller(j)->getBFactorScale();
			double countingSigma = reflections[i]->miller(j)->getRawCountingSigma();
			double rejectFlags = reflections[i]->miller(j)->getRejectionFlags();

			if (intensity != intensity)
			{
				continue;
			}

			num++;

			int h = reflections[i]->miller(j)->getH();
			int k = reflections[i]->miller(j)->getK();
			int l = reflections[i]->miller(j)->getL();

            
			fdata[0] = h;
			fdata[1] = k;
			fdata[2] = l;
            fdata[3] = reflections[i]->miller(j)->getRawestIntensity(); //intensity/bfactor
            fdata[4] = sigma;
			fdata[5] = partiality;
			fdata[6] = reflections[i]->miller(j)->getWavelength();
			fdata[7] = reflections[i]->miller(j)->getCorrectedX();
			fdata[8] = reflections[i]->miller(j)->getCorrectedY();
			fdata[9] = countingSigma;
			fdata[10] = rejectFlags;

			ccp4_lwrefl(mtzout, fdata, colout, columns, num);
		}
	}

	if (plusAmbiguity)
	{
		resetFlip();
	}


	MtzPut(mtzout, " ");
	MtzFree(mtzout);

	LogLevel shouldAnnounce = announce ? LogLevelNormal : LogLevelDebug;

	std::ostringstream logged;
	logged << "Written " << num << " refls to file " << newFilename << std::endl;
	Logger::mainLogger->addStream(&logged, shouldAnnounce);

	delete [] fdata;
}


MtzManager::~MtzManager(void)
{
	nearbyMillers.clear();
	vector<MillerPtr>().swap(nearbyMillers);

	clearReflections();
}

void MtzManager::description(void)
{
	logged << "Filename: " << getFilename() << std::endl;
	logged << "Number of reflections: " << reflectionCount() << std::endl;
	logged << "Number of accepted Millers: " << accepted() << std::endl;
	logged << "Average intensity: " << this->averageIntensity() << std::endl;

	sendLog();
}

void MtzManager::setSigmaToUnity()
{
	for (int i = 0; i < reflections.size(); i++)
	{
		for (int j = 0; j < reflections[i]->millerCount(); j++)
		{
			reflections[i]->miller(j)->setSigma(1);
		}
	}
}

int MtzManager::accepted()
{
	int acceptedCount = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			if (reflection(i)->miller(j)->accepted())
				acceptedCount++;
		}
	}

	return acceptedCount;
}

int MtzManager::rejectOverlaps()
{
	int count = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		count += reflection(i)->checkOverlaps();
	}

	return count;
}

int MtzManager::removeStrongSpots(std::vector<SpotPtr> *spots, bool actuallyDelete)
{
	int count = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		ReflectionPtr ref = reflection(i);

		count += ref->checkSpotOverlaps(spots, actuallyDelete);
	}

	return count;
}

int MtzManager::ambiguityCount()
{
	int count = reflection(0)->ambiguityCount();

	return count;
}

void MtzManager::incrementActiveAmbiguity()
{
	int count = ambiguityCount();

	activeAmbiguity++;
	activeAmbiguity = (activeAmbiguity % count);

	for (int i = 0; i < reflectionCount(); i++)
	{
		reflection(i)->setActiveAmbiguity(activeAmbiguity);
	}
}

void MtzManager::flipToActiveAmbiguity()
{
	for (int i = 0; i < reflectionCount(); i++)
	{
		reflection(i)->setFlip(activeAmbiguity);
	}

	this->sortReflections();
}


void MtzManager::resetFlip()
{
	for (int i = 0; i < reflectionCount(); i++)
	{
		reflection(i)->resetFlip();
	}

	this->sortReflections();
}

void MtzManager::setActiveAmbiguity(int newAmbiguity)
{
	activeAmbiguity = newAmbiguity;

	for (int i = 0; i < reflectionCount(); i++)
	{
		reflection(i)->setActiveAmbiguity(newAmbiguity);
	}
}

double MtzManager::maxResolution()
{
	double maxResolution = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		if (reflection(i)->getResolution() > maxResolution)
		{
			if (std::isfinite(maxResolution))
			{
				maxResolution = reflection(i)->getResolution();
			}
		}
	}

	return maxResolution;
}

std::string MtzManager::getParamLine()
{
	std::vector<double> unitCell = getUnitCell();
	std::ostringstream line;
	line << "params ";

	line << mosaicity << " " << spotSize << " " << wavelength << " ";
	line << bFactor << " " << scale << " ";
	line << unitCell[0] << " " << unitCell[1] << " " << unitCell[2];

	return line.str();
}

void MtzManager::setParamLine(std::string line)
{
	std::vector<std::string> components = FileReader::split(line, ' ');

	if (components.size() < 8)
	{
		return;
	}

	std::vector<double> unitCell = getUnitCell();

	mosaicity = atof(components[1].c_str());
	spotSize = atof(components[2].c_str());
	wavelength = atof(components[3].c_str());
	bFactor = atof(components[4].c_str());
	scale = atof(components[5].c_str());
	unitCell[0] = atof(components[6].c_str());
	unitCell[1] = atof(components[7].c_str());
	unitCell[2] = atof(components[8].c_str());

	setInitialValues = true;
}

void MtzManager::recalculateWavelengths()
{
	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			reflection(i)->miller(j)->getWavelength();
		}
	}
}

std::string MtzManager::parameterHeaders()
{
	std::ostringstream summary;

	summary << "Filename,Correl,Rsplit,Partcorrel,";
	summary << "Refcount,Strong,Mosaicity,Wavelength,Bandwidth,";

	summary << "hRot,kRot,lRot,";

	summary << "rlpSize,exp,cellA,cellB,cellC,scale,bfactor,avInt";

	return summary.str();
}

std::string MtzManager::writeParameterSummary()
{
	std::vector<double> unitCell = getUnitCell();

	std::ostringstream summary;
	summary << getFilename() << ","
	<< getRefCorrelation() << ","
	<< getLastRSplit() << ","
	<< getRefPartCorrel() << ","
	<< accepted() << ","
	<< getTotalReflections() << ","
	<< getMosaicity() << ","
	<< getWavelength() << ","
	<< getBandwidth() << ","
	<< hRot << "," << kRot << "," << lRot << ","
	<< getSpotSize() << ","
	<< getExponent() << ","
	<< unitCell[0] << ","
	<< unitCell[1] << ","
	<< unitCell[2] << ","
	<< getScale()<< ","
    << getBFactor() << ","
    << averageIntensity() ;
	return summary.str();
}



void MtzManager::millersToDetector()
{
	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);

			DetectorPtr det = miller->getDetector();

			if (!det)
			{
				miller->positionOnDetector(NULL, NULL, false);
				det = miller->getDetector();

				if (!det)
				{
					continue;
				}
			}

			det->addMillerCarefully(miller);
		}
	}
}

// MARK: IOMRefiner


bool MtzManager::checkMaxResolution()
{
    bool findMaxResolution = false;
    double maxIntegratedResolution = FileParser::getKey("MAX_INTEGRATED_RESOLUTION",0.0);
    double maxRefinedResolution = FileParser::getKey("MAX_REFINED_RESOLUTION",MAX_OPTIMISATION_RESOLUTION);
    if (maxRefinedResolution < -0.01)
    {
        findMaxResolution = true;
    }
    return findMaxResolution;
}


void MtzManager::calculateNearbyMillers()
{
	double wavelength = getWavelength();

	double minBandwidth = wavelength * (1 - testBandwidth * 2);
	double maxBandwidth = wavelength * (1 + testBandwidth * 2);

	double sphereThickness = FileParser::getKey("SPHERE_THICKNESS", 0.01);

	double minSphere = 1 / wavelength * (1 - sphereThickness);
	double maxSphere = 1 / wavelength * (1 + sphereThickness);

	double minSphereSquared = pow(minSphere, 2);
	double maxSphereSquared = pow(maxSphere, 2);

    bool checkRes = checkMaxResolution();
    logged << "checkRes: " << checkRes << std::endl;
    sendLog();
    
    maxResolutionAll = FileParser::getKey("MAX_REFINED_RESOLUTION",MAX_OPTIMISATION_RESOLUTION);
    
    if (checkRes)
    {
        maxResolutionAll = getImagePtr()->getFoundMaxResolution();
    }
    
    
    //bool checkpointer = !getImagePtr();
    
    
    logged << "Firstly, setting maximum resolution to found resolution " << maxResolutionAll << std::endl;
    std::cout << maxResolutionAll << std::endl;
    logged << "HELLOS1! Setting maximum resolution to resolution " << maxResolutionAll << std::endl;
    sendLog();
    
	nearbyMillers.clear();

	int maxMillers[3];

	Miller::rotateMatrixHKL(hRot, kRot, lRot, matrix, &rotatedMatrix);

	rotatedMatrix->maxMillers(maxMillers, maxResolutionAll);

	logged << "Integrating to maximum Miller indices: (" << maxMillers[0]
	<< ", " << maxMillers[1] << ", " << maxMillers[2] << ") to resolution " << maxResolutionAll << " Å." << std::endl;

	int overRes = 0;
	int underRes = 0;
	double maxDistSquared = pow(1 / maxResolutionAll, 2);
	double minResolutionSquared = pow(1 / minResolutionAll, 2);

	CCP4SPG *spg = getSpaceGroup();
	if (!spg)
	{
		logged << "Error: Trying to integrate reflections, but could not load missing SPACE_GROUP information." << std::endl;
		sendLogAndExit();
	}

	for (int h = -maxMillers[0]; h < maxMillers[0]; h++)
	{
		for (int k = -maxMillers[1]; k < maxMillers[1]; k++)
		{
			for (int l = -maxMillers[2]; l < maxMillers[2]; l++)
			{
				if (ccp4spg_is_sysabs(getSpaceGroup(), h, k, l))
					continue;

				vec hkl = new_vector(h, k, l);
				rotatedMatrix->multiplyVector(&hkl);

				if (hkl.l > 0.01)
					continue;

				vec beam = new_vector(0, 0, -1 / wavelength);
				vec beamToMiller = vector_between_vectors(beam, hkl);

				double sphereRadiusSquared = length_of_vector_squared(beamToMiller);

				if (sphereRadiusSquared < minSphereSquared ||
					sphereRadiusSquared > maxSphereSquared)
				{
					double ewaldSphere = getEwaldSphereNoMatrix(hkl);

					if (ewaldSphere < minBandwidth || ewaldSphere > maxBandwidth)
					{
						continue;
					}
				}

				double res = length_of_vector_squared(hkl);

				if (res > maxDistSquared)
				{
					overRes++;
					continue;
				}

				if (minResolutionAll > 0 && res < minResolutionSquared)
				{
					underRes++;
					continue;
				}

				if (h == 0 && k == 0 && l == 0)
					continue;

				MillerPtr newMiller = MillerPtr(new Miller(NULL, h, k, l));
				newMiller->setMatrix(rotatedMatrix);
				nearbyMillers.push_back(newMiller);
			}
		}
	}

	logged << "Rejected " << overRes << " due to being over resolution edge of "
	<< maxResolutionAll << " Å." << std::endl;

	if (minResolutionAll > 0)
		logged << "Rejected " << underRes <<
		" due to being under resolution edge of " << minResolutionAll << " Å." << std::endl;
    

	sendLog();
}

void MtzManager::checkNearbyMillers()
{
	MatrixPtr matrix = getMatrix();
    maxResolutionAll = FileParser::getKey("MAX_REFINED_RESOLUTION",MAX_OPTIMISATION_RESOLUTION);
    
    bool checkRes = checkMaxResolution();
    if (checkRes)
    {
        maxResolutionAll = getImagePtr()->getFoundMaxResolution();
    }

	std::vector<double> unitCell = getUnitCell();

	if (!(unitCell.size() == 6))
	{
		logged << "Error: Trying to integrate reflections, but need the UNIT_CELL dimension information." << std::endl;
		sendLogAndExit();
	}


	lockUnitCellDimensions();
	matrix->changeOrientationMatrixDimensions(unitCell);

	clearReflections();

	double wavelength = getWavelength();

	double maxD = 1 / maxResolutionAll;
    
	if (maxResolutionAll == 0)
		maxD = FLT_MAX;
    
	logged << "(" << getBasename() << ") testing " << nearbyMillers.size()
	<< " reflections close to the Ewald sphere with wavelength "
	<< wavelength << std::endl;

	int cutResolution = 0;
	int partialityTooLow = 0;
	int unacceptableIntensity = 0;

	logged << "Rotating by " << hRot << ", " << kRot << ", " << lRot << std::endl;
	sendLog(LogLevelDetailed);

	Miller::rotateMatrixHKL(hRot, kRot, lRot, matrix, &rotatedMatrix);

	for (int i = 0; i < nearbyMillers.size(); i++)
	{
		MillerPtr miller = nearbyMillers[i];
		miller->setMatrix(rotatedMatrix);
		miller->setImage(getImagePtr());

		double d = miller->resolution();

		if (d > maxD)
		{
            //Does this mean that the res is not cut?!
			cutResolution++;
			continue;
		}

		miller->crossesBeamRoughly(rotatedMatrix, 0.0, testSpotSize * 2, wavelength, testBandwidth);

		if (i == 0)
		{
			logged << "Calculated partiality with parameters: hRot " << hRot << ", kRot " << kRot << ", spot size " << testSpotSize << ", wavelength " << wavelength << ", bandwidth " << bandwidth << std::endl;
			sendLog(LogLevelDebug);
		}

		if (miller->getPartiality() <= 0.05)
		{
			partialityTooLow++;
			continue;
		}

		miller->recalculateWavelength();
		miller->setImage(getImagePtr());
		miller->setMtzParent(this);

		miller->positionOnDetector(NULL, NULL, false);

		if (!miller->hasDetector())
		{
			continue;
		}

		addMiller(miller);
	}

	logged << "Using wavelength " << std::setprecision(9) << getImagePtr()->getWavelength() << " Å." << std::endl;
	logged << "Beyond resolution cutoff: " << cutResolution << std::endl;
	logged << "Partiality equal to 0: " << partialityTooLow << std::endl;
	logged << "Reflections left: " << millerCount() << std::endl;

	sendLog(LogLevelDetailed);
}

void MtzManager::integrateChosenMillers(bool quick)
{
	std::vector<double> unitCell = this->getUnitCell();
	MatrixPtr newMat = Matrix::matrixFromUnitCell(unitCell);
	matrix->setUnitCell(newMat);

	Miller::rotateMatrixHKL(hRot, kRot, lRot, matrix, &rotatedMatrix);

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);
			miller->setMatrix(rotatedMatrix);
			miller->integrateIntensity(quick);
			miller->recalculateWavelength();
			miller->setMatrix(matrix);

			if (!miller->accepted())
			{
				reflection(i)->removeMiller(j);
				j--;
			}
		}

		if (reflection(i)->millerCount() == 0)
		{
			removeReflection(i);
			i--;
		}
	}
}

bool MtzManager::millerWithinBandwidth(MillerPtr miller)
{
	double minBandwidth = getWavelength() * (1 - testBandwidth * 2);
	double maxBandwidth = getWavelength() * (1 + testBandwidth * 2);

	if (FileParser::hasKey("WAVELENGTH_RANGE"))
	{
		std::vector<double> range = FileParser::getKey("WAVELENGTH_RANGE", std::vector<double>());

		if (range.size() >= 2)
		{
			minBandwidth = range[0];
			maxBandwidth = range[1];
		}
	}

	double wavelength = miller->getWavelength();

	return (wavelength > minBandwidth && wavelength < maxBandwidth);
}

double MtzManager::getReflectionWavelengthStdev()
{
	std::vector<double> values;

	bool considerDistances = FileParser::getKey("REFINE_DISTANCE_TO_EWALD", false);

	bool unbalanced = FileParser::getKey("UNBALANCED_REFLECTIONS", false);
	double mean = getWavelength();
	double sum = 0;
	double count = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);
			if (miller->reachesThreshold())
			{
				double wavelength = miller->getWavelength();

				if (wavelength != wavelength)
				{
					continue;
				}

				if (!considerDistances)
				{
					values.push_back(miller->getWavelength());
				}
				else
				{
					miller->setMatrix(rotatedMatrix);
					double val = miller->distanceToEwald(getWavelength());
					miller->setMatrix(getMatrix());
					sum += val;
					count++;
				}
			}
		}
	}

	if (considerDistances)
	{
		return sum / count;
	}

	double stdev = 0;

	if (unbalanced)
	{
		double mean = getImagePtr()->getWavelength();
		stdev = bitty_deviation(&values, NULL, mean);
	}
	else
	{
		stdev = bitty_deviation(&values, NULL);
	}

	return stdev;
}

int MtzManager::getTotalReflections()
{
	int count = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);

			if (miller->reachesThreshold() && millerWithinBandwidth(miller))
			{
				count++;
			}
		}
	}

	return count;
}


double MtzManager::hkScoreFinalise(void *object)
{
	MtzManager *mtz = static_cast<MtzManager *>(object);
	mtz->integrateChosenMillers(true);
	return mtz->hkScore(true);
}

double MtzManager::hkScoreStdevWrapper(void *object)
{
	MtzManager *mtz = static_cast<MtzManager *>(object);
	mtz->integrateChosenMillers(true);
	return mtz->hkScore(false);
}

double MtzManager::lScoreWrapper(void *object)
{
	MtzManager *mtz = static_cast<MtzManager *>(object);
	mtz->integrateChosenMillers();
	return mtz->lScore(true);
}

double MtzManager::hkScore(bool finalise)
{
	double stdev = getReflectionWavelengthStdev();
	double refTotal = getTotalReflections();

	if (Logger::getPriorityLevel() >= LogLevelDetailed)
	{
		getWavelengthHistogram();
	}

	return stdev;

	if (finalise)
	{
		lastTotal = refTotal;
		lastStdev = stdev;
	}

	double totalWeight = 20;
	double stdevWeight = 15;

	double totalChange = refTotal / lastTotal - 1;
	double totalStdev = stdev / lastStdev - 1;

	if (lastTotal == 0 || lastStdev == 0)
	{
		return 1;
	}

	double score = (totalStdev * stdevWeight - totalChange * totalWeight);

	if (!finalise)
	{
		logged << "Score is " << score << " from " << refTotal << " reflections" << std::endl;
		sendLog(LogLevelDetailed);
	}

	return score;
}

double MtzManager::lScore(bool silent)
{
	double averageShift = 0;
	int count = 0;

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);

			if (miller->reachesThreshold())
			{
				std::pair<double, double> shift = miller->getShift();
				double shiftDistance = pow(shift.first, 2) + pow(shift.second, 2);

				averageShift -= Detector::lookupCache(shiftDistance);
				count++;
			}
		}
	}

	averageShift /= count;

	logged << "Average shift: " << averageShift << std::endl;
	sendLog(LogLevelDetailed);

	return averageShift;
}

void MtzManager::dropMillers()
{
	nearbyMillers.clear();
	vector<MillerPtr>().swap(nearbyMillers);
}


void MtzManager::getWavelengthHistogram(vector<double> *wavelengths,
										vector<int> *frequencies, bool silent)
{
	bool populateVectors = (wavelengths != NULL && frequencies != NULL);

	if (populateVectors)
	{
		wavelengths->clear();
		frequencies->clear();
	}

	if (!silent)
	{
		logged << "Wavelength histogram for " << this->getFilename() << std::endl;
		sendLog();
	}

	double wavelength = getWavelength();
	vector<double> wavelengthRange = FileParser::getKey("WAVELENGTH_RANGE", vector<double>(2, 0));

	double spread = testBandwidth * 2;
	double interval = (wavelength * spread * 2) / 20;
	double total = 0;

	double minLength = wavelength * (1 - spread);
	double maxLength = wavelength * (1 + spread);

	if (wavelengthRange[0] != 0 || wavelengthRange[1] != 0)
	{
		minLength = wavelengthRange[0];
		maxLength = wavelengthRange[1];
		interval = (maxLength - minLength) / 20;
	}

	for (double i = minLength; i < maxLength; i += interval)
	{
		if (populateVectors)
		{
			wavelengths->push_back(i);
		}

		double frequency = 0;
		double totalInInterval = 0;

		for (int j = 0; j < reflectionCount(); j++)
		{
			for (int k = 0; k < reflection(j)->millerCount(); k++)
			{
				MillerPtr miller = reflection(j)->miller(k);

				double ewald = miller->getWavelength();

				if (ewald < i || ewald > i + interval)
					continue;

				bool strong = miller->reachesThreshold();

				double weight = 1;

				totalInInterval += weight;

				if (strong)
				{
					frequency += weight;
				}

				if (frequency < 0)
				{
					frequency = 0;
				}
			}
		}

		if (!silent)
		{
			logged << std::setprecision(4) << i << "\t";

			for (int i=0; i < frequency; i++)
			{
				logged << ".";
			}

			logged << std::endl;
		}

		if (populateVectors)
		{
			frequencies->push_back(frequency);
		}
		total += frequency;
	}

	if (!silent)
	{
		logged << std::endl;
		logged << "Total strong reflections: " << total << std::endl;
		sendLog();
	}
}


void MtzManager::refineOrientationMatrix(bool force)
{
	bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);

	if (fixUnitCell)
	{
		vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
		setUnitCell(unitCell);
	}

	this->calculateNearbyMillers();
	checkNearbyMillers();

	lastStdev = getReflectionWavelengthStdev();
	lastTotal = getTotalReflections();

	int oldSearchSize = searchSize;
	int bigSize = FileParser::getKey("METROLOGY_SEARCH_SIZE_BIG", 3);
	bool refineOffset = FileParser::getKey("REFINE_EACH_DETECTOR_DISTANCE", false);

	bool refineUnitCellA = FileParser::getKey("OPTIMISING_UNIT_CELL_A", false);
	bool refineUnitCellB = FileParser::getKey("OPTIMISING_UNIT_CELL_B", false);
	bool refineUnitCellC = FileParser::getKey("OPTIMISING_UNIT_CELL_C", false);

	double stepSizeUnitCellA = FileParser::getKey("STEP_SIZE_UNIT_CELL_A", 0.2);
	double stepSizeUnitCellB = FileParser::getKey("STEP_SIZE_UNIT_CELL_B", 0.2);
	double stepSizeUnitCellC = FileParser::getKey("STEP_SIZE_UNIT_CELL_C", 0.2);

	bool refineCell = refineUnitCellA || refineUnitCellB || refineUnitCellC;

	// FIXME: read support for unit cell dimensions

	searchSize = oldSearchSize + bigSize;

	if (refineOrientations || force)
	{
		RefinementStrategyPtr lStrategy = RefinementStepSearchPtr(new RefinementStepSearch());
		lStrategy->setVerbose(Logger::getPriorityLevel() >= LogLevelDetailed);
		lStrategy->setEvaluationFunction(lScoreWrapper, this);
		lStrategy->setJobName("Refining in-detector-plane angle for " + getFilename());
		lStrategy->addParameter(this, getLRot, setLRot, initialStep, orientationTolerance, "lRot");

		if (refineOffset)
		{
			lStrategy->setJobName("Refining in-detector-plane angle & distance offset for " + getFilename());
			lStrategy->addParameter(&*this->getImagePtr(), Image::getDistanceOffset, Image::setDistanceOffset, 1.0, 0.001, "distance_offset");
		}

		lStrategy->refine();

		if (refineOffset)
		{
			logged << "Detector distance offset for " << getBasename() << " is "
			<< Image::getDistanceOffset(&*getImagePtr()) << " pixels." << std::endl;
			sendLog();
		}

		searchSize = oldSearchSize;

		calculateNearbyMillers();
		checkNearbyMillers();

		hkScoreFinalise(this);
		RefinementStrategyPtr hkStrategy = RefinementStrategy::userChosenStrategy();
		hkStrategy->setVerbose(Logger::getPriorityLevel() >= LogLevelDetailed);
		hkStrategy->setEvaluationFunction(hkScoreStdevWrapper, this);
		hkStrategy->setJobName("Refining angles for " + getFilename());
		hkStrategy->addParameter(this, getHRot, setHRot, initialStep, orientationTolerance, "hRot");
		hkStrategy->addCoupledParameter(this, getKRot, setKRot, initialStep, orientationTolerance, "kRot");


		if (refineCell)
		{
			if (refineUnitCellA)
			{
				hkStrategy->addParameter(this, getUnitCellAStatic, setUnitCellAStatic, stepSizeUnitCellA, 0.001, "unit_cell_a");
			}

			if (refineUnitCellB)
			{
				hkStrategy->addParameter(this, getUnitCellBStatic, setUnitCellBStatic, stepSizeUnitCellB, 0.001, "unit_cell_b");
			}

			if (refineUnitCellC)
			{
				hkStrategy->addParameter(this, getUnitCellCStatic, setUnitCellCStatic, stepSizeUnitCellC, 0.001, "unit_cell_c");
			}
		}

		hkStrategy->refine();
	}

	calculateNearbyMillers();
	checkNearbyMillers();
	integrateChosenMillers();
}


bool MtzManager::isGoodSolution()
{
	bool good = false;
	double goodSolutionSumRatio = FileParser::getKey("GOOD_SOLUTION_SUM_RATIO", 6.5);
	int goodSolutionHighestPeak = FileParser::getKey("GOOD_SOLUTION_HIGHEST_PEAK", 17);
	int badSolutionHighestPeak = FileParser::getKey("BAD_SOLUTION_HIGHEST_PEAK", 10);
	int minimumReflections = FileParser::getKey("MINIMUM_REFLECTION_CUTOFF", 30);
	std::ostringstream details;

	logged << "(" << getFilename() << ") Standard deviation: " << lastStdev << std::endl;
	sendLog(LogLevelNormal);

	vector<double> wavelengths;
	vector<int> frequencies;

	getWavelengthHistogram(&wavelengths, &frequencies, true);

	double totalMean = 0;
	double totalStdev = 0;
	int maxFrequency = 0;

	histogram_gaussian(&wavelengths, &frequencies, totalMean, totalStdev);

	for (int i = 0; i < frequencies.size(); i++)
	{
		if (frequencies[i] > maxFrequency)
		{
			maxFrequency = frequencies[i];
		}
	}

	std::sort(frequencies.begin(), frequencies.end(), std::greater<int>());

	double highSum = 0;
	std::vector<double> lowOnly;
	int cutoff = 3;

	for (int i = 0; i < frequencies.size(); i++)
	{
		if (i < cutoff)
		{
			highSum += frequencies[i];
			continue;
		}

		lowOnly.push_back(frequencies[i]);
	}

	highSum /= cutoff;

	double stdevLow = standard_deviation(&lowOnly, NULL, 0);

	if (highSum > stdevLow * goodSolutionSumRatio)
	{
		good = true;
		details << "(" << getFilename() << ") Sum ratio is sufficiently high (" << highSum << " vs " << stdevLow << ")" << std::endl;
	}

	if (frequencies.size() && frequencies[0] > goodSolutionHighestPeak)
	{
		details << "(" << getFilename() << ") Highest peak is high enough (" << frequencies[0] << " vs " << goodSolutionHighestPeak << ")" << std::endl;
		good = true;
	}

	if (frequencies.size() && frequencies[0] < badSolutionHighestPeak)
	{
		details << "(" << getFilename() << ") Highest peak is too low (" << frequencies[0] << " vs " << badSolutionHighestPeak << ")" << std::endl;
		good = false;
	}

	if (getTotalReflections() < minimumReflections)
	{
		details << "(" << getFilename() << ") However, not enough reflections (" << getTotalReflections() << " vs " << minimumReflections << ")" << std::endl;
		good = false;
	}

	Logger::mainLogger->addStream(&details, LogLevelNormal);

	return good;
}


void MtzManager::fakeSpots()
{
	this->calculateNearbyMillers();
	this->checkNearbyMillers();
	double partialCutoff = FileParser::getKey("PARTIAL_CUTOFF", 0.2);

	for (int i = 0; i < reflectionCount(); i++)
	{
		for (int j = 0; j < reflection(i)->millerCount(); j++)
		{
			MillerPtr miller = reflection(i)->miller(j);
			miller->recalculatePartiality(rotatedMatrix, mosaicity, spotSize,
										  getWavelength(), bandwidth, exponent);

			if (miller->getPartiality() > partialCutoff)
			{
				double x = 0;
				double y = 0;
				miller->positionOnDetector(&x, &y);
				vec intersection;
				double xSpot = 0;
				double ySpot = 0;

				DetectorPtr myDet;
				myDet = Detector::getMaster()->intersectionForMiller(miller, &intersection);

				if (!myDet)
				{
					continue;
				}

				myDet->intersectionToSpotCoord(intersection, &xSpot, &ySpot);

				if (x == 0 && y == 0)
				{
					continue;
				}

				SpotPtr spot = SpotPtr(new Spot(getImagePtr()));
				spot->setXY(xSpot, ySpot);
				spot->setFake();

				getImagePtr()->addSpotIfNotMasked(spot);
			}
		}
	}
}

std::string MtzManager::getOrientationMatrixListScript()
{
	std::ostringstream script;

	script << "rlp_size " << spotSize << std::endl;
	script << "mosaicity " << mosaicity << std::endl;
	script << "bin " << _bin << std::endl;

	updateLatestMatrix();
	std::string description = rotatedMatrix->description(true);
	script << description << std::endl;

	return script.str();
}

void MtzManager::applyScaleFactorsForBins(int binCount)
{
	double highRes = maxResolution();
	double lowRes = 0;
//Where is partiality applied?
	vector<double> bins;
	StatisticsManager::generateResolutionBins(lowRes, highRes, binCount, &bins);

	for (int shell = 0; shell < bins.size() - 1; shell++)
	{
		double low = bins[shell];
		double high = bins[shell + 1];

		vector<ReflectionPtr> refReflections, imgReflections;

		this->findCommonReflections(referenceManager, imgReflections, refReflections,
									NULL);

		double weights = 0;
		double refMean = 0;
		double imgMean = 0;
		int count = 0;

		for (int i = 0; i < imgReflections.size(); i++)
		{
			if (!imgReflections[i]->anyAccepted())
				continue;

			if (imgReflections[i]->betweenResolutions(low, high))
			{
				weights++;
				double refIntensity = refReflections[i]->meanIntensity();
				double imgIntensity = imgReflections[i]->meanIntensity();

				if (refIntensity != refIntensity || imgIntensity != imgIntensity)
					continue;

				refMean += refIntensity;
				imgMean += imgIntensity;
				count++;
			}
		}

		refMean /= weights;
		imgMean /= weights;

		double ratio = refMean / imgMean;

		applyScaleFactor(ratio, low, high, 1);
	}
}

void MtzManager::calcXYOffset()
{
	int frameOffset = FileParser::getKey("FRAMES_OFFSET", 0);
	int frameRepeat = FileParser::getKey("FRAMES_PER_ROW", -1);

	if (frameRepeat < 0)
	{
		return;
	}

	int frameNumber = getImagePtr()->getFrameNumber();

	int minusOffset = frameNumber - frameOffset;
	double yRowFrac = (double)minusOffset / (double)frameRepeat;
	int yRow = yRowFrac;
	int restOfX = (minusOffset % frameRepeat);

	bool odd = (yRow % 2 ? false : true);

	if (odd)
	{
		restOfX = frameRepeat - restOfX;
	}

	_xPos = restOfX;
	_yPos = yRow;
}

