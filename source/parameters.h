/*
 * parameters.h
 *
 *  Created on: 17 Nov 2014
 *      Author: helenginn
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#define PARTIAL_CUTOFF 0.30
#define MULTIPLIER 512
#define OFFSET 256

#define PARAM_HROT 0
#define PARAM_KROT 1
#define PARAM_MOS 2
#define PARAM_SPOT_SIZE 3
#define PARAM_WAVELENGTH 4
#define PARAM_BANDWIDTH 5
#define PARAM_EXPONENT 6
#define PARAM_B_FACTOR 7
#define PARAM_SCALE_FACTOR 8
#define PARAM_NUM 9

#define INITIAL_BANDWIDTH 0.0013
#define INITIAL_MOSAICITY 0.0
#define INITIAL_SPOT_SIZE 0.0001
#define INITIAL_EXPONENT 1.5

#define BANDWIDTH_STEP 0.0003
#define EXPONENT_STEP 0.1
#define SPOT_STEP 0.00002
#define ROT_STEP 0.06
#define MEAN_STEP 0.005
#define MOS_STEP 0.001

#define USE_FIXED_WAVELENGTH false

#define OPTIMISED_ROT false
#define OPTIMISED_BANDWIDTH true
#define OPTIMISED_MOSAICITY true
#define OPTIMISED_EXPONENT true
#define OPTIMISED_WAVELENGTH false
#define OPTIMISED_SPOT_SIZE false

#define ROT_TOLERANCE 0.0001
#define MEAN_TOLERANCE 0.00001
#define BANDWIDTH_TOLERANCE 0.00001
#define EXPO_TOLERANCE 0.005
#define SPOT_SIZE_TOLERANCE 0.00001
#define MOS_TOLERANCE 0.001

#define MAX_OPTIMISATION_RESOLUTION 1.4
#define MAX_SPOT_SIZE_OPT_RESOLUTION 3.5
#define FURTHER_OPTIMISATION false
#define DEFAULT_SCORE_TYPE ScoreTypeMinimizeRSplit
#define CORRELATION_THRESHOLD 0.9
#define INITIAL_CORRELATION_THRESHOLD 0.8
#define THRESHOLD_SWAP 1

#define BEAM_CENTRE_X 882
#define BEAM_CENTRE_Y 882
#define DEFAULT_DETECTOR_DISTANCE 90.3
#define DEFAULT_WAVELENGTH 1.317
#define INTENSITY_THRESHOLD 12.0

#define POLARISATION_CORRECTION false
#define HORIZONTAL_POLARISATION_FACTOR 0.0

#define MINIMUM_REFLECTION_CUTOFF 70
#define REJECTING_MILLERS true
#define MIN_MILLER_COUNT 5
#define OUTLIER_REJECTION_SIGMA 1.8
#define SCALING_STRATEGY ScalingTypeReference

#define OVER_PRED_BANDWIDTH 0.03
#define OVER_PRED_SPOT_SIZE 0.0002
#define INDEXING_ORIENTATION_TOLERANCE 0.001
#define INITIAL_ORIENTATION_STEP 1.0
#define MM_PER_PIXEL 0.11
#define METROLOGY_SEARCH_SIZE 5
#define SHOEBOX_FOREGROUND_RADIUS 1
#define SHOEBOX_NEITHER_RADIUS 2
#define SHOEBOX_BACKGROUND_RADIUS 3
#define MAX_INTEGRATED_RESOLUTION 1.4
#define WINDOW_SIZE 4
#define SIGMA_RESOLUTION_CUTOFF 0.0

#define PIXEL_LEAK 2.5
#define SHOEBOX_BANDWIDTH_MULTIPLIER 0.33

#define SUPER_GAUSSIAN_STEP 0.0001
#define MAX_SUPER_GAUSSIAN 10

class Miller;

#include <sstream>
#include <memory>
#include <vector>
#include <map>
#include <boost/variant.hpp>
#include <string>

typedef enum
{
	MaskForeground = 2, MaskBackground = 1, MaskNeither = 0
} Mask;

typedef enum
{
	RFactorNone, RFactorTypeMerge, RFactorTypePim, RFactorTypeMeas,
} RFactorType;

class MtzManager;
class Panel;
class Logger;
class Image;
class Shoebox;
class Matrix;
class Indexer;

typedef std::shared_ptr<Miller> MillerPtr;
typedef std::shared_ptr<Shoebox>ShoeboxPtr;
typedef std::shared_ptr<Panel>PanelPtr;
typedef std::shared_ptr<MtzManager>MtzPtr;
typedef std::shared_ptr<std::ostringstream> StreamPtr;
typedef std::shared_ptr<Logger>LoggerPtr;
typedef std::shared_ptr<Image>ImagePtr;
typedef std::shared_ptr<Matrix>MatrixPtr;
typedef std::shared_ptr<Indexer>IndexerPtr;

typedef boost::variant<double, double, std::string, bool, int,
		std::vector<double>, std::vector<int> > ParameterVariant;
typedef std::map<std::string, ParameterVariant> ParametersMap;
typedef void (*ParserFunction)(ParametersMap *, std::string, std::string);
typedef std::map<std::string, ParserFunction> ParserMap;

typedef double (StatisticsFunction)(MtzManager *, MtzManager *, int, int *,
		double *, double, double, bool);
typedef double (RFactorFunction)(RFactorType, MtzManager *, int *, double *,
		double, double);

typedef enum
{
	WeightTypeAverage,
	WeightTypePartiality,
	WeightTypePartialitySigma,
	WeightTypeISigI,
	WeightTypePartialityCorrelation
} WeightType;

typedef enum
{
	PartialityModelNone, PartialityModelSimple, PartialityModelScaled, PartialityModelFixed
} PartialityModel;

#endif /* PARAMETERS_H_ */
