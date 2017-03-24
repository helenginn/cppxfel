/*
 * parameters.h
 *
 *  Created on: 17 Nov 2014
 *      Author: helenginn
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#define PARTIAL_CUTOFF 0.2
#define MULTIPLIER 600
#define OFFSET 300
#define DISTANCE_BUFFER 0.025

#define PARAM_HROT 0
#define PARAM_KROT 1
#define PARAM_MOS 2
#define PARAM_SPOT_SIZE 3
#define PARAM_WAVELENGTH 4
#define PARAM_BANDWIDTH 5
#define PARAM_EXPONENT 6
#define PARAM_B_FACTOR 7
#define PARAM_SCALE_FACTOR 8
#define PARAM_UNIT_CELL_A 9
#define PARAM_UNIT_CELL_B 10
#define PARAM_UNIT_CELL_C 11
#define PARAM_NUM 12

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

#define ROT_TOLERANCE 0.0002
#define MEAN_TOLERANCE 0.0001
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
#define INITIAL_ORIENTATION_STEP 0.5
#define MM_PER_PIXEL 0.11
#define METROLOGY_SEARCH_SIZE 5
#define SHOEBOX_FOREGROUND_PADDING 1
#define SHOEBOX_NEITHER_PADDING 2
#define SHOEBOX_BACKGROUND_PADDING 3
#define MAX_INTEGRATED_RESOLUTION 1.4
#define WINDOW_SIZE 4
#define SIGMA_RESOLUTION_CUTOFF 0.0

#define PIXEL_LEAK 2.5
#define SHOEBOX_BANDWIDTH_MULTIPLIER 0.33
#define REFERENCE_WEAK_REFLECTION 0

#define SUPER_GAUSSIAN_STEP 0.0001
#define MAX_SUPER_GAUSSIAN 10

class Miller;

#define member_size(type, member) sizeof(((type *)0)->member)

#include "parameters.h"
#include <sstream>
#include <memory>
#include <vector>
#include <map>
#include <boost/variant.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <mutex>

typedef std::pair<double, double> Coord;

using std::vector;


typedef struct
{
    double h;
    double k;
    double l;
} vec;

typedef enum
{
    RotationModeHorizontalVertical = 0,
    RotationModeUnitCellABC = 1,
} RotationMode;

typedef enum
{
	MaskForeground = 2, MaskBackground = 1, MaskNeither = 0
} Mask;

typedef enum
{
	RFactorNone, RFactorTypeMerge, RFactorTypePim, RFactorTypeMeas,
} RFactorType;

typedef enum
{
    MinimizationMethodStepSearch = 0,
    MinimizationMethodNelderMead = 1,
} MinimizationMethod;

typedef enum
{
    RejectReasonNone = 0, // not rejected
    RejectReasonMerge = 1 << 0,
    RejectReasonPartiality = 1 << 1,
    RejectReasonCorrelation = 1 << 2,
} RejectReason;

typedef enum
{
    ScalingTypeAverage = 0,
    ScalingTypeReference = 1,
    ScalingTypeReferenceLeastSquares = 2,
    ScalingTypeMinimizeRMerge = 3,
    ScalingTypeBFactor = 4,
    ScalingTypeResolutionShells = 5
} ScalingType;

typedef enum
{
    GeometryFormatCppxfel = 0,
    GeometryFormatCrystFEL = 1,
    GeometryFormatPanelList = 2,
} GeometryFormat;

class Detector;
class SolventMask;
class Logger;
class Image;
class IndexManager;
class Hdf5Image;
class Shoebox;
class Matrix;
class IOMRefiner;
class MtzManager;
class Spot;
class SpotVector;
class IndexingSolution;
class UnitCellLattice;
class Beam;
class GaussianBeam;
class SpectrumBeam;
class RefinementGridSearch;
class RefinementStepSearch;
class RefinementStrategy;
class Hdf5ManagerProcessing;
class Hdf5ManagerCheetahLCLS;
class Hdf5ManagerCheetahSacla;
class Hdf5ManagerCheetah;
class Hdf5Crystal;
class PNGFile;
class TextManager;
class CSV;
class SpotFinderQuick;
class SpotFinder;
class Reflection;
class NelderMead;

typedef boost::shared_ptr<SpectrumBeam> SpectrumBeamPtr;
typedef boost::shared_ptr<RefinementStepSearch> RefinementStepSearchPtr;
typedef boost::shared_ptr<RefinementGridSearch> RefinementGridSearchPtr;
typedef boost::shared_ptr<RefinementStrategy> RefinementStrategyPtr;
typedef boost::shared_ptr<NelderMead> NelderMeadPtr;
typedef boost::shared_ptr<Beam> BeamPtr;
typedef boost::shared_ptr<GaussianBeam> GaussianBeamPtr;
typedef boost::shared_ptr<Miller> MillerPtr;
typedef boost::weak_ptr<Miller> MillerWeakPtr;
typedef boost::shared_ptr<Reflection> ReflectionPtr;
typedef boost::weak_ptr<Reflection> ReflectionWeakPtr;
typedef boost::shared_ptr<Shoebox>ShoeboxPtr;
typedef boost::shared_ptr<Spot> SpotPtr;
typedef boost::shared_ptr<Detector>DetectorPtr;
typedef boost::weak_ptr<Detector>DetectorWeakPtr;
typedef boost::shared_ptr<SolventMask> SolventMaskPtr;
typedef boost::shared_ptr<MtzManager> MtzPtr;
typedef boost::shared_ptr<Hdf5Crystal> Hdf5CrystalPtr;
typedef boost::shared_ptr<std::ostringstream> StreamPtr;
typedef boost::shared_ptr<Logger>LoggerPtr;
typedef boost::shared_ptr<Image> ImagePtr;
typedef boost::weak_ptr<Image>ImageWeakPtr;
typedef boost::shared_ptr<Hdf5Image> Hdf5ImagePtr;
typedef boost::shared_ptr<Matrix>MatrixPtr;
typedef boost::shared_ptr<IOMRefiner>IOMRefinerPtr;
typedef boost::weak_ptr<IOMRefiner>IOMRefinerWeakPtr;
typedef boost::shared_ptr<SpotVector> SpotVectorPtr;
typedef boost::shared_ptr<IndexManager> IndexManagerPtr;
typedef boost::weak_ptr<IndexManager> IndexManagerWeakPtr;
typedef boost::shared_ptr<IndexingSolution> IndexingSolutionPtr;
typedef boost::shared_ptr<std::mutex> MutexPtr;
typedef boost::shared_ptr<UnitCellLattice> UnitCellLatticePtr;
typedef boost::shared_ptr<Hdf5ManagerProcessing> Hdf5ManagerProcessingPtr;
typedef std::shared_ptr<PNGFile> PNGFilePtr;
typedef std::shared_ptr<CSV> CSVPtr;
typedef std::shared_ptr<TextManager> TextManagerPtr;
typedef std::shared_ptr<SpotFinder> SpotFinderPtr;
typedef std::shared_ptr<Hdf5ManagerCheetahSacla> Hdf5ManagerCheetahSaclaPtr;
typedef std::shared_ptr<Hdf5ManagerCheetahLCLS> Hdf5ManagerCheetahLCLSPtr;
typedef std::shared_ptr<Hdf5ManagerCheetah> Hdf5ManagerCheetahPtr;

typedef boost::variant<double, double, std::string, bool, int,
  vector<double>, vector<int>, vector<std::string> > ParameterVariant;
typedef std::map<std::string, ParameterVariant> ParametersMap;
typedef void (*ParserFunction)(ParametersMap *, std::string, std::string);
typedef std::map<std::string, ParserFunction> ParserMap;

typedef double (StatisticsFunction)(MtzManager *, MtzManager *, int, int *,
		double *, double, double, bool, bool);
typedef double (RFactorFunction)(RFactorType, MtzManager *, int *, double *,
		double, double, bool);
typedef std::pair<vec, double> VectorDistance;
typedef std::pair<std::pair<SpotVectorPtr, VectorDistance>, double> Match;

typedef double (*Getter)(void *);
typedef void (*Setter)(void *, double newValue);

typedef std::map<int, std::pair<int, int> > PowderHistogram;


typedef enum
{
    LogLevelNormal = 0,
    LogLevelDetailed = 1,
    LogLevelDebug = 2
} LogLevel;


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
    PseudoScoreTypeInvalid,
    PseudoScoreTypeAll,
    PseudoScoreTypeIntraPanel,
    PseudoScoreTypeInterPanel,
    PseudoScoreTypeAllInterPanel,
    PseudoScoreTypeAngleConsistency,
    PseudoScoreTypeBeamCentre,
} PseudoScoreType;

typedef enum
{
    GeometryScoreTypeInterMiller,
    GeometryScoreTypeIntraMiller,
    GeometryScoreTypeIntrapanel,
    GeometryScoreTypeInterpanel,
    GeometryScoreTypeIntrapanelParent,
    GeometryScoreTypeBeamCentre,
    
} GeometryScoreType;

typedef enum
{
    PseudoScoreWeightingAxisH = 0,
    PseudoScoreWeightingAxisK = 1,
    PseudoScoreWeightingAxisNone = 2,
} PseudoScoreWeightingAxis;

typedef enum
{
	PartialityModelNone, PartialityModelSimple, PartialityModelScaled, PartialityModelFixed, PartialityModelBinary
} PartialityModel;

#endif /* PARAMETERS_H_ */
