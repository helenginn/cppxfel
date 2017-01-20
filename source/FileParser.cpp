/*
 * FileParser.cpp
 *
 *  Created on: 4 Jan 2015
 *      Author: helenginn
 */

#include "FileParser.h"

#include "FileReader.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <locale>
#include <stdio.h>
#include "misc.h"

#define FILE_PARSER_CPP_

#include "MtzRefiner.h"

ParametersMap FileParser::parameters;
std::ostringstream FileParser::log;
int FileParser::threadsFound = 0;
char FileParser::splitCharMajor = ' ';
char FileParser::splitCharMinor = ' ';
std::map<std::string, std::string> FileParser::deprecatedList;

int FileParser::getMaxThreads()
{
    if (threadsFound != 0)
        return threadsFound;
    
    std::ostringstream logged;
    
    char *nslots;
    nslots = getenv("NSLOTS");
    int maxThreads = 0;
    
    if (nslots != NULL)
    {
        maxThreads = atoi(nslots);
        logged << "Using environment variable NSLOTS: " << maxThreads << " threads." << std::endl;
    }
    
    if (maxThreads == 0)
    {
        maxThreads = getKey("MAX_THREADS", MAX_THREADS);
        
        if (hasKey("MAX_THREADS"))
        {
            logged << "Getting number of threads from user input: " << maxThreads << " threads." << std::endl;
        }
        else
        {
            logged << "Using default number of threads: " << maxThreads << " threads." << std::endl;
        }
    }
    
    Logger::mainLogger->addStream(&logged);
    
    threadsFound = maxThreads;
    
    return maxThreads;
}

bool FileParser::hasKey(std::string key)
{
	return (parameters.count(key) > 0);
}


void FileParser::simpleFloat(ParametersMap *map, std::string command,
		std::string rest)
{
	double theFloat = atof(rest.c_str());
    
	log << "Setting double " << command << " to " << theFloat << std::endl;

	(*map)[command] = theFloat;
}

void FileParser::simpleBool(ParametersMap *map, std::string command,
		std::string rest)
{
	bool on = (rest == "ON" || rest == "on" ? 1 : 0);

	log << "Setting bool " << command << " to " << on << std::endl;

	(*map)[command] = on;
}

void FileParser::simpleString(ParametersMap *map, std::string command,
		std::string rest)
{
    trim(rest);
    
	(*map)[command] = rest;

	log << "Setting string " << command << " to " << rest << std::endl;

}

void FileParser::simpleInt(ParametersMap *map, std::string command,
		std::string rest)
{
	int theInt = atoi(rest.c_str());

	log << "Setting int " << command << " to " << theInt << std::endl;

	(*map)[command] = theInt;
}

void FileParser::stringVector(ParametersMap *map, std::string command,
                              std::string rest)
{
    trim(rest);
    
    vector<std::string> stringVector = FileReader::split(rest, splitCharMinor);
    
    log << "Setting " << command << " to ";
    
    for (int i = 0; i < stringVector.size(); i++)
    {
        log << stringVector[i] << " ";
    }
    
    log << std::endl;
    
    (*map)[command] = stringVector;
}

void FileParser::doubleVector(ParametersMap *map, std::string command,
                              std::string rest)
{
	vector<std::string> components = FileReader::split(rest, splitCharMinor);
	vector<double> doubleVector;

	log << "Setting " << command << " to ";

	for (int i = 0; i < components.size(); i++)
	{
        if (components[i].length())
        {
            double theFloat = atof(components[i].c_str());
            log << theFloat << " ";
            doubleVector.push_back(theFloat);
        }
	}

	log << std::endl;

	(*map)[command] = doubleVector;
}

void FileParser::intVector(ParametersMap *map, std::string command,
		std::string rest)
{
	vector<std::string> components = FileReader::split(rest, splitCharMinor);
	vector<int> intVector;

	log << "Setting " << command << " to ";

	for (int i = 0; i < components.size(); i++)
    {
        if (components[i].length())
        {
            int theInt = atoi(components[i].c_str());
            log << theInt << " ";
            intVector.push_back(theInt);
        }
	}

	log << std::endl;

	(*map)[command] = intVector;
}

void FileParser::generateDeprecatedList()
{
    deprecatedList["MAX_MILLER_INDEX_TRIAL"] = "Please specify your maximum reciprocal distance using MAX_RECIPROCAL_DISTANCE.";
    deprecatedList["MINIMUM_TRUST_DISTANCE"] = "This value is automatically calculated from INITIAL_RLP_SIZE. Please use this instead.";
    deprecatedList["RECIPROCAL_TOLERANCE"] = "This value is now provided by INITIAL_RLP_SIZE. Please use this instead.";
    deprecatedList["ROTATION_MODE"] = "This option has been removed as it provided no benefit.";
    deprecatedList["INITIAL_GRID_SEARCH"] = "This option has been removed as it provided no benefit.";
    deprecatedList["LANDSCAPE_DIVISIONS"] = "This option has been removed as it provided no benefit.";
    deprecatedList["PENALTY_WEIGHT"] = "This option has been removed as it provided no benefit.";
    deprecatedList["PENALTY_RESOLUTION"] = "This option has been removed as it provided no benefit.";
    deprecatedList["MAX_RESOLUTION_ALL"] = "This option has been renamed to MAX_REFINED_RESOLUTION.";
    deprecatedList["MAX_RESOLUTION_RLP_SIZE"] = "This option has been removed as it provided no benefit.";
    deprecatedList["REFINE_B_FACTOR"] = "This option has been replaced by assigning SCALING_STRATEGY 4. For now.";
    deprecatedList["HDF5_OUTPUT_FILE"] = "This option is not supported yet.";
    deprecatedList["REJECT_IF_SPOT_COUNT"] = "This option was badly worded. Please change this to REJECT_OVER_SPOT_COUNT.";
}

void FileParser::generateFunctionList()
{
	parserMap = ParserMap();

    parserMap["VERBOSITY_LEVEL"] = simpleInt;
    parserMap["MAX_THREADS"] = simpleInt;
    
	// Refinement parameters
    parserMap["MINIMUM_CYCLES"] = simpleInt;
    parserMap["MAXIMUM_CYCLES"] = simpleInt;
    parserMap["STOP_REFINEMENT"] = simpleBool;
    parserMap["OLD_MERGE"] = simpleBool;
    parserMap["MERGE_MEDIAN"] = simpleBool;
    parserMap["READ_REFINED_MTZS"] = simpleBool;
    
    parserMap["SET_SIGMA_TO_UNITY"] = simpleBool;
    parserMap["APPLY_UNREFINED_PARTIALITY"] = simpleBool;
    parserMap["BINARY_PARTIALITY"] = simpleBool;
    parserMap["MINIMIZATION_METHOD"] = simpleInt;
    parserMap["NELDER_MEAD_CYCLES"] = simpleInt;
    parserMap["MEDIAN_WAVELENGTH"] = simpleBool;
    parserMap["WAVELENGTH_RANGE"] = doubleVector;
    parserMap["WAVELENGTH_FROM_REF_COUNT"] = simpleInt;
    parserMap["EXCLUSION_BY_CC_HALF"] = simpleBool;
    parserMap["ACCEPTABLE_UNIT_CELL_TOLERANCE"] = simpleFloat;
    parserMap["ALLOW_TRUST"] = simpleBool;
    parserMap["EXCLUDE_OWN_REFLECTIONS"] = simpleBool;
    parserMap["PARTIALITY_CUTOFF"] = simpleFloat;
	parserMap["DEFAULT_TARGET_FUNCTION"] = simpleInt;
    parserMap["TARGET_FUNCTIONS"] = intVector;
    parserMap["RLP_MODEL"] = simpleInt;
	parserMap["CORRELATION_THRESHOLD"] = simpleFloat;
    parserMap["PARTIALITY_CORRELATION_THRESHOLD"] = simpleFloat;
	parserMap["MAX_REFINED_RESOLUTION"] = simpleFloat;
    parserMap["MERGE_TO_RESOLUTION"] = simpleFloat;
    parserMap["MIN_REFINED_RESOLUTION"] = simpleFloat;
    parserMap["INITIAL_CORRELATION_THRESHOLD"] = simpleFloat;
	parserMap["THRESHOLD_SWAP"] = simpleInt;
	parserMap["OUTLIER_REJECTION_SIGMA"] = simpleFloat;
	parserMap["OUTLIER_REJECTION"] = simpleBool;
	parserMap["CORRELATION_REJECTION"] = simpleBool;
	parserMap["PARTIALITY_REJECTION"] = simpleBool;
    parserMap["POLARISATION_CORRECTION"] = simpleBool;
    parserMap["POLARISATION_FACTOR"] = simpleFloat;
    parserMap["REFINEMENT_INTENSITY_THRESHOLD"] = simpleFloat;
    parserMap["TRUST_INDEXING_SOLUTION"] = simpleBool;
    parserMap["CUSTOM_AMBIGUITY"] = doubleVector;
    parserMap["R_FACTOR_THRESHOLD"] = simpleFloat;
    parserMap["REINITIALISE_WAVELENGTH"] = simpleBool;
    parserMap["PARTIALITY_SLICES"] = simpleInt;
    parserMap["MAX_SLICES"] = simpleInt;
    parserMap["CAREFUL_RESOLUTION"] = simpleFloat;
    parserMap["SMOOTH_FUNCTION"] = simpleBool;
    parserMap["NORMALISE_PARTIALITIES"] = simpleBool;
    parserMap["REPLACE_REFERENCE"] = simpleBool;
    parserMap["REFINE_ENERGY_SPECTRUM"] = simpleBool;
    parserMap["FREE_MILLER_LIST"] = simpleString;
    parserMap["FREE_MILLER_PROPORTION"] = simpleFloat;
    parserMap["USE_HDF5_WAVELENGTH"] = simpleBool;
    parserMap["FREE_ELECTRON_LASER"] = simpleInt;
    parserMap["CHEETAH_DATA_ADDRESSES"] = simpleString;
    parserMap["CHEETAH_ID_ADDRESSES"] = simpleString;
    
	parserMap["INITIAL_WAVELENGTH"] = simpleFloat;
	parserMap["INITIAL_BANDWIDTH"] = simpleFloat;
	parserMap["INITIAL_MOSAICITY"] = simpleFloat;
	parserMap["INITIAL_EXPONENT"] = simpleFloat;
	parserMap["INITIAL_RLP_SIZE"] = simpleFloat;

	parserMap["STEP_SIZE_WAVELENGTH"] = simpleFloat;
	parserMap["STEP_SIZE_BANDWIDTH"] = simpleFloat;
	parserMap["STEP_SIZE_MOSAICITY"] = simpleFloat;
	parserMap["STEP_SIZE_EXPONENT"] = simpleFloat;
	parserMap["STEP_SIZE_ORIENTATION"] = simpleFloat;
    parserMap["STEP_SIZE_ORIENTATION_ABC"] = simpleFloat;
    parserMap["STEP_SIZE_RLP_SIZE"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_A"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_B"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_C"] = simpleFloat;

	parserMap["TOLERANCE_WAVELENGTH"] = simpleFloat;
	parserMap["TOLERANCE_BANDWIDTH"] = simpleFloat;
	parserMap["TOLERANCE_MOSAICITY"] = simpleFloat;
	parserMap["TOLERANCE_EXPONENT"] = simpleFloat;
	parserMap["TOLERANCE_ORIENTATION"] = simpleFloat;
	parserMap["TOLERANCE_RLP_SIZE"] = simpleFloat;

	parserMap["OPTIMISING_WAVELENGTH"] = simpleBool;
	parserMap["OPTIMISING_BANDWIDTH"] = simpleBool;
	parserMap["OPTIMISING_MOSAICITY"] = simpleBool;
	parserMap["OPTIMISING_EXPONENT"] = simpleBool;
	parserMap["OPTIMISING_ORIENTATION"] = simpleBool;
	parserMap["OPTIMISING_RLP_SIZE"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_A"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_B"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_C"] = simpleBool;

    parserMap["HDF5_SOURCE_FILES"] = stringVector;
    parserMap["HDF5_OUTPUT_FILE"] = simpleString;
    parserMap["DUMP_IMAGES"] = simpleBool;
    
    parserMap["ORIENTATION_MATRIX_LIST"] = simpleString;
    parserMap["SECOND_MATRIX_LIST"] = simpleString;
    parserMap["OUTPUT_DIRECTORY"] = simpleString;
    parserMap["OUTPUT_INDIVIDUAL_CYCLES"] = simpleBool;
    parserMap["MATRIX_LIST_VERSION"] = simpleFloat;
	parserMap["INITIAL_MTZ"] = simpleString;
    parserMap["IMAGE_LIMIT"] = simpleInt;
    parserMap["IMAGE_SKIP"] = simpleInt;
    parserMap["NEW_MATRIX_LIST"] = simpleString;
    
    parserMap["RECALCULATE_WAVELENGTHS"] = simpleBool;
    parserMap["RECALCULATE_SIGMA"] = simpleBool;
    parserMap["MERGE_ANOMALOUS"] = simpleBool;
    parserMap["SCALING_STRATEGY"] = simpleInt;
	parserMap["MINIMUM_REFLECTION_CUTOFF"] = simpleInt;
    parserMap["MINIMUM_MULTIPLICITY"] = simpleInt;
    parserMap["FAST_MERGE"] = simpleBool;
    parserMap["CORRECTED_PARTIALITY_MODEL"] = simpleBool;
    parserMap["REJECT_BELOW_SCALE"] = simpleFloat;

	// Indexing parameters

    parserMap["DETECTOR_GAIN"] = simpleFloat;
    parserMap["BITS_PER_PIXEL"] = simpleInt;
	parserMap["SPACE_GROUP"] = simpleInt;
	parserMap["INTEGRATION_WAVELENGTH"] = simpleFloat;
	parserMap["DETECTOR_DISTANCE"] = simpleFloat;
    parserMap["BEAM_CENTRE"] = doubleVector;
    parserMap["MM_PER_PIXEL"] = simpleFloat;
    parserMap["DETECTOR_SIZE"] = doubleVector;
	parserMap["OVER_PRED_BANDWIDTH"] = simpleFloat;
	parserMap["OVER_PRED_RLP_SIZE"] = simpleFloat;
    parserMap["REFINE_ORIENTATIONS"] = simpleBool;
    parserMap["INDEXING_ORIENTATION_TOLERANCE"] = simpleFloat;
    parserMap["LOW_INTENSITY_PENALTY"] = simpleBool;
	parserMap["INTENSITY_THRESHOLD"] = simpleFloat;
    parserMap["ABSOLUTE_INTENSITY"] = simpleBool;
	parserMap["METROLOGY_SEARCH_SIZE"] = simpleInt;
    parserMap["METROLOGY_SEARCH_SIZE_BIG"] = simpleInt;
    parserMap["METROLOGY_MOVE_THRESHOLD"] = simpleFloat;
    parserMap["FOCUS_ON_PEAK_SIZE"] = simpleInt;
	parserMap["SHOEBOX_FOREGROUND_PADDING"] = simpleInt;
	parserMap["SHOEBOX_NEITHER_PADDING"] = simpleInt;
	parserMap["SHOEBOX_BACKGROUND_PADDING"] = simpleInt;
    parserMap["SHOEBOX_MAKE_EVEN"] = simpleBool;
    parserMap["MIN_INTEGRATED_RESOLUTION"] = simpleFloat;
    parserMap["MAX_INTEGRATED_RESOLUTION"] = simpleFloat;
	parserMap["UNIT_CELL"] = doubleVector;
    parserMap["FIX_UNIT_CELL"] = simpleBool;
    parserMap["ADD_MASK"] = intVector;
    parserMap["INITIAL_ORIENTATION_STEP"] = simpleFloat;
    parserMap["COMPLEX_SHOEBOX"] = simpleBool;
    parserMap["SHOEBOX_BANDWIDTH_MULTIPLIER"] = simpleFloat;
    parserMap["PIXEL_LEAK"] = simpleFloat;
    parserMap["MILLER_INDIVIDUAL_WAVELENGTHS"] = simpleBool;
    parserMap["UNBALANCED_REFLECTIONS"] = simpleBool;
    parserMap["ORIENTATION_SCORE"] = simpleInt;
    parserMap["ORIENTATION_CORRECTION"] = doubleVector;
    parserMap["IMAGE_MASKED_VALUE"] = simpleInt;
    parserMap["IMAGE_IGNORE_UNDER_VALUE"] = simpleInt;
    parserMap["SPHERE_THICKNESS"] = simpleFloat;
    parserMap["SIGMA_RESOLUTION_CUTOFF"] = simpleFloat;
    parserMap["PIXEL_COUNT_CUTOFF"] = simpleInt;
    parserMap["EXPECTED_SPOTS"] = simpleInt;
    parserMap["REFINE_AGAINST_PSEUDO"] = simpleBool;
    parserMap["INDEXING_SLICE_ANGLE"] = simpleFloat;
    parserMap["REFINE_UNIT_CELL_A"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_B"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_C"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_ALPHA"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_BETA"] = simpleBool;
    parserMap["REFINE_UNIT_CELL_GAMMA"] = simpleBool;
    parserMap["STEP_UNIT_CELL_A"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_B"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_C"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_ALPHA"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_BETA"] = simpleFloat;
    parserMap["STEP_UNIT_CELL_GAMMA"] = simpleFloat;
    parserMap["FROM_DIALS"] = simpleBool;
    parserMap["DO_NOT_REJECT_REFLECTIONS"] = simpleBool;
    parserMap["REFINE_IN_PLANE_OF_DETECTOR"] = simpleBool;
    parserMap["FIT_BACKGROUND_AS_PLANE"] = simpleBool;
    parserMap["SKIP_BAD_PIXELS"] = simpleBool;
    parserMap["ROUGH_CALCULATION"] = simpleBool;
    parserMap["BINS"] = simpleInt;

    parserMap["MINIMUM_SPOTS_EXPLAINED"] = simpleInt;
    parserMap["MINIMUM_TRUST_ANGLE"] = simpleFloat;
    parserMap["SOLUTION_ANGLE_SPREAD"] = simpleFloat;
    parserMap["REJECT_CLOSE_SPOTS"] = simpleBool;
    parserMap["THOROUGH_SOLUTION_SEARCHING"] = simpleBool;
    parserMap["MAX_SEARCH_NUMBER_MATCHES"] = simpleInt;
    parserMap["MAX_SEARCH_NUMBER_SOLUTIONS"] = simpleInt;
    parserMap["ACCEPT_ALL_SOLUTIONS"] = simpleBool;
    parserMap["INDEXING_MIN_RESOLUTION"] = simpleFloat;
    parserMap["SPOTS_PER_LATTICE"] = simpleInt;
    parserMap["RECIPROCAL_TOLERANCE"] = simpleFloat;
    parserMap["GOOD_SOLUTION_ST_DEV"] = simpleFloat;
    parserMap["BAD_SOLUTION_ST_DEV"] = simpleFloat;
    parserMap["GOOD_SOLUTION_SUM_RATIO"] = simpleFloat;
    parserMap["GOOD_SOLUTION_HIGHEST_PEAK"] = simpleInt;
    parserMap["SOLUTION_ATTEMPTS"] = simpleInt;
    parserMap["ONE_INDEXING_CYCLE_ONLY"] = simpleBool;
    parserMap["NEW_INDEXING_METHOD"] = simpleBool;
    parserMap["MAX_RECIPROCAL_DISTANCE"] = simpleFloat;
    parserMap["ALWAYS_FILTER_SPOTS"] = simpleBool;
    parserMap["MINIMUM_NEIGHBOURS"] = simpleInt;
    parserMap["MINIMUM_SOLUTION_NETWORK_COUNT"] = simpleInt;
    parserMap["SCRAMBLE_SPOT_VECTORS"] = simpleBool;
    parserMap["NETWORK_TRIAL_LIMIT"] = simpleInt;
    parserMap["INDEXING_TIME_LIMIT"] = simpleInt;
    parserMap["MAX_LATTICES_PER_IMAGE"] = simpleInt;
    parserMap["CHECKING_COMMON_SPOTS"] = simpleBool;
    parserMap["EXCLUDE_WEAKEST_SPOT_FRACTION"] = simpleFloat;
    parserMap["REJECT_UNDER_SPOT_COUNT"] = simpleInt;
    parserMap["REJECT_OVER_SPOT_COUNT"] = simpleInt;
    parserMap["POWDER_PATTERN_STEP"] = simpleFloat;
    parserMap["POWDER_PATTERN_STEP_ANGLE"] = simpleFloat;
    parserMap["BAD_SOLUTION_HIGHEST_PEAK"] = simpleInt;
    parserMap["LOW_MEMORY_MODE"] = simpleBool;
    parserMap["GEOMETRY_FORMAT"] = simpleInt;
    
    // force spot finding
    
    parserMap["IMAGE_MIN_SPOT_INTENSITY"] = simpleFloat;
    parserMap["IMAGE_MIN_CORRELATION"] = simpleFloat;
    parserMap["IMAGE_SPOT_PROBE_HEIGHT"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_BACKGROUND"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_PADDING"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_BG_PADDING"] = simpleInt;
    parserMap["PROBE_DISTANCES"] = doubleVector;
    parserMap["RECIPROCAL_UNIT_CELL"] = doubleVector;
    parserMap["FORCE_SPOT_FINDING"] = simpleBool;
    parserMap["FORCE_RECENTRE_SPOTS"] = simpleBool;
    parserMap["FORCE_RESTART_POST_REFINEMENT"] = simpleBool;
    parserMap["SPOT_FINDING_MIN_PIXELS"] = simpleInt;
    parserMap["SPOT_FINDING_SIGNAL_TO_NOISE"] = simpleFloat;
    parserMap["SPOT_FINDING_MAX_PIXELS"] = simpleInt;
    parserMap["SPOT_FINDING_ALGORITHM"] = simpleInt;
    parserMap["SPOTS_ARE_RECIPROCAL_COORDINATES"] = simpleBool;
    
    parserMap["IGNORE_MISSING_IMAGES"] = simpleBool;
    parserMap["FROM_TIFFS"] = simpleBool;
    
    parserMap["PIXEL_TOLERANCE"] = simpleFloat;
    parserMap["MINIMUM_CIRCLE_SPOTS"] = simpleInt;
    parserMap["MAX_UNIT_CELL"] = simpleFloat;
    parserMap["RANDOM_SEED"] = simpleInt;
    
    parserMap["PNG_TOTAL"] = simpleInt;
    parserMap["PNG_THRESHOLD"] = simpleFloat;
    parserMap["PNG_SHOEBOX"] = simpleBool;
    parserMap["PNG_ALL_LATTICES"] = simpleBool;
    parserMap["PNG_HEIGHT"] = simpleInt;
    parserMap["ENABLE_IMAGE_CSVS"] = simpleBool;
    
	parserMap["PANEL_LIST"] = simpleString;
    parserMap["DETECTOR_LIST"] = simpleString;
    parserMap["SKIP_LINES"] = simpleInt;

    
    parserMap["MILLER_INDEX"] = intVector;
}

ParserFunction FileParser::splitLine(std::string line, std::string &command,
		std::string &rest)
{
    int space_index = (int)line.find_first_of(splitCharMajor);
    
    if (space_index == std::string::npos)
    {
        command = "NULL";
        return NULL;
    }
    
	command = line.substr(0, space_index);

	std::ostringstream stream;

	std::locale theLocale;
	for (std::string::size_type j = 0; j < command.length(); ++j)
		stream << std::toupper(command[j], theLocale);

	std::string upperCommand = stream.str();

	rest = line.substr(space_index + 1, std::string::npos);

    if (deprecatedList.count(command) > 0)
    {
        logged << "Deprecated command: " << command << std::endl;
        sendLog();
        logged << deprecatedList[command] << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelNormal, true);
    }
    
    // make a list of permitted words soon
    else if (parserMap.count(upperCommand) == 0 && upperCommand != "PANEL" && upperCommand != "MASK" && upperCommand != "SOLVENT_MASK")
    {
        logged << "Error: do not understand command " << upperCommand << std::endl;
        sendLog();
    }
    
    ParserFunction function = parserMap[upperCommand];
	command = upperCommand;

	return function;
}

bool FileParser::checkSpaces(std::string line)
{
    int space_index = (int)line.find_first_of(splitCharMajor);

	if (space_index == std::string::npos)
	{
		log << "Warning: " << line << " has no assignment" << std::endl;
		return false;
	}

	return true;
}

FileParser::FileParser(void)
{
    
}

FileParser::FileParser(std::string name, std::vector<std::string> someExtras)
{
	this->filename = name;
	generateFunctionList();
    generateDeprecatedList();
    
    extras = someExtras;
}

FileParser::~FileParser()
{
	// TODO Auto-generated destructor stub
}

