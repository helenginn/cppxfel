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
#include <iomanip>

#define FILE_PARSER_CPP_

#include "MtzRefiner.h"

ParametersMap FileParser::parameters;
std::ostringstream FileParser::log;
int FileParser::threadsFound = 0;
char FileParser::splitCharMajor = ' ';
char FileParser::splitCharMinor = ' ';
ParserMap FileParser::parserMap;
std::map<std::string, std::string> FileParser::deprecatedList;
std::map<std::string, std::string> FileParser::helpMap;
std::map<std::string, CodeMap> FileParser::codeMaps;
CategoryTree FileParser::categoryTree;

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

void FileParser::printAllCommands()
{
    std::ostringstream logged;
    logged << "------------------------------------" << std::endl;
    logged << "--- List of all cppxfel commands ---" << std::endl;
    logged << "------------------------------------" << std::endl << std::endl;
    
    std::vector<std::string> categorised;
    
    for (CategoryTree::iterator treeit = categoryTree.begin(); treeit != categoryTree.end(); treeit++)
    {
        std::string treeString = treeit->first;
        
        logged << treeString << std::endl << std::endl;
        
        for (CategoryMap::iterator mapit = treeit->second.begin(); mapit != treeit->second.end(); mapit++)
        {
            std::string mapString = mapit->first;
            
            logged << "    " << mapString << std::endl << std::endl;
            
            for (int i = 0; i < mapit->second.size(); i++)
            {
                std::string command = mapit->second[i];
                
                if (!parserMap.count(command))
                {
                    logged << "(!! cannot find " << command << " in parameter list!)" << std::endl;
                }
                
                logged << "        ";
                
                bool helpExists = helpMap.count(command);
                logged << std::setw(35) << command << (helpExists ? " --help exists." : "") << std::endl;

                categorised.push_back(command);
            }
            
            logged << std::endl;
        }
    }

	Logger::log(logged);
    
    logged << std::endl << "Uncategorised commands below:" << std::endl << std::endl;
    
    for (ParserMap::iterator it = parserMap.begin(); it != parserMap.end(); it++)
    {
        std::vector<std::string>::iterator done_it;
        done_it = std::find(categorised.begin(), categorised.end(), it->first);
        
        if (*done_it == it->first)
        {
            continue;
        }
        
        bool helpExists = helpMap.count(it->first);
        logged << std::setw(35) << it->first << (helpExists ? " --help exists." : "") << std::endl;
    }
    
    Logger::log(logged);
}

void FileParser::printCommandInfo(std::string unsanitisedCommand)
{
    std::ostringstream commandStr;
    std::locale loc;
    for (int i = 0; i < unsanitisedCommand.length(); i++)
    {
        commandStr << std::toupper(unsanitisedCommand[i], loc);
    }
    std::string command = commandStr.str();
    
    std::ostringstream logged;
    
    if (helpMap.count(command)) {
        logged << "There is a help statement for this command." << std::endl << std::endl;
        logged << command << ": " << helpMap[command] << std::endl << std::endl;
    }
    else
    {
        if (parserMap.count(command))
        {
            logged << "Sorry, there doesn't appear to be a help statement for this command, though it does exist..." << std::endl << std::endl;
        }
        else
        {
            logged << "The command " << command << " is not understood by cppxfel." << std::endl << std::endl;
            Logger::log(logged);
            return;
        }

    }
    
    if (codeMaps.count(command))
    {
        logged << "Allowed options for " << command << ":" << std::endl;
        CodeMap map = codeMaps[command];
        
        for (CodeMap::iterator it = map.begin(); it != map.end(); it++)
        {
            logged << "   - " << it->first << std::endl;
        }
        
        logged << std::endl;
        
    }
    else
    {
        std::string expectation = "[unknown]";
        
        if (parserMap[command] == simpleBool)
        {
            expectation = "a boolean (\"ON\" or \"OFF\")";
        }
        else if (parserMap[command] == simpleFloat)
        {
            expectation = "a double (a potentially non-integer number)";
        }
        else if (parserMap[command] == simpleInt)
        {
            expectation = "an integer";
        }
        else if (parserMap[command] == simpleString)
        {
            expectation = "a string (text)";
        }
        else if (parserMap[command] == doubleVector)
        {
            expectation = "a series of potentially non-integer numbers separated by spaces.";
        }
        else if (parserMap[command] == intVector)
        {
            expectation = "a series of integers separated by spaces.";
        }
        else if (parserMap[command] == stringVector)
        {
            expectation = "one or more strings (text) separated by spaces.";
        }
        
        logged << "Command " << command << " is expecting " << expectation << std::endl << std::endl;
    }
    
    Logger::log(logged);
}

void FileParser::simpleInt(ParametersMap *map, std::string command,
		std::string rest)
{
    if (codeMaps.count(command) && codeMaps[command].size() > 0)
    {
        CodeMap codeMap = codeMaps[command];
        std::locale loc;
        std::ostringstream lowerCaseStream;
        for (int i = 0; i < rest.length(); i++)
            lowerCaseStream << std::tolower(rest[i], loc);
        std::string lowered = lowerCaseStream.str();
        
        if (!codeMap.count(lowered))
        {
            std::ostringstream logged;
            logged << "You used the code \"" << lowered << "\" to set parameter " << command << std::endl;
            logged << "This is not one of the allowed options. Help says:" << std::endl << std::endl;
            Logger::log(logged);
            printCommandInfo(command);
            logged << "Please edit your log file, for example, to:" << std::endl;
            logged << command << " " << codeMap.begin()->first << std::endl << std::endl;
            staticLogAndExit(logged);
            return;
        }

        int theInt = codeMap[lowered];
        (*map)[command] = theInt;
        
        return;
    }
    
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
    deprecatedList["REFINE_B_FACTOR"] = "This option has been replaced by assigning SCALING_STRATEGY b_factor. For now.";
    deprecatedList["HDF5_OUTPUT_FILE"] = "This option is not supported yet.";
    deprecatedList["REJECT_IF_SPOT_COUNT"] = "This option was badly worded. Please change this to REJECT_OVER_SPOT_COUNT.";
    deprecatedList["ORIENTATION_SCORE"] = "This option has been removed as it provided no benefit.";

    deprecatedList["REFINE_UNIT_CELL_A"] = "This option has been merged with OPTIMISING_UNIT_CELL_A, please use this instead.";
    deprecatedList["REFINE_UNIT_CELL_B"] = "This option has been merged with OPTIMISING_UNIT_CELL_B, please use this instead.";
    deprecatedList["REFINE_UNIT_CELL_C"] = "This option has been merged with OPTIMISING_UNIT_CELL_C, please use this instead.";
    deprecatedList["REFINE_UNIT_CELL_ALPHA"] = "This option has been merged with OPTIMISING_UNIT_CELL_ALPHA, please use this instead.";
    deprecatedList["REFINE_UNIT_CELL_BETA"] = "This option has been merged with OPTIMISING_UNIT_CELL_BETA, please use this instead.";
    deprecatedList["REFINE_UNIT_CELL_GAMMA"] = "This option has been merged with OPTIMISING_UNIT_CELL_GAMMA, please use this instead.";
    deprecatedList["ACCEPTABLE_UNIT_CELL_TOLERANCE"] = "This option has been removed as it provided no benefit.";
    deprecatedList["STEP_UNIT_CELL_A"] = "This option was badly worded. Please change this to STEP_SIZE_UNIT_CELL_A.";
    deprecatedList["STEP_UNIT_CELL_B"] = "This option was badly worded. Please change this to STEP_SIZE_UNIT_CELL_B.";
    deprecatedList["STEP_UNIT_CELL_C"] = "This option was badly worded. Please change this to STEP_SIZE_UNIT_CELL_C.";
    deprecatedList["DO_NOT_REJECT_REFLECTIONS"] = "This option was badly worded. Please change this to REJECT_OVERLAPPING_REFLECTIONS (setting may need swapping over - default is now ON).";
    deprecatedList["PANEL_LIST"] = "This option has been deprecated as the panels definitions were too basic. This can now be specified under DETECTOR_LIST and using GEOMETRY_FORMAT panel_list. However it is highly recommended that you switch to CrystFEL or ideally cppxfel format.";
    deprecatedList["FREE_MILLER_LIST"] = "This option never worked properly and has been disabled for the time being.";
    deprecatedList["RECALCULATE_SIGMA"] = "This option has been removed. The sigmas are recalculated automatically and stored in the CSIGI column of the MTZ file.";
    
    deprecatedList["PARTIALITY_SLICES"] = "The slicing of reciprocal lattice points is now determined automatically on a resolution-independent basis and so this option has been removed.";
    deprecatedList["MAX_SLICES"] = "The slicing of reciprocal lattice points is now determined automatically on a resolution-independent basis and so this option has been removed.";
    deprecatedList["CAREFUL_RESOLUTION"] = "The slicing of reciprocal lattice points is now determined automatically on a resolution-independent basis and so this option has been removed.";
    deprecatedList["RLP_MODEL"] = "This parameter made little to no difference and so has been removed.";
    deprecatedList["OLD_MERGE"] = "The old merging system has been removed as the new merging system is now identical in nature.";
	deprecatedList["THOROUGH_SOLUTION_SEARCHING"] = "This option has been removed because the network method is better than the cluster method, so the latter has been removed.";
	deprecatedList["MAX_SEARCH_NUMBER_MATCHES"] = "This option has been removed because the network method is better than the cluster method, so the latter has been removed.";
	deprecatedList["MAX_SEARCH_NUMBER_SOLUTIONS"] = "This option has been removed because the network method is better than the cluster method, so the latter has been removed.";
	deprecatedList["NEW_INDEXING_METHOD"] = "This option has been removed because the network method is better than the cluster method, so the latter has been removed.";
	deprecatedList["MINIMUM_NEIGHBOURS"] = "This option has been removed because the network method is better than the cluster method, so the latter has been removed.";

}

void FileParser::generateCodeList()
{
    {
        CodeMap codeMap;
        codeMap["cppxfel"] = 0;
        codeMap["crystfel"] = 1;
        codeMap["panel_list"] = 2;
        codeMaps["GEOMETRY_FORMAT"] = codeMap;
    }
    {
        CodeMap codeMap;
        codeMap["step_search"] = 0;
        codeMap["nelder_mead"] = 1;
        codeMap["grid_search"] = 2;
        codeMaps["MINIMIZATION_METHOD"] = codeMap;
    }
    {
        CodeMap codeMap;
        codeMap["average"] = 0;
        codeMap["reference"] = 1;
        codeMap["b_factor"] = 2;
        codeMaps["SCALING_STRATEGY"] = codeMap;
    }
    {
        CodeMap codeMap;
        codeMap["normal"] = 0;
        codeMap["detailed"] = 1;
        codeMap["debug"] = 2;
        codeMaps["VERBOSITY_LEVEL"] = codeMap;
    }
    {
        CodeMap codeMap;
        codeMap["rfactor"] = 0;
        codeMap["correl"] = 1;
        codeMap["part"] = 2;
        codeMap["r+int"] = 11;
        codeMaps["DEFAULT_TARGET_FUNCTION"] = codeMap;
    }
    {
        CodeMap codeMap;
        codeMap["uniform"] = 0;
        codeMap["gaussian"] = 1;
        codeMaps["RLP_MODEL"] = codeMap;
    }
    {
        CodeMap codeMap;
        codeMap["lcls"] = 0;
        codeMap["sacla"] = 1;
        codeMap["european_xfel"] = 2;
        codeMap["swissfel"] = 3;
        codeMaps["FREE_ELECTRON_LASER"] = codeMap;
    }
    {
        CodeMap codeMap;
        codeMap["blob"] = 0;
        codeMap["peakfinder6"] = 1;
        codeMaps["SPOT_FINDING_ALGORITHM"] = codeMap;
    }
}

void FileParser::generateCategories()
{
    CategoryMap spotFinding;
    std::vector<std::string> easySpotFindings;
    easySpotFindings.push_back("IMAGE_MIN_SPOT_INTENSITY");
    easySpotFindings.push_back("IMAGE_MIN_CORRELATION");
    easySpotFindings.push_back("IMAGE_SPOT_PROBE_PADDING");
    easySpotFindings.push_back("IMAGE_SPOT_PROBE_BG_PADDING");
    easySpotFindings.push_back("FORCE_SPOT_FINDING");
    easySpotFindings.push_back("REJECT_CLOSE_SPOTS");
    easySpotFindings.push_back("SPOT_FINDING_ALGORITHM");
    spotFinding["Basic parameters"] = easySpotFindings;
    
    std::vector<std::string> hardSpotFindings;
    hardSpotFindings.push_back("IMAGE_SPOT_PROBE_BACKGROUND");
    hardSpotFindings.push_back("IMAGE_SPOT_PROBE_HEIGHT");
    hardSpotFindings.push_back("EXCLUDE_WEAKEST_SPOT_FRACTION");
    hardSpotFindings.push_back("IMAGE_ACCEPTABLE_COMBO");
    spotFinding["Expert parameters"] = hardSpotFindings;
    
    categoryTree["Spot finding"] = spotFinding;
    
    CategoryMap indexing;
    std::vector<std::string> indexingSolutionChecks;
    indexingSolutionChecks.push_back("GOOD_SOLUTION_HIGHEST_PEAK");
    indexingSolutionChecks.push_back("GOOD_SOLUTION_ST_DEV");
    indexingSolutionChecks.push_back("GOOD_SOLUTION_SUM_RATIO");
    indexingSolutionChecks.push_back("BAD_SOLUTION_HIGHEST_PEAK");
    indexingSolutionChecks.push_back("BAD_SOLUTION_HIGHEST_PEAK");
    indexingSolutionChecks.push_back("BAD_SOLUTION_ST_DEV");
    indexingSolutionChecks.push_back("MINIMUM_SPOTS_EXPLAINED");
    
    std::vector<std::string> easyIndexingParameters;
    easyIndexingParameters.push_back("INITIAL_RLP_SIZE");
    easyIndexingParameters.push_back("MAX_RECIPROCAL_DISTANCE");
    easyIndexingParameters.push_back("INDEXING_TIME_LIMIT");
    easyIndexingParameters.push_back("SOLUTION_ATTEMPTS");

    std::vector<std::string> hardIndexingParameters;
    hardIndexingParameters.push_back("SOLUTION_ANGLE_SPREAD");
    hardIndexingParameters.push_back("MINIMUM_TRUST_ANGLE");
    hardIndexingParameters.push_back("MINIMUM_SOLUTION_NETWORK_COUNT");
    hardIndexingParameters.push_back("ACCEPT_ALL_SOLUTIONS");
	hardIndexingParameters.push_back("MAX_LATTICES_PER_IMAGE");
    hardIndexingParameters.push_back("CHECKING_COMMON_SPOTS");
    hardIndexingParameters.push_back("NETWORK_TRIAL_LIMIT");

	std::vector<std::string> powderPatternSpecific;
	powderPatternSpecific.push_back("ALWAYS_FILTER_SPOTS");
	powderPatternSpecific.push_back("SPOTS_PER_LATTICE");
	powderPatternSpecific.push_back("POWDER_PATTERN_STEP");
	powderPatternSpecific.push_back("REJECT_OVER_SPOT_COUNT");
	powderPatternSpecific.push_back("REJECT_UNDER_SPOT_COUNT");

    indexing["Indexing solution checks"] = indexingSolutionChecks;
    indexing["Basic indexing parameters"] = easyIndexingParameters;
    indexing["Expert indexing parameters"] = hardIndexingParameters;
	indexing["Powder pattern specific"] = powderPatternSpecific;

    categoryTree["Indexing"] = indexing;
    
    CategoryMap hdf5Imports;
    std::vector<std::string> easyHdf5Imports;
    easyHdf5Imports.push_back("HDF5_SOURCE_FILES");
    easyHdf5Imports.push_back("IMAGE_LIMIT");
    easyHdf5Imports.push_back("IMAGE_SKIP");
    easyHdf5Imports.push_back("ORIENTATION_MATRIX_LIST");

    std::vector<std::string> hardHdf5Imports;
    hardHdf5Imports.push_back("HDF5_AS_FLOAT");
    hardHdf5Imports.push_back("CHEETAH_DATA_ADDRESSES");
    hardHdf5Imports.push_back("HDF5_MASK_ADDRESS");
    hardHdf5Imports.push_back("CHEETAH_ID_ADDRESSES");
    hardHdf5Imports.push_back("BITS_PER_PIXEL");
    hardHdf5Imports.push_back("MATRIX_LIST_VERSION");
    hardHdf5Imports.push_back("FREE_ELECTRON_LASER");
	hardHdf5Imports.push_back("USE_HDF5_WAVELENGTH");

    hdf5Imports["Basic image import parameters"] = easyHdf5Imports;
    hdf5Imports["Expert image import parameters"] = hardHdf5Imports;
    
    categoryTree["Importing images and metadata"] = hdf5Imports;

    CategoryMap pngOutput;
    std::vector<std::string> easyPNG;
    easyPNG.push_back("PNG_SHOEBOX");
    easyPNG.push_back("PNG_STRONG_ONLY");
    easyPNG.push_back("PNG_TOTAL");
    easyPNG.push_back("PNG_HEIGHT");
    easyPNG.push_back("DRAW_GEOMETRY_PNGS");
    
    std::vector<std::string> hardPNG;
    hardPNG.push_back("PNG_THRESHOLD");
    hardPNG.push_back("PNG_ALL_LATTICES");

    pngOutput["Basic PNG parameters"] = easyPNG;
    pngOutput["Expert PNG parameters"] = hardPNG;
    
    categoryTree["Output PNG files"] = pngOutput;

    CategoryMap integration;
    std::vector<std::string> easyShoeboxes;
    easyShoeboxes.push_back("SHOEBOX_BACKGROUND_PADDING");
    easyShoeboxes.push_back("SHOEBOX_FOREGROUND_PADDING");
    easyShoeboxes.push_back("SHOEBOX_NEITHER_PADDING");
    easyShoeboxes.push_back("METROLOGY_SEARCH_SIZE");
    
    std::vector<std::string> hardShoeboxes;
    hardShoeboxes.push_back("SHOEBOX_MAKE_EVEN");
    hardShoeboxes.push_back("COMPLEX_SHOEBOX");
    hardShoeboxes.push_back("PIXEL_LEAK");
    hardShoeboxes.push_back("SHOEBOX_BANDWIDTH_MULTIPLIER");
    
    std::vector<std::string> generalIntegration;
    generalIntegration.push_back("MIN_INTEGRATED_RESOLUTION");
    generalIntegration.push_back("MAX_INTEGRATED_RESOLUTION");
    generalIntegration.push_back("MM_PER_PIXEL");
    generalIntegration.push_back("OVER_PRED_BANDWIDTH");
    generalIntegration.push_back("OVER_PRED_RLP_SIZE");
    generalIntegration.push_back("SPACE_GROUP");
    generalIntegration.push_back("UNIT_CELL");
    generalIntegration.push_back("INTENSITY_THRESHOLD");
    generalIntegration.push_back("REFINE_ORIENTATIONS");
    generalIntegration.push_back("ROUGH_CALCULATION");
    generalIntegration.push_back("INTEGRATION_WAVELENGTH");

    std::vector<std::string> hardIntegration;
    hardIntegration.push_back("REFINE_IN_PLANE_OF_DETECTOR");
    hardIntegration.push_back("SPHERE_THICKNESS");
    hardIntegration.push_back("ABSOLUTE_INTENSITY");
    hardIntegration.push_back("UNBALANCED_REFLECTIONS");
    hardIntegration.push_back("IGNORE_MISSING_IMAGES");
    
    integration["Basic integration window parameters"] = easyShoeboxes;
    integration["Expert integration window parameters"] = hardShoeboxes;
    integration["General integration parameters"] = generalIntegration;
    integration["Expert integration parameters"] = hardIntegration;
    
    categoryTree["Integration parameters"] = integration;

    CategoryMap postRefinement;
    std::vector<std::string> postRefOptimisers;
    postRefOptimisers.push_back("OPTIMISING_BANDWIDTH");
    postRefOptimisers.push_back("OPTIMISING_EXPONENT");
    postRefOptimisers.push_back("OPTIMISING_ORIENTATION");
    postRefOptimisers.push_back("OPTIMISING_RLP_SIZE");
    postRefOptimisers.push_back("OPTIMISING_WAVELENGTH");
    
    std::vector<std::string> postRefStepSizes;
    postRefStepSizes.push_back("STEP_SIZE_BANDWIDTH");
    postRefStepSizes.push_back("STEP_SIZE_EXPONENT");
    postRefStepSizes.push_back("STEP_SIZE_ORIENTATION");
    postRefStepSizes.push_back("STEP_SIZE_RLP_SIZE");
    postRefStepSizes.push_back("STEP_SIZE_WAVELENGTH");
    
    std::vector<std::string> postRefTolerances;
    postRefTolerances.push_back("TOLERANCE_BANDWIDTH");
    postRefTolerances.push_back("TOLERANCE_EXPONENT");
    postRefTolerances.push_back("TOLERANCE_ORIENTATION");
    postRefTolerances.push_back("TOLERANCE_RLP_SIZE");
    postRefTolerances.push_back("TOLERANCE_WAVELENGTH");
    
    std::vector<std::string> initialValues;
    postRefTolerances.push_back("INITIAL_BANDWIDTH");
    postRefTolerances.push_back("INITIAL_EXPONENT");
    postRefTolerances.push_back("INITIAL_RLP_SIZE");
    postRefTolerances.push_back("INITIAL_WAVELENGTH");
    
    std::vector<std::string> postRefGeneral;
    postRefGeneral.push_back("REFINEMENT_INTENSITY_THRESHOLD");
    postRefGeneral.push_back("PARTIALITY_CUTOFF");
    postRefGeneral.push_back("OUTPUT_INDIVIDUAL_CYCLES");
    postRefGeneral.push_back("DEFAULT_TARGET_FUNCTION");
    postRefGeneral.push_back("BINARY_PARTIALITY");
    postRefGeneral.push_back("INITIAL_MTZ");
    
    postRefinement["On/off optimisation switches"] = postRefOptimisers;
    postRefinement["Step sizes for parameters"] = postRefStepSizes;
    postRefinement["Post-refinement tolerance values"] = postRefTolerances;
    postRefinement["Post-refinement initial values"] = initialValues;
    postRefinement["General options"] = postRefGeneral;
    
    categoryTree["Post-refinement"] = postRefinement;
    
    CategoryMap others;
    std::vector<std::string> generalOptions;
    generalOptions.push_back("MAX_THREADS");
    generalOptions.push_back("OUTPUT_DIRECTORY");
    generalOptions.push_back("VERBOSITY_LEVEL");
    generalOptions.push_back("MAXIMUM_CYCLES");
    others["General options"] = generalOptions;

    std::vector<std::string> nonOperational;
    nonOperational.push_back("STEP_SIZE_UNIT_CELL_A");
    nonOperational.push_back("STEP_SIZE_UNIT_CELL_ALPHA");
    nonOperational.push_back("STEP_SIZE_UNIT_CELL_B");
    nonOperational.push_back("STEP_SIZE_UNIT_CELL_BETA");
    nonOperational.push_back("STEP_SIZE_UNIT_CELL_C");
    nonOperational.push_back("STEP_SIZE_UNIT_CELL_GAMMA");
    nonOperational.push_back("STEP_SIZE_MOSAICITY");
    nonOperational.push_back("INITIAL_MOSAICITY");
    nonOperational.push_back("TOLERANCE_MOSAICITY");
    nonOperational.push_back("OPTIMISING_MOSAICITY");
    nonOperational.push_back("NORMALISE_PARTIALITIES");
    nonOperational.push_back("MILLER_INDEX");
    nonOperational.push_back("FREE_MILLER_LIST");
    nonOperational.push_back("MILLER_INDEX");
    nonOperational.push_back("FIT_BACKGROUND_AS_PLANE");
    nonOperational.push_back("FREE_MILLER_PROPORTION");
    others["Non-operational but upcoming, or not entirely abandoned"] = nonOperational;

    
    categoryTree["Other"] = others;

}

void FileParser::generateHelpList()
{
    helpMap["VERBOSITY_LEVEL"] = "Sets the level of output to the terminal. Default normal.";
    helpMap["MAX_THREADS"] = "Sets maximum number of threads to use during threaded functionality. Generally set equal to the number of cores on the machine (unless you know each core can support more than one thread).";
    helpMap["MATRIX_LIST_VERSION"] = "Mostly for backwards compatibility. V2 is suitable for pickle-based data, V3 is auto-picked up and suitable for Cheetah output from the all-*.dat file.";
    helpMap["ORIENTATION_MATRIX_LIST"] = "The name of the file containing list of images (and/or their metadata and any determined crystal orientation matrices). If using HDF5 format, this list is optional (images can be loaded directly from HDF5). However it is available if you wish to load in only certain images from the HDF5 file, or reprocess indexed results. These are also created as a result of indexing, integration and so on. However if you wish to load in raw pixel streams (4-byte or 2-byte streams from .img files), this is necessary to know which files to load.\n\nA very basic version of the contents of the file, just to load the images (e.g. images.dat) would be:\n\nimage LCLS_2013_Mar16_r0003_094913_1368f\nimage LCLS_2013_Mar16_r0003_094916_13a9a\nimage LCLS_2013_Mar16_r0003_094918_13d9d\nimage LCLS_2013_Mar16_r0003_094921_1410c\nimage LCLS_2013_Mar16_r0003_094923_1443f\nimage LCLS_2013_Mar16_r0003_094925_14757\n\nThis should omit the file extension (.img) if loading from raw streams.";
    helpMap["INITIAL_MTZ"] = "MTZ used as a starting reference data set in place of an initial merge. Useful for small numbers of images or to break the indexing ambiguity.";
    helpMap["IMAGE_LIMIT"] = "Only attempt to load x images from the list of images specified in ORIENTATION_MATRIX_LIST or number of images found in HDF5_SOURCE_FILES. Default 0 (no limit).";
    helpMap["NEW_MATRIX_LIST"] = "Filename for the set of matrix lists (all-x) which may be produced from a round of integration or post-refinement. No default. At the moment other forms (refine-x, integrate-x and merge-x) are supported for back compatibility.";
    helpMap["SPACE_GROUP"] = "Space group used for indexing or integration (will update existing files if changed). Note for indexing, it is best to use the point group symmetry with no screw axes.";
    
    helpMap["MINIMUM_CYCLES"] = "Minimum number of cycles of post-refinement to exceed before finishing (but may perform more cycles if not converged). Default 6.";
    helpMap["MAXIMUM_CYCLES"] = "Integer x – maximum number of cycles of post-refinement to execute even if not converged. Default 0 (no maximum).";
    
    helpMap["STOP_REFINEMENT"] = "If set to OFF, post-refinement will continue indefinitely. Default ON.";
    helpMap["MINIMIZATION_METHOD"] = "Minimization method used for various minimization events throughout the software. Grid search NOT recommended for normal use but for debugging purposes.";
    helpMap["NELDER_MEAD_CYCLES"] = "If using Nelder Mead, specify how many cycles are carried out (convergence criteria not implemented).";
    helpMap["MEDIAN_WAVELENGTH"] = "Calculate starting X-ray beam wavelength for post-refinement of an image using the median excitation wavelength of all strong reflections. Otherwise a mean average is used. Default OFF.";
    helpMap["WAVELENGTH_RANGE"] = "x, y – start and end for range of wavelengths to consider when calculating the starting X-ray beam wavelength. Can ignore extreme outliers. Default 0 0 (not applied).";
    helpMap["REFINEMENT_INTENSITY_THRESHOLD"] = "Double x - Intensity threshold x in absolute terms to define ‘strong’ reflections in initial wavelength determination. Default 200.";
    helpMap["ALLOW_TRUST"] = "If an image correlates well with the reference data set, fix the indexing ambiguity chosen in the future to reduce computation time. Default ON";
    helpMap["TRUST_INDEXING_SOLUTION"] = "Do not attempt to check alternative indexing solutions if set to ON. Good if the indexing ambiguity has been resolved by some other means. Default OFF.";
    helpMap["PARTIALITY_CUTOFF"] = "If reflections are calculated with a cutoff below a certain partiality they are not included in target function calculation or merging. Default 0.2.";
    helpMap["SCALING_STRATEGY"] = "number representing the strategy for scaling individual crystals on each merging cycle. Default reference.";
    helpMap["MINIMUM_REFLECTION_CUTOFF"] = "If a crystal refines to have fewer than x reflections then it is not included in the final merge. Default 30.";
    helpMap["DEFAULT_TARGET_FUNCTION"] = "Target function used for post-refinement minimisation.";
    helpMap["TARGET_FUNCTIONS"] = "Not recommended right now.";
    helpMap["CORRELATION_THRESHOLD"] = "Threshold for merging images; image which correlate less than x with reference will not be included in refinement. Also see INITIAL_CORRELATION_THRESHOLD. Default 0.9.";
    helpMap["INITIAL_CORRELATION_THRESHOLD"] = "For a certain number of rounds of refinement determined by THRESHOLD_SWAP, a lower correlation threshold of x will be used. Also see CORRELATION_THRESHOLD. Default 0.8.";
    helpMap["THRESHOLD_SWAP"] = "This specifies the number of rounds of refinement before swapping from INITIAL_CORRELATION_THRESHOLD to CORRELATION_THRESHOLD. Default 2.";
    helpMap["PARTIALITY_CORRELATION_THRESHOLD"] = "Threshold for merging images, image whose partialities correlate less than x with the proportion of merged intensities of the reference will not be included in refinement. Not recommended. Default 0.";
    helpMap["R_FACTOR_THRESHOLD"] = "In addition to CORRELATION_THRESHOLD, images with an R factor of more than x with the reference data set are rejected during merging. Default is 0 (not in use).";
    helpMap["MAX_REFINED_RESOLUTION"] = "Do not use reflections in post-refinement beyond x Å resolution (but these will be included in the merge). Default 1.4 Å.";
    helpMap["MIN_REFINED_RESOLUTION"] = "Do not refine using reflections below x Å resolution (but these will be included in the merge). Default 0 (no minimum).";
    helpMap["OUTLIER_REJECTION"] = "Master switch for rejection of outliers based on standard deviations from mean. Default ON.";
    helpMap["OUTLIER_REJECTION_SIGMA"] = "Number of standard deviations away from the mean intensity of a merged reflection beyond which an observation is rejected during merging. Default 1.8.";
    helpMap["CORRELATION_REJECTION"] = "Rejection on a per image basis if individual reflections correlate poorly with the image. Good for unindexed multiple lattices or bad-pixel detectors. Default ON.";
    helpMap["POLARISATION_CORRECTION"] = "If switched on, polarisation factor is applied. Default OFF. Does this even still work?";
    helpMap["POLARISATION_FACTOR"] = "Number between 0 for fully horizontal polarisation, 1 for fully vertical polarisation.";
    helpMap["REPLACE_REFERENCE"] = "If ON, reference is updated on each merging cycle (recommended). Otherwise reference is not replaced. Only generally used to measure degree of overfitting from refining a small number of images. Default ON.";
    helpMap["INITIAL_WAVELENGTH"] = "Does not need to be set, but if set, wavelength is always set to this value and not calculated per image. Possibly useful for very low resolution structures where individual spots deviate far from the true wavelength. Default 0 (not set).";
    helpMap["INITIAL_BANDWIDTH"] = "Initial bandwidth standard deviation for beam model. Multiply this by the wavelength to get the standard deviation in Å. Default 0.0013 (calibrated for a SASE pulse).";
    helpMap["INITIAL_MOSAICITY"] = "Initial mosaicity for rlp size determining how much the rlp increases with resolution in degrees. Defined according to Rossmann et al 1979 in degrees. Default 0º.";
    helpMap["INITIAL_EXPONENT"] = "Initial exponent for super-Gaussian distribution of beam wavelengths. Default 1.5.";
    helpMap["INITIAL_RLP_SIZE"] = "Size of theoretical Miller index (0, 0, 0) in reciprocal space distance (Å-1). Default 0.0001 Å-1. Also used for determining the tolerance for inter-spot vectors during TakeTwo indexing.";
    helpMap["STEP_SIZE_WAVELENGTH"] = "Initial step size or initial displacement for step search or Nelder-Mead algorithm respectively for wavelength. Default 0.005 Å.";
    helpMap["STEP_SIZE_BANDWIDTH"] = "Initial step size or initial displacement for step search or Nelder-Mead algorithm respectively for bandwidth. Default 0.0003.";
    helpMap["STEP_SIZE_MOSAICITY"] = "Initial step size or initial displacement for step search or Nelder-Mead algorithm respectively for mosaicity. Default 0.001.";
    helpMap["STEP_SIZE_EXPONENT"] = "Initial step size or initial displacement for step search or Nelder-Mead algorithm respectively for the super-Gaussian exponent. Default 0.1.";
    helpMap["STEP_SIZE_ORIENTATION"] = "Initial step size or initial displacement for step search or Nelder-Mead algorithm respectively for horizontal/vertical rotation displacement (degrees). Default 0.06º.";
    helpMap["STEP_SIZE_RLP_SIZE"] = "Initial step size or initial displacement for step search or Nelder-Mead algorithm respectively for base rlp size. Default 2 × 10-5 Å-1.";
    helpMap["TOLERANCE_WAVELENGTH"] = "Step size search algorithm only – when step search size falls below this tolerance, stop parameter refinement; wavelength. Default 1 × 10-5 Å.";
    helpMap["TOLERANCE_BANDWIDTH"] = "Step size search algorithm only – when step search size falls below this tolerance, stop parameter refinement; bandwidth. Default 1 × 10-5.";
    helpMap["TOLERANCE_MOSAICITY"] = "Step size search algorithm only – when step search size falls below this tolerance, stop parameter refinement; mosaicity. Default 1 × 10-3º.";
    helpMap["TOLERANCE_EXPONENT"] = "Step size search algorithm only – when step search size falls below this tolerance, stop parameter refinement; super-Gaussian exponent. Default 5 × 10-3.";
    helpMap["TOLERANCE_ORIENTATION"] = "Step size search algorithm only – when step search size falls below this tolerance, stop parameter refinement; horizontal/vertical orientation displacement (degrees). Default 1 × 10-4º.";
    helpMap["TOLERANCE_RLP_SIZE"] = "Step size search algorithm only – when step search size falls below this tolerance, stop parameter refinement; base rlp size. Default 1 × 10-5 Å-1.";
    helpMap["OPTIMISING_WAVELENGTH"] = "Sets whether wavelength is a refined parameter. Default ON.";
    helpMap["OPTIMISING_BANDWIDTH"] = "Sets whether bandwidth is a refined parameter. Default OFF.";
    helpMap["OPTIMISING_MOSAICITY"] = "Sets whether mosaicity is a refined parameter. Default OFF.";
    helpMap["OPTIMISING_EXPONENT"] = "Sets whether exponent is a refined parameter. Default OFF.";
    helpMap["OPTIMISING_ORIENTATION"] = "Sets whether horizontal/vertical orientation displacement is a refined parameter. Default ON";
    helpMap["OPTIMISING_RLP_SIZE"] = "Sets whether rlp size is a refined parameter. Default ON.";
    
    helpMap["INTEGRATION_WAVELENGTH"] = "Wavelength of X-ray beam in Å. No default. Can be overridden by value in HDF5 file.";
    helpMap["DETECTOR_DISTANCE"] = "Detector distance from crystal in mm. No default. Can be overridden by value in cppxfel-geometry file.";
    helpMap["BEAM_CENTRE"] = "Coordinate position at which the beam hits the image in pixels. Default 0, 0 px to match CSPAD detector.";
    helpMap["MM_PER_PIXEL"] = "Height/width of a pixel on the detector in mm. Default 0.11 mm to match CSPAD detector. (MPCCD is 0.05 if you need this value.)";
    helpMap["DETECTOR_SIZE"] = "Int x, y – width and height of the detector. Default 1765, 1765 to match CSPAD detector. Can be overridden by values in HDF5 file.";
    helpMap["OVER_PRED_BANDWIDTH"] = "Integration requires an over-estimated bandwidth to ‘catch’ all the reflections and more. This represents the % bandwidth to either side of the mean wavelength which should be used. Default 0.03.";
    helpMap["OVER_PRED_RLP_SIZE"] = "Reciprocal lattice point radius in Å-1 at the Miller index position (0, 0, 0) (i.e. without consideration of mosaicity). Default 2 × 10-4 Å-1. If you aren't picking up all your spots at low resolution, please use this.";
    helpMap["REFINE_ORIENTATIONS"] = "If set to ON, using the INTEGRATE command will also refine orientation matrices according to the protocol in Ginn et al Nat Comms 2015. Default ON.";
    helpMap["INDEXING_ORIENTATION_TOLERANCE"] = "Initial orientation matrix refinement will finish when the step search angle decreases below x degrees. Default 1 × 10-3º.";
    helpMap["INTENSITY_THRESHOLD"] = "Absolute (counts) or signal-to-noise (I/sigI) threshold at which a spot is considered to be a strong spot for initial orientation matrix refinement. Default 12 (assuming a detector gain of 1.0) and also dependent on ABSOLUTE_INTENSITY.";
    helpMap["ABSOLUTE_INTENSITY"] = "If set to ON, intensity threshold is measured in counts, whereas if set to OFF this is a I/sigI threshold for spots to be considered strong spots for initial orientation matrix refinement. Default OFF.";
    helpMap["METROLOGY_SEARCH_SIZE"] = "Metrology search size can be set to an integer number of pixels from which the predicted spot position may deviate in order to find the highest pixel value. Integration will be performed centred on this highest pixel value. Default 0 but needs changing.";
    helpMap["SHOEBOX_FOREGROUND_PADDING"] = "For a simple shoebox centred on a coordinate, this represents the extra number of pixels in both axes which will be integrated as part of the signal foreground. This overrides the background and neither flags. Default 1.";
    helpMap["SHOEBOX_NEITHER_PADDING"] = "For a simple shoebox centred on a coordinate, this represents the extra number of pixels from the centre in both axes for which pixels will not be considered as part of foreground or background signal. This overrides the background flags. Default 2.";
    helpMap["SHOEBOX_BACKGROUND_PADDING"] = "For a simple shoebox centred on a coordinate, this represents the extra number of pixels from the centre in both axes for which pixels will be used to estimate the background. Default 3.";
    helpMap["COMPLEX_SHOEBOX"] = "A complex shoebox will be generated to attempt to take into account the shape of the spot at high resolution, which may be distorted by the wavelength spread of the beam. May be buggy in combination with HDF5.";
    helpMap["MAX_INTEGRATED_RESOLUTION"] = "If set, integration will only occur up to the specified resolution. Default 1.4 Å.";
    helpMap["UNIT_CELL"] = "Unit cell parameters for the crystal of interest.";
    helpMap["FIX_UNIT_CELL"] = "If set to ON, the unit cell found in the ORIENTATION_MATRIX_LIST file will have its unit cell reset to the values specified in UNIT_CELL, preserving orientation, before any further integration or post-refinement.";
    helpMap["INITIAL_ORIENTATION_STEP"] = "For initial orientation matrix refinement, this is the angle in degrees that the algorithm steps in order to search for a new minimum. Step search algorithm is used. Default 1.0º.";
    helpMap["ORIENTATION_CORRECTION"] = "Sometimes the indexing solutions from other programs must be corrected by a systematic error in the horizontal and vertical rotations due to a component of the experiment not modelled by cppxfel. If this systematic shift is known, it can be supplied in degrees (x in the horizontal axis and y in the vertical axis). These values may be deduced from a table of corrections produced at the end of an INTEGRATE command. Default 0,0.";
    helpMap["IMAGE_MASKED_VALUE"] = "If there is a pixel value used as a masking value in the images provided, this can be set using this command. If IMAGE_MASKED_VALUE is not set, the default is not to mask any pixel value.";
    helpMap["SPHERE_THICKNESS"] = "This is the thickness in Å-1 around the mean Ewald radius to consider integrating. The subset of reflections which are chosen for integration are then determined by OVER_PRED_RLP_SIZE, so SPHERE_THICKNESS is a prior scanning step. Default 0.02 Å-1. Increasing this will improve spot prediction at low resolution but will reduce the overall speed of the program.";
    helpMap["PIXEL_COUNT_CUTOFF"] = "If, during integration, a pixel exceeds the value x, this reflection is rejected. Can be used to protect against non-linear gain of the detector at high pixel counts. Default is not to use.";
    helpMap["OPTIMISING_UNIT_CELL_A"] = "Sets whether unit cell dimension A is a refined parameter in either integration or post-refinement. Default OFF, maybe turn on for certain space groups. Temporarily not working!!!";
    helpMap["OPTIMISING_UNIT_CELL_B"] = "Sets whether unit cell dimension B is a refined parameter in either integration or post-refinement. Default OFF, maybe turn on for certain space groups. Temporarily not working!!!";
    helpMap["OPTIMISING_UNIT_CELL_C"] = "Sets whether unit cell dimension C is a refined parameter in either integration or post-refinement. Default OFF, maybe turn on for certain space groups. Temporarily not working!!!";
    helpMap["OPTIMISING_UNIT_CELL_ALPHA"] = "If set to ON, this unit cell dimension is refined in either integration or post-refinement, using the step size specified by STEP_SIZE_UNIT_CELL_ALPHA.  Temporarily not working!!!";
    helpMap["OPTIMISING_UNIT_CELL_BETA"] = "If set to ON, this unit cell dimension is refined in either integration or post-refinement, using the step size specified by STEP_SIZE_UNIT_CELL_BETA.  Temporarily not working!!!";
    helpMap["OPTIMISING_UNIT_CELL_GAMMA"] = "If set to ON, this unit cell dimension is refined in either integration or post-refinement, using the step size specified by STEP_SIZE_UNIT_CELL_GAMMA.  Temporarily not working!!!";
    helpMap["STEP_SIZE_UNIT_CELL_A"] = "this determines the step size used during either initial orientation matrix refinement or post-refinement for unit cell dimension A. Default 0.2 Å.";
    helpMap["STEP_SIZE_UNIT_CELL_B"] = "this determines the step size used during either initial orientation matrix refinement or post-refinement for unit cell dimension B. Default 0.2 Å.";
    helpMap["STEP_SIZE_UNIT_CELL_C"] = "this determines the step size used during either initial orientation matrix refinement or post-refinement for unit cell dimension C. Default 0.2 Å.";
    helpMap["STEP_SIZE_UNIT_CELL_ALPHA"] = "this determines the step size used during either initial orientation matrix refinement or post-refinement for unit cell angle alpha. Default 0.2º.";
    helpMap["STEP_SIZE_UNIT_CELL_BETA"] = "this determines the step size used during either initial orientation matrix refinement or post-refinement for unit cell angle beta. Default 0.2º.";
    helpMap["STEP_SIZE_UNIT_CELL_GAMMA"] = "this determines the step size used during either initial orientation matrix refinement or post-refinement for unit cell angle gamma. Default 0.2º.";
    helpMap["FROM_DIALS"] = "If you are using a matrix from DIALS, then this should be set to ON. Not applicable to MATRIX_LIST_VERSION 3.0.";
    helpMap["REJECT_OVERLAPPING_REFLECTIONS"] = "Usually, if there are multiple lattices on the same image then overlapping shoeboxes lead to reflections being rejected. This also applies to overlapping reflections from the same lattice, if the over-prediction is too wide. This may not always be a desired feature, so setting this value to OFF will disable it. Default ON.";
    helpMap["REFINE_IN_PLANE_OF_DETECTOR"] = "In some cases, the rotation of the orientation matrix around the beam axis is not exact. Setting this value to ON will refine the orientation matrix around this axis. Default ON. Unlikely to want to change this.";
    helpMap["FIT_BACKGROUND_AS_PLANE"] = "To fit the background pixels as a plane rather than simple summation, switch this to ON. Individual pixels may also be rejected as outliers. Default OFF. Not sure if currently functional.";
    helpMap["ROUGH_CALCULATION"] = "If switched to ON, some recalculations are skipped during integration, which vastly increases the speed of integration. Comes with a mild warning that the integration solutions may be compromised - in practice this appears to be fine in most cases, but may be preferable to switch off (will slow down significantly). Default ON.";
    
    helpMap["READ_REFINED_MTZS"] = "If reading from matrix version 3.0, READ_REFINED_MTZS set to ON can be used to force the refined images to be read in. Default OFF.";
    helpMap["IGNORE_MISSING_IMAGES"] = "If image data happens to be unavailable, but you have spot-finding results or other metadata, this will try to continue to run regardless. If it tries to load image data, this is horrendous. Default OFF.";
    helpMap["BINARY_PARTIALITY"] = "If a varying partiality between 0 and 1 is not working for you, maybe a BINARY_PARTIALITY would work better.";
    helpMap["DETECTOR_LIST"] = "Path to a file containing detector information. Should use one of the formats specified in GEOMETRY_FORMAT.";
    helpMap["GEOMETRY_FORMAT"] = "Read in with the detector information from cppxfel or CrystFEL format. Can also load in panel_list format for backwards compatibility.";
    helpMap["HDF5_SOURCE_FILES"] = "HDF5 files from which image data should be found. Supports glob strings (e.g. run*.h5 for SACLA HDF5 files from cheetah-dispatcher).";
    helpMap["USE_HDF5_WAVELENGTH"] = "Use the wavelengths stored within HDF5 files. If you do not trust these values, disable this in order to default to the value of INTEGRATION_WAVELENGTH.";
    helpMap["HDF5_AS_FLOAT"] = "HDF5 file image data should be interpreted as float (rudimentary conversion to ints!)";
    
    helpMap["FREE_ELECTRON_LASER"] = "Which free electron laser did this data come from? This is used for interpreting HDF5 files. Only LCLS and SACLA currently supported.";
    helpMap["OUTPUT_DIRECTORY"] = "Path to a directory into which almost all processing files will be deposited, except for spot-finding results.";
    
    helpMap["SPOT_FINDING_ALGORITHM"] = "Choose which algorithm is used for spot-finding. Blob-finding is highly recommended.";
    helpMap["IMAGE_MIN_SPOT_INTENSITY"] = "Do not attempt to perform spot-finding calculation on a pixel intensity below this value.";
    helpMap["IMAGE_MIN_CORRELATION"] = "Minimum correlation between 'model' spot and 'real' spot below which a spot will be rejected.";
    helpMap["IMAGE_SPOT_PROBE_PADDING"] = "Model spot will have a signal area of dimensions (2 * IMAGE_SPOT_PROBE_PADDING + 1).";
    helpMap["IMAGE_SPOT_PROBE_BG_PADDING"] = "Model spot will have a background area of this pixel width padding around the signal region.";
    helpMap["IMAGE_SPOT_PROBE_HEIGHT"] = "Height of signal in model spot. Having a lower value will make the correlation more sensitive, I think. Default 100.";
    helpMap["FORCE_SPOT_FINDING"] = "If set to ON, spot-finding will be re-attempted on each run of cppxfel even if spots have already been saved to disk. Default OFF (reload from disk).";

    helpMap["SOLUTION_ATTEMPTS"] = "Maximum number of lattices which should be attempted to be indexed by the TakeTwo algorithm before stopping. Default 1.";
    helpMap["INDEXING_TIME_LIMIT"] = "Maximum number of seconds after which cppxfel will give up on indexing a lattice.";
    helpMap["MAX_RECIPROCAL_DISTANCE"] = "Maximum distance between two potential reciprocal lattice points used for TakeTwo indexing.";
    helpMap["REJECT_UNDER_SPOT_COUNT"] = "Do not index or use this image for generating powder patterns if underneath this spot count. Default 0.";
    helpMap["REJECT_OVER_SPOT_COUNT"] = "Do not index or use this image for generating powder patterns if over this spot count. Default 4000.";

    helpMap["PNG_TOTAL"] = "Specifies how many PNG files should be written to disk by cppxfel when WRITE_PNGS is called in the commands section.";
    helpMap["PNG_SHOEBOX"] = "If writing PNGs when WRITE_PNGS is called, whether to draw the integration shoeboxes (currently might be buggy for complex/non-symmetrical shoeboxes).";
    helpMap["IMAGE_IGNORE_UNDER_VALUE"] = "If a value on an image goes below this value, it is masked out.";

    
    helpMap["SWEEP_DETECTOR_DISTANCE"] = "If you wish to refine geometry, specifying two distances in mm for this parameter will perform a sweep of intervening detector distances before continuing with geometry refinement. At the moment it will refine against the pseudo-powder pattern. Not completely tested, default is not to be set.";
    helpMap["EXPECTED_GEOMETRY_MOVEMENT"] = "Default increment for geometry refinement. The default is 0.02, but if you suspect convergence is too slow or you are driving into a local minimum too quickly, try increasing or decreasing this value respectively. Helen is not yet aware if it is worth changing this value.";
    
    
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
    //parserMap["MERGE_MEDIAN"] = simpleBool;
    parserMap["READ_REFINED_MTZS"] = simpleBool;
    
    parserMap["SET_SIGMA_TO_UNITY"] = simpleBool;
    parserMap["APPLY_UNREFINED_PARTIALITY"] = simpleBool;
    parserMap["BINARY_PARTIALITY"] = simpleBool;
    parserMap["MINIMIZATION_METHOD"] = simpleInt;
    parserMap["NELDER_MEAD_CYCLES"] = simpleInt;
    parserMap["MEDIAN_WAVELENGTH"] = simpleBool;
    parserMap["WAVELENGTH_RANGE"] = doubleVector;
    parserMap["ALLOW_TRUST"] = simpleBool;
    parserMap["PARTIALITY_CUTOFF"] = simpleFloat;
	parserMap["DEFAULT_TARGET_FUNCTION"] = simpleInt;
    parserMap["TARGET_FUNCTIONS"] = intVector;
   // parserMap["RLP_MODEL"] = simpleInt;
	parserMap["CORRELATION_THRESHOLD"] = simpleFloat;
    parserMap["PARTIALITY_CORRELATION_THRESHOLD"] = simpleFloat;
	parserMap["MAX_REFINED_RESOLUTION"] = simpleFloat;
    parserMap["MERGE_TO_RESOLUTION"] = simpleFloat;
    parserMap["MIN_REFINED_RESOLUTION"] = simpleFloat; // simplify all these resolutions?
    parserMap["INITIAL_CORRELATION_THRESHOLD"] = simpleFloat;
	parserMap["THRESHOLD_SWAP"] = simpleInt;
	parserMap["OUTLIER_REJECTION_SIGMA"] = simpleFloat;
	parserMap["OUTLIER_REJECTION"] = simpleBool;
	parserMap["CORRELATION_REJECTION"] = simpleBool;
	parserMap["PARTIALITY_REJECTION"] = simpleBool;
    parserMap["POLARISATION_CORRECTION"] = simpleBool;
    parserMap["POLARISATION_FACTOR"] = simpleFloat;
    parserMap["REFINEMENT_INTENSITY_THRESHOLD"] = simpleFloat; // merge with intensity threshold?
    parserMap["TRUST_INDEXING_SOLUTION"] = simpleBool;
    parserMap["CUSTOM_AMBIGUITY"] = doubleVector;
    parserMap["R_FACTOR_THRESHOLD"] = simpleFloat;
    parserMap["REINITIALISE_WAVELENGTH"] = simpleBool;
    parserMap["SMOOTH_FUNCTION"] = simpleBool;
    parserMap["NORMALISE_PARTIALITIES"] = simpleBool;
    parserMap["REPLACE_REFERENCE"] = simpleBool;
    parserMap["FREE_MILLER_LIST"] = simpleString;
    parserMap["FREE_MILLER_PROPORTION"] = simpleFloat;
    parserMap["USE_HDF5_WAVELENGTH"] = simpleBool;
    parserMap["FREE_ELECTRON_LASER"] = simpleInt;
    parserMap["CHEETAH_DATA_ADDRESSES"] = simpleString;
    parserMap["CHEETAH_ID_ADDRESSES"] = simpleString;
    parserMap["HDF5_MASK_ADDRESS"] = simpleString;
    parserMap["HDF5_AS_FLOAT"] = simpleBool;
    
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
    parserMap["OPTIMISING_UNIT_CELL_ALPHA"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_BETA"] = simpleBool;
    parserMap["OPTIMISING_UNIT_CELL_GAMMA"] = simpleBool;

    parserMap["HDF5_SOURCE_FILES"] = stringVector;
 //   parserMap["HDF5_OUTPUT_FILE"] = simpleString;
    parserMap["DUMP_IMAGES"] = simpleBool;
    
    parserMap["ORIENTATION_MATRIX_LIST"] = simpleString;
    parserMap["OUTPUT_DIRECTORY"] = simpleString;
    parserMap["OUTPUT_INDIVIDUAL_CYCLES"] = simpleBool;
    parserMap["MATRIX_LIST_VERSION"] = simpleFloat;
	parserMap["INITIAL_MTZ"] = simpleString;
    parserMap["IMAGE_LIMIT"] = simpleInt;
    parserMap["IMAGE_SKIP"] = simpleInt;
    parserMap["NEW_MATRIX_LIST"] = simpleString;
    
    parserMap["RECALCULATE_WAVELENGTHS"] = simpleBool;
    parserMap["MERGE_ANOMALOUS"] = simpleBool;
    parserMap["SCALING_STRATEGY"] = simpleInt;
	parserMap["MINIMUM_REFLECTION_CUTOFF"] = simpleInt;
    parserMap["MINIMUM_MULTIPLICITY"] = simpleInt;
    parserMap["REJECT_BELOW_SCALE"] = simpleFloat;

	// Indexing parameters

 //   parserMap["DETECTOR_GAIN"] = simpleFloat;
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
  	parserMap["INTENSITY_THRESHOLD"] = simpleFloat;
    parserMap["ABSOLUTE_INTENSITY"] = simpleBool;
	parserMap["METROLOGY_SEARCH_SIZE"] = simpleInt;
    parserMap["METROLOGY_SEARCH_SIZE_BIG"] = simpleInt;
    parserMap["SHOEBOX_FOREGROUND_PADDING"] = simpleInt;
	parserMap["SHOEBOX_NEITHER_PADDING"] = simpleInt;
	parserMap["SHOEBOX_BACKGROUND_PADDING"] = simpleInt;
    parserMap["SHOEBOX_MAKE_EVEN"] = simpleBool;
    parserMap["MIN_INTEGRATED_RESOLUTION"] = simpleFloat;
    parserMap["MAX_INTEGRATED_RESOLUTION"] = simpleFloat;
	parserMap["UNIT_CELL"] = doubleVector;
    parserMap["FIX_UNIT_CELL"] = simpleBool;
    parserMap["INITIAL_ORIENTATION_STEP"] = simpleFloat;
    parserMap["COMPLEX_SHOEBOX"] = simpleBool;
    parserMap["SHOEBOX_BANDWIDTH_MULTIPLIER"] = simpleFloat;
    parserMap["PIXEL_LEAK"] = simpleFloat;
    parserMap["MILLER_INDIVIDUAL_WAVELENGTHS"] = simpleBool;
    parserMap["UNBALANCED_REFLECTIONS"] = simpleBool;
    parserMap["ORIENTATION_CORRECTION"] = doubleVector;
    parserMap["IMAGE_MASKED_VALUE"] = simpleInt;
    parserMap["IMAGE_IGNORE_UNDER_VALUE"] = simpleInt;
    parserMap["SPHERE_THICKNESS"] = simpleFloat;
    parserMap["PIXEL_COUNT_CUTOFF"] = simpleInt;
    parserMap["BAD_PIXELS"] = simpleString;
    parserMap["STEP_SIZE_UNIT_CELL_A"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_B"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_C"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_ALPHA"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_BETA"] = simpleFloat;
    parserMap["STEP_SIZE_UNIT_CELL_GAMMA"] = simpleFloat;
    parserMap["FROM_DIALS"] = simpleBool;
    parserMap["REJECT_OVERLAPPING_REFLECTIONS"] = simpleBool;
    parserMap["REFINE_IN_PLANE_OF_DETECTOR"] = simpleBool;
    parserMap["FIT_BACKGROUND_AS_PLANE"] = simpleBool;
    parserMap["ROUGH_CALCULATION"] = simpleBool;
    
    parserMap["MINIMUM_SPOTS_EXPLAINED"] = simpleInt;
    parserMap["MINIMUM_TRUST_ANGLE"] = simpleFloat;
    parserMap["SOLUTION_ANGLE_SPREAD"] = simpleFloat;
    parserMap["REJECT_CLOSE_SPOTS"] = simpleBool;
    parserMap["ACCEPT_ALL_SOLUTIONS"] = simpleBool;
    parserMap["INDEXING_MIN_RESOLUTION"] = simpleFloat;
    parserMap["SPOTS_PER_LATTICE"] = simpleInt;
    parserMap["GOOD_SOLUTION_ST_DEV"] = simpleFloat;
    parserMap["GOOD_SOLUTION_SUM_RATIO"] = simpleFloat;
    parserMap["GOOD_SOLUTION_HIGHEST_PEAK"] = simpleInt;
	parserMap["BAD_SOLUTION_ST_DEV"] = simpleFloat;
	parserMap["BAD_SOLUTION_HIGHEST_PEAK"] = simpleInt;
	parserMap["SOLUTION_ATTEMPTS"] = simpleInt;
    parserMap["MAX_RECIPROCAL_DISTANCE"] = simpleFloat;
    parserMap["MAXIMUM_ANGLE_DISTANCE"] = simpleFloat;
    parserMap["ALWAYS_FILTER_SPOTS"] = simpleBool;
    parserMap["MINIMUM_SOLUTION_NETWORK_COUNT"] = simpleInt;
    parserMap["NETWORK_TRIAL_LIMIT"] = simpleInt;
    parserMap["INDEXING_TIME_LIMIT"] = simpleInt;
    parserMap["MAX_LATTICES_PER_IMAGE"] = simpleInt;
    parserMap["CHECKING_COMMON_SPOTS"] = simpleBool;
    parserMap["EXCLUDE_WEAKEST_SPOT_FRACTION"] = simpleFloat;
    parserMap["REJECT_UNDER_SPOT_COUNT"] = simpleInt;
    parserMap["REJECT_OVER_SPOT_COUNT"] = simpleInt;
    parserMap["POWDER_PATTERN_STEP"] = simpleFloat;
    parserMap["POWDER_PATTERN_STEP_ANGLE"] = simpleFloat;
    parserMap["LOW_MEMORY_MODE"] = simpleBool;
    parserMap["GEOMETRY_FORMAT"] = simpleInt;
    
    // force spot finding
    
    parserMap["IMAGE_MIN_SPOT_INTENSITY"] = simpleFloat;
    parserMap["IMAGE_MIN_CORRELATION"] = simpleFloat;
    parserMap["IMAGE_ACCEPTABLE_COMBO"] = doubleVector;
    parserMap["IMAGE_SPOT_PROBE_HEIGHT"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_BACKGROUND"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_PADDING"] = simpleInt;
    parserMap["IMAGE_SPOT_PROBE_BG_PADDING"] = simpleInt;
    parserMap["PROBE_DISTANCES"] = doubleVector;
    parserMap["FORCE_SPOT_FINDING"] = simpleBool;
 //  parserMap["FORCE_RECENTRE_SPOTS"] = simpleBool;
 //   parserMap["FORCE_RESTART_POST_REFINEMENT"] = simpleBool;
    parserMap["SPOT_FINDING_MIN_PIXELS"] = simpleInt;
    parserMap["SPOT_FINDING_SIGNAL_TO_NOISE"] = simpleFloat;
    parserMap["SPOT_FINDING_MAX_PIXELS"] = simpleInt;
    parserMap["SPOT_FINDING_ALGORITHM"] = simpleInt;
 //   parserMap["SPOTS_ARE_RECIPROCAL_COORDINATES"] = simpleBool;
    
    parserMap["IGNORE_MISSING_IMAGES"] = simpleBool;
    parserMap["FROM_TIFFS"] = simpleBool;
    
//    parserMap["RANDOM_SEED"] = simpleInt;
    
    parserMap["PNG_TOTAL"] = simpleInt;
    parserMap["PNG_THRESHOLD"] = simpleFloat;
    parserMap["PNG_SHOEBOX"] = simpleBool;
    parserMap["PNG_STRONG_ONLY"] = simpleBool;
    parserMap["PNG_ALL_LATTICES"] = simpleBool;
    parserMap["PNG_HEIGHT"] = simpleInt;
    parserMap["DRAW_GEOMETRY_PNGS"] = simpleBool;
    parserMap["SWEEP_DETECTOR_DISTANCE"] = doubleVector;
    parserMap["ENABLE_IMAGE_CSVS"] = simpleBool;
    
    parserMap["TRUST_GLOBAL_GEOMETRY"] = simpleBool;
    parserMap["TRUST_QUADRANT_GEOMETRY"] = simpleBool;
    parserMap["TRUST_LOCAL_GEOMETRY"] = simpleBool;
    parserMap["EXPECTED_GEOMETRY_MOVEMENT"] = simpleFloat;
    
    parserMap["DETECTOR_LIST"] = simpleString;

    
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

    if (deprecatedList.count(upperCommand) > 0)
    {
        logged << "Deprecated command: " << command << std::endl;
        logged << deprecatedList[command] << std::endl;
        sendLogAndExit();
        return NULL;
    }
    
    // make a list of permitted words soon
    else if (parserMap.count(upperCommand) == 0 && upperCommand != "PANEL" && upperCommand != "MASK" && upperCommand != "SOLVENT_MASK")
    {
        logged << "Error: do not understand command " << upperCommand << std::endl;
        sendLogAndExit();
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
		logged << "Warning: " << line << " has no assignment" << std::endl;
        sendLog();
		return false;
	}

	return true;
}

FileParser::FileParser(void)
{
    generateCategories();
    generateFunctionList();
    generateHelpList();
    generateCodeList();
    generateDeprecatedList();
}

FileParser::FileParser(std::string name, std::vector<std::string> someExtras)
{
	this->filename = name;
	generateFunctionList();
    generateCodeList();
    generateDeprecatedList();
    
    extras = someExtras;
}

FileParser::~FileParser()
{
	// TODO Auto-generated destructor stub
}

