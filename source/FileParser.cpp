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
#include "MtzRefiner.h"

ParametersMap FileParser::parameters;
std::ostringstream FileParser::log;
int FileParser::threadsFound = 0;

int FileParser::getMaxThreads()
{
    if (threadsFound != 0)
        return threadsFound;
    
    ostringstream logged;
    
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

void FileParser::doubleVector(ParametersMap *map, std::string command,
		std::string rest)
{
	std::vector<std::string> components = FileReader::split(rest, ' ');
	std::vector<double> doubleVector;

	log << "Setting " << command << " to ";

	for (int i = 0; i < components.size(); i++)
	{
		double theFloat = atof(components[i].c_str());
		log << theFloat << " ";
		doubleVector.push_back(theFloat);
	}

	log << std::endl;

	(*map)[command] = doubleVector;
}

void FileParser::intVector(ParametersMap *map, std::string command,
		std::string rest)
{
	std::vector<std::string> components = FileReader::split(rest, ' ');
	std::vector<int> intVector;

	log << "Setting " << command << " to ";

	for (int i = 0; i < components.size(); i++)
	{
		int theInt = atoi(components[i].c_str());
		log << theInt << " ";
		intVector.push_back(theInt);
	}

	log << std::endl;

	(*map)[command] = intVector;
}

void FileParser::generateFunctionList()
{
	parserMap = ParserMap();

    parserMap["VERBOSITY_LEVEL"] = simpleInt;
    parserMap["MAX_THREADS"] = simpleInt;
    
	// Refinement parameters
	parserMap["REMOVE_WEDGE"] = simpleFloat;

    parserMap["MINIMUM_CYCLES"] = simpleInt;
	parserMap["STOP_REFINEMENT"] = simpleBool;

    parserMap["ALLOW_TRUST"] = simpleBool;
    parserMap["EXCLUDE_OWN_REFLECTIONS"] = simpleBool;
    parserMap["PARTIALITY_CUTOFF"] = simpleFloat;
	parserMap["DEFAULT_TARGET_FUNCTION"] = simpleInt;
	parserMap["USE_PARTIALITY_FUNCTION"] = simpleBool;
	parserMap["CORRELATION_THRESHOLD"] = simpleFloat;
	parserMap["MAX_RESOLUTION_ALL"] = simpleFloat;
	parserMap["MAX_RESOLUTION_RLP_SIZE"] = simpleFloat;
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
    parserMap["REFINE_B_FACTOR"] = simpleBool;
    parserMap["INITIAL_GRID_SEARCH"] = simpleBool;
    
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

	parserMap["ORIENTATION_MATRIX_LIST"] = simpleString;
    parserMap["MATRIX_LIST_VERSION"] = simpleFloat;
	parserMap["INITIAL_MTZ"] = simpleString;
    parserMap["IMAGE_LIMIT"] = simpleInt;

    parserMap["RECALCULATE_SIGMA"] = simpleBool;
    parserMap["MERGE_ANOMALOUS"] = simpleBool;
    parserMap["FAKE_ANOMALOUS"] = simpleBool;
    parserMap["SCALING_STRATEGY"] = simpleInt;
	parserMap["MINIMUM_REFLECTION_CUTOFF"] = simpleInt;
    parserMap["APPLY_INFLATION"] = simpleBool;

	// Indexing parameters

	parserMap["SPACE_GROUP"] = simpleInt;
	parserMap["INTEGRATION_WAVELENGTH"] = simpleFloat;
	parserMap["DETECTOR_DISTANCE"] = simpleFloat;
    parserMap["BEAM_CENTRE"] = doubleVector;
    parserMap["MM_PER_PIXEL"] = simpleFloat;
	parserMap["OVER_PRED_BANDWIDTH"] = simpleFloat;
	parserMap["OVER_PRED_RLP_SIZE"] = simpleFloat;
    parserMap["REFINE_ORIENTATIONS"] = simpleBool;
    parserMap["REFINE_DISTANCES"] = simpleBool;
    parserMap["INDEXING_ORIENTATION_TOLERANCE"] = simpleFloat;
	parserMap["INTENSITY_THRESHOLD"] = simpleFloat;
    parserMap["ABSOLUTE_INTENSITY"] = simpleBool;
	parserMap["METROLOGY_SEARCH_SIZE"] = simpleInt;
	parserMap["SHOEBOX_FOREGROUND_RADIUS"] = simpleInt;
	parserMap["SHOEBOX_NEITHER_RADIUS"] = simpleInt;
	parserMap["SHOEBOX_BACKGROUND_RADIUS"] = simpleInt;
    parserMap["COMPLEX_SHOEBOX"] = simpleBool;
    parserMap["MAX_INTEGRATED_RESOLUTION"] = simpleFloat;
	parserMap["UNIT_CELL"] = doubleVector;
    parserMap["FIX_UNIT_CELL"] = simpleBool;
    parserMap["ADD_MASK"] = intVector;
    parserMap["INITIAL_ORIENTATION_STEP"] = simpleFloat;
    parserMap["SHOEBOX_BANDWIDTH_MULTIPLIER"] = simpleFloat;
    parserMap["PIXEL_LEAK"] = simpleFloat;
    parserMap["ORIENTATION_SCORE"] = simpleInt;
    parserMap["IMAGE_MASKED_VALUE"] = simpleInt;
    parserMap["SPHERE_THICKNESS"] = simpleFloat;
    parserMap["SIGMA_RESOLUTION_CUTOFF"] = simpleFloat;
    
	parserMap["PANEL_LIST"] = simpleString;
    parserMap["SKIP_LINES"] = simpleInt;
}

ParserFunction FileParser::splitLine(std::string line, std::string &command,
		std::string &rest)
{
	int space_index = (int)line.find_first_of(" ");

	command = line.substr(0, space_index);

	std::ostringstream stream;

	std::locale loc;
	for (std::string::size_type j = 0; j < command.length(); ++j)
		stream << std::toupper(command[j], loc);

	std::string upperCommand = stream.str();

	rest = line.substr(space_index + 1, std::string::npos);


	ParserFunction function = parserMap[upperCommand];
	command = upperCommand;

	return function;
}

bool FileParser::checkSpaces(std::string line)
{
	int space_index = (int)line.find_first_of(" ");

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

FileParser::FileParser(std::string name)
{
	std::cout << "Initialising parser" << std::endl;
    
	this->filename = name;
	generateFunctionList();
}

FileParser::~FileParser()
{
	// TODO Auto-generated destructor stub
}

