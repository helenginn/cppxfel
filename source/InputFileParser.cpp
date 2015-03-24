/*
 * InputFileParser.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#include "InputFileParser.h"
#include "FileReader.h"
#include "MtzRefiner.h"
#include <sstream>
#include "Logger.h"

InputFileParser::InputFileParser(std::string filename) : FileParser(filename)
{
	// TODO Auto-generated constructor stub
}

InputFileParser::~InputFileParser()
{
	// TODO Auto-generated destructor stub
}

void InputFileParser::parse()
{
	parameters = ParametersMap();

	std::string fileContents = FileReader::get_file_contents(filename.c_str());
	std::vector<std::string> fileLines = FileReader::split(fileContents, '\n');

	bool foundCommands = false;
	int continueFrom = 0;

	for (int i = 0; i < fileLines.size(); i++)
	{
		std::string line = fileLines[i];

		if (line.length() == 0)
			continue;

		if (line.substr(0, 1) == "#")
			continue;

		if (line == "COMMANDS")
		{
			log << "Found COMMANDS section" << std::endl;
			continueFrom = i;
			foundCommands = true;
			break;
		}

		if (!checkSpaces(line))
			continue;

		std::string command; std::string rest;
		ParserFunction function = this->splitLine(line, command, rest);

		if (parserMap.count(command) == 0)
		{
            std::cout << "Warning: Line \"" << line << "\" not recognised."
					<< std::endl;
			exit(0);
		}
        
		function(&parameters, command, rest);
	}
    
    Logger::mainLogger->addStream(&log);
    log.str("");

	MtzRefiner refiner = MtzRefiner();

	if (foundCommands)
	{
		for (int i = continueFrom; i < fileLines.size(); i++)
		{
			std::string line = fileLines[i];
            bool understood = false;
            
			if (line.length() == 0)
				continue;

			if (line == "INTEGRATE")
			{
                understood = true;
                refiner.integrate(false);
			}
            
            if (line == "REFINE_DETECTOR_GEOMETRY")
            {
                understood = true;
                refiner.refineDetectorGeometry();
            }

			if (line == "REFINE_PARTIALITY")
			{
                understood = true;
                refiner.refine();
			}
            
            if (line == "REFINE_METROLOGY")
            {
                understood = true;
                refiner.refineMetrology();
            }

			if (line == "MERGE")
			{
                understood = true;
                refiner.merge();
			}
            
            if (line == "LOAD_X_FILES")
            {
                understood = true;
                refiner.xFiles();
            }
            
            if (line == "REMOVE_SIGMA_VALUES")
            {
                understood = true;
                refiner.removeSigmaValues();
            }
            
            if (line == "DISPLAY_INDEXING_HANDS")
            {
                understood = true;
                refiner.displayIndexingHands();
            }
            
            if (line == "CORRELATION_PLOT")
            {
                understood = true;
                refiner.correlationAndInverse();
            }
            
            if (line == "POLARISATION_GRAPH")
            {
                understood = true;
                refiner.polarisationGraph();
            }
            
            if (line == "LOAD_MTZ_FILES")
            {
                understood = true;
                refiner.readMatricesAndMtzs();
            }
            
            if (line == "DETECTOR_GAINS")
            {
                understood = true;
                refiner.plotDetectorGains();
            }
            
            if (line == "LOAD_INITIAL_MTZ")
            {
                understood = true;
                refiner.loadInitialMtz(true);
            }
            
            if (line == "REFINE_WITH_SYMMETRY")
            {
                understood = true;
                refiner.refineSymmetry();
            }
            
            if (!understood)
            {
                log << "Skipping line " << line << std::endl;
            }
            
            if (understood)
            {
                log << "Executed line " << line << std::endl;
            }
            
            Logger::mainLogger->addStream(&log);
            log.str("");
		}
	}
	else
	{
		log << "No commands issued; defaulting to REFINE_PARTIALITY."
				<< std::endl;
        Logger::mainLogger->addStream(&log);
        log.str("");
		refiner.refine();
	}
}
