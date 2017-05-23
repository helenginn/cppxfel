/*
 * InputFileParser.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#include "InputFileParser.h"
#include "FileReader.h"
#include "MtzRefiner.h"
#include "Miller.h"
#include <sstream>
#include "Hdf5ManagerProcessing.h"
#include "Logger.h"
//#include <boost/python.hpp>

InputFileParser::InputFileParser(std::string filename, std::vector<std::string> someExtras) : FileParser(filename, someExtras)
{
}

InputFileParser::~InputFileParser()
{
	// TODO Auto-generated destructor stub
}

void InputFileParser::refine(int maxCycles)
{
    setKey("MAXIMUM_CYCLES", maxCycles);
    refiner->refine();
}

void InputFileParser::parseFromPython()
{
    parse(true);
}

void InputFileParser::integrate()
{
    if (refiner == NULL)
        parse(true);
    
    refiner->integrate();
}

int InputFileParser::processOptions(std::vector<std::string> lines)
{
    int continueFrom = -1;

    for (int i = 0; i < lines.size(); i++)
    {
        std::string line = lines[i];
        
        if (line.length() == 0)
            continue;
        
        if (line.substr(0, 1) == "#")
            continue;
        
        if (line == "COMMANDS")
        {
            continueFrom = i;
            break;
        }
        
        if (!checkSpaces(line))
            continue;
        
        if (line.length() == 0)
            continue;
        
        std::string command; std::string rest;
        
        ParserFunction function = this->splitLine(line, command, rest);
        
        if (parserMap.count(command) == 0)
        {
            logged << "Warning: Line \"" << line << "\" not recognised."
            << std::endl;
            sendLogAndExit();
        }
        
        if (function)
            function(&parameters, command, rest);
    }

    return continueFrom;
}

void InputFileParser::parse(bool fromPython)
{
 //   Logger::mainLogger = LoggerPtr(new Logger());
 //   boost::thread thr = boost::thread(Logger::awaitPrintingWrapper, Logger::mainLogger);

    parameters = ParametersMap();

    if (!FileReader::exists(filename))
    {
        logged << "Input file " << filename << " seems to be missing." << std::endl;
        sendLogAndExit();
    }
    
	std::string fileContents = FileReader::get_file_contents(filename.c_str());
	vector<std::string> fileLines = FileReader::split(fileContents, '\n');

	if (fromPython)
    {
  //      log << "If you are running initial orientation matrix refinement using image files, please load individual experiments from DIALS indexing results. If you are running post-refinement on individual MTZ files, please supply the appropriate orientation matrix list file on the keyword ORIENTATION_MATRIX_LIST. This command will be ignored for image files. If you want to load image files please do so from your Python script." << std::endl;
        
        log << "Imported commands from parameter file. Waiting for commands." << std::endl;
    }
    
    Logger::mainLogger->addStream(&log);
    log.str("");
    
    if (fromPython)
    {
        refiner->setFromPython(true);
        return;
    }
    
    splitCharMajor = ' ';
    splitCharMinor = ' ';
    
    int continueFrom = processOptions(fileLines);
    
    splitCharMajor = '=';
    splitCharMinor = ',';
    
    log << " --- Overwriting commands from command line ---" << std::endl;
    processOptions(extras);
    log << " --- Finished command line options ---" << std::endl;

    Logger::mainLogger->addStream(&log);
    log.str("");
    
    refiner = boost::shared_ptr<MtzRefiner>(new MtzRefiner());
    Miller::setupStaticVariables();
    Hdf5ManagerProcessing::setupProcessingManager();
    
    int seed = FileParser::getKey("RANDOM_SEED", 0);
    srand((unsigned int)seed);
    
	if (continueFrom != -1)
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
                refiner->integrate();
			}
            
            if (line == "INTEGRATE_SPOTS")
            {
                understood = true;
                refiner->integrateSpots();
            }
            
			if (line == "REFINE_PARTIALITY" || line == "POST_REFINE")
			{
                understood = true;
                refiner->refine();
			}
            
            if (line == "REPORT_GEOMETRY")
            {
                understood = true;
                refiner->reportMetrology();
            }
            
            if (line == "REFINE_METROLOGY" || line == "REFINE_GEOMETRY")
            {
                understood = true;
                refiner->refineMetrology(false);
            }

			if (line == "REFINE_UNIT_CELL")
			{
				understood = true;
				refiner->refineUnitCell();
			}


			if (line == "MERGE")
			{
                understood = true;
                refiner->merge(true);
			}
            
            if (line == "DISPLAY_INDEXING_HANDS")
            {
                understood = true;
                refiner->displayIndexingHands();
            }
            
            if (line == "CORRELATION_PLOT")
            {
                understood = true;
                refiner->correlationAndInverse();
            }
            
            if (line == "LOAD_MTZ_FILES")
            {
                understood = true;
                refiner->loadPanels(false);
                refiner->readMatricesAndMtzs();
            }
            
            if (line == "LOAD_IMAGE_FILES")
            {
                understood = true;
                refiner->loadImageFiles();
            }
            
            if (line == "LOAD_INITIAL_MTZ")
            {
                understood = true;
                refiner->loadInitialMtz(true);
            }
            
            if (line == "LINEAR_SCALING")
            {
                understood = true;
                refiner->linearScaling();
            }
            
            if (line == "WRITE_NEW_MATRICES")
            {
                understood = true;
                refiner->writeNewOrientations(true, true);
            }
            
            if (line == "INITIAL_MERGE")
            {
                understood = true;
                refiner->initialMerge();
            }

			if (line == "PIXEL_VS_FIDUCIAL")
			{
				understood = true;
				refiner->plotPixelValueVsFiducial();
			}

            if (line == "ORIENTATION_PLOT")
            {
                understood = true;
                refiner->orientationPlot();
            }
            
            if (line == "APPLY_PARTIALITY")
            {
                understood = true;
                refiner->applyUnrefinedPartiality();
            }
            
            if (line == "INDEX")
            {
                understood = true;
                refiner->index();
            }
            
            if (line == "MAXIMUM_IMAGE")
            {
                understood = true;
                refiner->maximumImage();
            }
            
            if (line == "POWDER_PATTERN")
            {
                understood = true;
                refiner->powderPattern();
            }
            
            if (line == "WRITE_PNGS")
            {
                understood = true;
                refiner->writePNGs();
            }
            
            if (line == "IMAGE_TO_DETECTOR_MAP")
            {
                understood = true;
                refiner->imageToDetectorMap();
            }
            
            if (line == "PLOT_INTENSITIES")
            {
                understood = true;
                refiner->plotIntensities();
            }
            
            if (line == "PLOT_INTEGRATION_WINDOWS")
            {
                understood = true;
                refiner->plotIntegrationWindows();
            }
            
            if (line == "FAKE_SPOTS")
            {
                understood = true;
                refiner->fakeSpots();
            }
            
            if (line == "TAKETWO_PNG")
            {
                understood = true;
                refiner->takeTwoPNG();
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
		refiner->refine();
	}
}

