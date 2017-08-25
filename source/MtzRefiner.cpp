/*
 * MtzRefiner.cpp
 *
 *  Created on: 4 Jan 2015
 *      Author: helenginn
 */

#include "MtzRefiner.h"
#include <vector>
#include "misc.h"
#include "Vector.h"
#include <boost/thread/thread.hpp>
#include "GraphDrawer.h"
#include "Hdf5Image.h"
#include "Image.h"
#include "Miller.h"
#include <fstream>
#include "AmbiguityBreaker.h"
#include "IndexManager.h"
#include "UnitCellLattice.h"
#include "CSV.h"

#include "FileParser.h"
#include "parameters.h"
#include "FileReader.h"

#include "Logger.h"
#include "MtzMerger.h"
#include "Hdf5ManagerProcessing.h"
#include "Detector.h"
#include "GeometryParser.h"
#include "GeometryRefiner.h"

int MtzRefiner::imageLimit;
int MtzRefiner::cycleNum;

MtzRefiner::MtzRefiner()
{
    // TODO Auto-generated constructor stub
    int logInt = FileParser::getKey("VERBOSITY_LEVEL", 0);
    Logger::mainLogger->changePriorityLevel((LogLevel)logInt);
    
    imageLimit = FileParser::getKey("IMAGE_LIMIT", 0);
    
    hasRefined = false;
    isPython = false;
    readRefinedMtzs = FileParser::getKey("READ_REFINED_MTZS", false);
    
    indexManager = NULL;
}

// MARK: Refinement

void MtzRefiner::cycleThreadWrapper(MtzRefiner *object, int offset)
{
    object->cycleThread(offset);
}

void MtzRefiner::cycleThread(int offset)
{
    int img_num = (int)images.size();

    std::vector<int> targets = FileParser::getKey("TARGET_FUNCTIONS", std::vector<int>());

    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < img_num; i += maxThreads)
    {
        std::ostringstream logged;

		ImagePtr image = images[i];

		for (int j = 0; j < image->mtzCount(); j++)
		{
			MtzPtr mtz = image->mtz(j);

			if (!mtz->isRejected())
			{
				bool silent = (targets.size() > 0);

				mtz->refinePartialities();

				if (targets.size() > 0)
				{
					ScoreType firstScore = mtz->getScoreType();
					for (int i = 0; i < targets.size(); i++)
					{
						silent = (i < targets.size() - 1);
						mtz->setDefaultScoreType((ScoreType)targets[i]);

						mtz->refinePartialities();
					}

					mtz->setDefaultScoreType(firstScore);
				}

				mtz->setRefineOrientations(false);
			}
		}
    }
}

void MtzRefiner::cycle()
{
    MtzManager::setReference(&*reference);

    time_t startcputime;
    time(&startcputime);

    boost::thread_group threads;
    
    int maxThreads = FileParser::getMaxThreads();
    
    std::ostringstream logged;
    logged << "Filename\tScore type\t\tCorrel\tRfactor\tPart correl\tHits" << std::endl;
    Logger::mainLogger->addStream(&logged);
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(cycleThreadWrapper, this, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    time_t endcputime;
    time(&endcputime);
    
    clock_t difference = endcputime - startcputime;
    double seconds = difference;
    
    int finalSeconds = (int) seconds % 60;
    int minutes = seconds / 60;
    
    std::cout << "N: Refinement cycle: " << minutes << " minutes, "
    << finalSeconds << " seconds." << std::endl;
}

void MtzRefiner::initialMerge()
{
    MtzManager *originalMerge = NULL;
	std::vector<MtzPtr> mtzManagers = getAllMtzs();

    AmbiguityBreaker breaker = AmbiguityBreaker(mtzManagers);
    breaker.run();
    referencePtr = breaker.getMergedMtz();
    originalMerge = &*referencePtr;
    reference = referencePtr;
    
    reference->writeToFile("initialMerge.mtz");
}

void MtzRefiner::redumpBins()
{
	for (BinList::iterator it = binList.begin(); it != binList.end(); it++)
	{
		it->second.clear();
	}

	std::vector<MtzPtr> mtzs = getAllMtzs();

	for (int i = 0; i < mtzs.size(); i++)
	{
		int bin = mtzs[i]->getBin();

		if (!binList.count(bin))
		{
			binList[bin] = std::vector<MtzPtr>();
		}

		binList[bin].push_back(mtzs[i]);
	}
}

std::vector<MtzPtr> MtzRefiner::getAllMtzs()
{
	std::vector<MtzPtr> allMtzs;

	for (int i = 0; i < images.size(); i++)
	{
		ImagePtr image = images[i];

		for (int j = 0; j < image->mtzCount(); j++)
		{
			MtzPtr mtz = image->mtz(j);
			allMtzs.push_back(mtz);
		}
	}

	return allMtzs;
}

void MtzRefiner::refine()
{
    MtzPtr originalMerge;
    
    bool initialExists = loadInitialMtz();
	loadImageFiles();

    if (!initialExists)
    {
        initialMerge();
    }

    reference->writeToFile("originalMerge.mtz");
    
    MtzManager::setReference(&*reference);
    
    refineCycle();

    hasRefined = true;
}

void MtzRefiner::refineCycle(bool once)
{
	std::vector<MtzPtr> mtzManagers = getAllMtzs();

    int i = 0;
    bool finished = false;
    
    int maximumCycles = FileParser::getKey("MAXIMUM_CYCLES", 6);
	bool stop = FileParser::getKey("STOP_REFINEMENT", true);
    bool outputIndividualCycles = FileParser::getKey("OUTPUT_INDIVIDUAL_CYCLES", false);

	while (!finished)
    {
        if (outputIndividualCycles)
        {
            FileParser::setKey("CYCLE_NUMBER", i);
        }
        
        cycleNum = i;
        cycle();

		merge(i);

		if (once)
			finished = true;

		if (i == maximumCycles - 1 && stop)
			finished = true;

		i++;
	}
}


// MARK: Loading data

bool MtzRefiner::loadInitialMtz(bool force)
{
    bool hasInitialMtz = FileParser::hasKey("INITIAL_MTZ");

    if (reference && !force)
        return true;

    logged << "Initial MTZ has "
    << (hasInitialMtz ? "" : "not ") << "been provided." << std::endl;
    sendLog();
    
    if (hasInitialMtz)
    {
        std::string referenceFile = FileParser::getKey("INITIAL_MTZ",
                                                       std::string(""));
        
        reference = MtzPtr(new MtzManager());
        
        reference->setFilename(referenceFile.c_str());
		reference->loadReflections();
		reference->setDefaultMatrix();
		reference->setSigmaToUnity();
        
        if (reference->reflectionCount() == 0)
        {
            logged << "Initial MTZ reference missing or reflection count is 0. Exiting." << std::endl;
            sendLogAndExit();
        }
        
        MtzManager::setReference(&*reference);
    }

    return hasInitialMtz;
}


int MtzRefiner::imageSkip(size_t totalCount)
{
	int skip = FileParser::getKey("IMAGE_SKIP", 0);

	if (skip > totalCount)
	{
		std::ostringstream logged;
		logged << "IMAGE_SKIP specified is beyond the image count of the available images." << std::endl;
		LoggableObject::staticLogAndExit(logged);
	}

	return skip;
}

int MtzRefiner::imageMax(size_t lineCount)
{
    int skip = imageSkip(lineCount);
    
    int end = (int)lineCount;
    
    if (imageLimit != 0)
	{
        end = imageLimit < lineCount ? imageLimit : (int)lineCount;
	}

    end += skip;
    
    if (end > lineCount)
	{
        end = (int)lineCount;
	}

    return end;
}

void MtzRefiner::readSingleImageV2(std::string *filename, vector<ImagePtr> *newImages, vector<MtzPtr> *newMtzs, int offset, bool v3, MtzRefiner *me)
{
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    
    vector<double> givenUnitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    
    std::vector<std::string> hdf5Sources = FileParser::getKey("HDF5_SOURCE_FILES", std::vector<std::string>());
    bool readFromHdf5 = hdf5Sources.size() > 0;
    vector<double> cellDims = FileParser::getKey("UNIT_CELL", vector<double>());
    
    bool setSigmaToUnity = FileParser::getKey("SET_SIGMA_TO_UNITY", true);
    
    bool ignoreMissing = FileParser::getKey("IGNORE_MISSING_IMAGES", false);
    bool lowMemoryMode = FileParser::getKey("LOW_MEMORY_MODE", false);
    
    if (readFromHdf5)
    {
        ignoreMissing = true;
    }
    
    const std::string contents = FileReader::get_file_contents( filename->c_str());
    vector<std::string> imageList = FileReader::split(contents, "\nimage ");
    
    int maxThreads = FileParser::getMaxThreads();
    
    int skip = imageSkip(imageList.size());
    int end = imageMax(imageList.size());
    
    if (skip > 0)
    {
        std::ostringstream logged;
        logged << "Skipping " << skip << " lines" << std::endl;
        Logger::mainLogger->addStream(&logged);
    }
    
    for (int i = offset + skip; i < end; i += maxThreads)
    {
        if (imageList[i].length() == 0)
            continue;
        
        vector<std::string> lines = FileReader::split(imageList[i], '\n');
        
        vector<std::string> components = FileReader::split(lines[0], ' ');
        
        if (components.size() <= 1)
            continue;
        
        std::string imgName = components[1];
        std::string imgNameOnly = components[1];
        
        imgName.erase(std::remove(imgName.begin(), imgName.end(), '\r'), imgName.end());
        imgName.erase(std::remove(imgName.begin(), imgName.end(), '\n'), imgName.end());
        
        if (newImages)
            imgName += ".img";
        else if (newMtzs)
            imgName += ".mtz";
        
        std::ostringstream logged;
        
        if (!FileReader::exists(imgName) && !ignoreMissing && !v3)
        {
            logged << "Skipping image " << imgName << std::endl;
            Logger::mainLogger->addStream(&logged);
            continue;
        }
        
        if ((readFromHdf5 && newImages != NULL) || v3)
        {
            Hdf5ManagerCheetahPtr manager = Hdf5ManagerCheetah::hdf5ManagerForImage(imgNameOnly);
            
            if (!manager)
            {
                continue;
            }
        }
        
        double usedWavelength = wavelength;
        double usedDistance = detectorDistance;
        
        bool fromDials = FileParser::getKey("FROM_DIALS", false);
        
        if (components.size() >= 3)
        {
            usedWavelength = atof(components[1].c_str());
            usedDistance = atof(components[2].c_str());
        }
        
        MatrixPtr unitCell;
        MatrixPtr newMatrix;
        double delay = 0;
		double rlpSize = -1;
		double mosaicity = -1;

        ImagePtr newImage;
        
        if (readFromHdf5)
        {
            Hdf5ImagePtr hdf5Image = Hdf5ImagePtr(new Hdf5Image(imgName, wavelength,
                                                                0));
            newImage = boost::static_pointer_cast<Image>(hdf5Image);
        }
        else
        {
            newImage = ImagePtr(new Image(imgName, wavelength,
                                          0));
        }
        
        bool hasSpots = false;
        std::string parentImage = "";
		int currentCrystal = -1;
		int bin = 0;

		for (int i = 1; i < lines.size(); i++)
		{
			vector<std::string> components = FileReader::split(lines[i], ' ');

			if (components.size() == 0)
				continue;

			if (components[0] == "spots")
			{
				std::string spotsFile = components[1];
				newImage->setSpotsFile(spotsFile);
				hasSpots = true;
			}

			if (components[0] == "bin")
			{
				std::string binString = components[1];
				bin = atoi(binString.c_str());
			}

			if (components[0] == "distance_offset")
			{
				float offset = atof(components[1].c_str());
				Image::setDistanceOffset(&*newImage, offset);
			}

			if (components[0] == "crystal")
			{
				currentCrystal = atoi(components[1].c_str());
			}

			if (components[0] == "rlp_size")
			{
				rlpSize = atof(components[1].c_str());
			}

			if (components[0] == "mosaicity")
			{
				mosaicity = atof(components[1].c_str());
			}

			if (components[0] == "matrix")
			{
				// individual matrices
				double matrix[9];
				readMatrix(matrix, lines[i]);

				newMatrix = MatrixPtr(new Matrix(matrix));

				if (newImages)
				{
					newImage->setUpCrystal(newMatrix);
				}
			}

			if (components[0] == "wavelength")
			{
				double newWavelength = wavelength;

				if (components.size() >= 2)
					newWavelength = atof(components[1].c_str());

				newImage->setWavelength(newWavelength);

				logged << "Setting wavelength for " << imgName << " to " << newWavelength << " Angstroms." << std::endl;
				Logger::log(logged);
			}

			if (components[0] == "unitcell")
			{
				// individual matrices
				double matrix[9];
				readMatrix(matrix, lines[i]);

				unitCell = MatrixPtr(new Matrix(matrix));
			}

			if (components[0] == "rotation")
			{
				// individual matrices
				double matrix[9];
				readMatrix(matrix, lines[i]);

				MatrixPtr rotation = MatrixPtr(new Matrix(matrix));

				if (unitCell)
				{
					if (newImages)
					{
						vector<double> correction = FileParser::getKey("ORIENTATION_CORRECTION", vector<double>());

						if (fromDials)
						{
							rotation->rotate(0, 0, M_PI / 2);
							rotation->components[1] *= -1;
							rotation->components[5] *= -1;
							rotation->components[9] *= -1;
						}

						if (correction.size() >= 2)
						{
							double rightRot = correction[0] * M_PI / 180;
							double upRot = correction[1] * M_PI / 180;
							double swivelRot = 0;

							if (correction.size() > 2)
								swivelRot = correction[2] * M_PI / 180;

							rotation->rotate(rightRot, upRot, swivelRot);
						}
					}

					newMatrix = MatrixPtr(new Matrix);
					newMatrix->setComplexMatrix(unitCell, rotation);

					if (!v3)
					{
						currentCrystal++;
					}

				//	if (v3)
					{
						MtzPtr newManager = MtzPtr(new MtzManager());
						std::string prefix = (me->readRefinedMtzs ? "ref-" : "");
						newManager->setFilename((prefix + "img-" + imgNameOnly + "_" + i_to_str(currentCrystal) + ".mtz").c_str());
						newManager->setMatrix(newMatrix);

						if (setSigmaToUnity)
							newManager->setSigmaToUnity();

						if (rlpSize > 0)
						{
							newManager->setSpotSize(rlpSize);
						}
						if (mosaicity > 0)
						{
							newManager->setSpotSize(rlpSize);
						}

						newManager->setTimeDelay(delay);
						newManager->setImage(newImage);
						newManager->calcXYOffset();

						newManager->loadReflections();
						newManager->setWavelength(newImage->getWavelength());

						if (newManager->reflectionCount() > 0 && !v3)
						{
							newImage->addMtz(newManager);

							if (newMtzs)
							{
								newMtzs->push_back(newManager);
							}
						}
						else if (v3)
						{
							newImage->addMtz(newManager);

							if (!me->binList.count(bin))
							{
								me->binList[bin] = std::vector<MtzPtr>();
							}

							me->binList[bin].push_back(newManager);
							newManager->setBin(bin);
						}
					}
				}
			}

			if (components[0] == "delay" && newMtzs)
			{
				if (components.size() > 1)
					delay = atof(components[1].c_str());
			}
		}

        if (newImages)
        {
            newImages->push_back(newImage);
        }

        if (newMtzs && !v3)
        {
            MtzPtr newManager = MtzPtr(new MtzManager());

            newManager->setFilename(imgName.c_str());

			newManager->setImage(newImage);

            if (!newMatrix)
            {
                logged << "Warning! Matrix for " << imgName << " is missing." << std::endl;
                Logger::setShouldExit();
                Logger::log(logged);
            }

            newManager->setMatrix(newMatrix);

            if (!lowMemoryMode)
            {
                newManager->loadReflections();
            }

            if (setSigmaToUnity)
                newManager->setSigmaToUnity();

            newManager->setTimeDelay(delay);

            if (newManager->reflectionCount() > 0 || lowMemoryMode)
            {
                newMtzs->push_back(newManager);
            }
        }
    }
}

void MtzRefiner::readDataFromOrientationMatrixList(std::string *filename, bool areImages, std::vector<ImagePtr> *targetImages)
{
    double version = FileParser::getKey("MATRIX_LIST_VERSION", 2.0);

    // thought: turn the vector concatenation into a templated function

    boost::thread_group threads;

    int maxThreads = FileParser::getMaxThreads();

    vector<vector<ImagePtr> > imageSubsets;
    vector<vector<MtzPtr> > mtzSubsets;
    imageSubsets.resize(maxThreads);
    mtzSubsets.resize(maxThreads);

    std::string contents;

    try
    {
        contents = FileReader::get_file_contents( filename->c_str());
    }
    catch (int errno)
    {
        logged << "Missing file " << filename << ", cannot continue." << std::endl;
        sendLogAndExit();
    }
    vector<std::string> imageList = FileReader::split(contents, "\nimage ");

    std::ostringstream logged;
    
    if (imageList[0].substr(0, 6) == "ersion")
    {
        std::string vString = imageList[0].substr(7, imageList[0].length() - 7);
        float inputVersion = atof(vString.c_str());
        
        logged << "Autodetecting matrix list version: " << inputVersion << std::endl;
        
        version = inputVersion;
    }
    
    Logger::log(logged);
    
    if (version > 2.99 && version < 3.99)
    {
        loadPanels();
    }

    for (int i = 0; i < maxThreads; i++)
    {
        if (version > 1.99 && version < 2.99)
        {
            vector<MtzPtr> *chosenMtzs = areImages ? NULL : &mtzSubsets[i];
            vector<ImagePtr> *chosenImages = areImages ? &imageSubsets[i] : NULL;
            boost::thread *thr = new boost::thread(readSingleImageV2, filename,
                                                   chosenImages, chosenMtzs, i, false, this);
            threads.add_thread(thr);
        }
        else if (version > 2.99 && version < 3.99)
        {
            boost::thread *thr = new boost::thread(readSingleImageV2, filename,
                                                   &imageSubsets[i], &mtzSubsets[i], i, true, this);
            threads.add_thread(thr);
        }
    }
    
    threads.join_all();
    
    int total = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
		total += imageSubsets[i].size();
    }
    
    if (targetImages == NULL)
    {
        targetImages = &images;
    }
	targetImages->reserve(total);
	int lastPos = 0;

	for (int i = 0; i < maxThreads; i++)
	{
		targetImages->insert(targetImages->begin() + lastPos,
							 imageSubsets[i].begin(), imageSubsets[i].end());
		lastPos += imageSubsets[i].size();
	}

}

void MtzRefiner::readMatricesAndImages(std::string *filename, bool areImages, std::vector<ImagePtr> *targetImages)
{
    if (targetImages == NULL && images.size() > 0)
        return;
    
    std::string hdf5OutputFile = FileParser::getKey("HDF5_OUTPUT_FILE", std::string(""));
    std::vector<std::string> hdf5ProcessingFiles = FileParser::getKey("HDF5_SOURCE_FILES", std::vector<std::string>());
    
    bool hdf5 = hdf5OutputFile.length() || hdf5ProcessingFiles.size();
    Hdf5ManagerCheetah::initialiseCheetahManagers();

    std::string aFilename = "";
    aFilename = FileParser::getKey("ORIENTATION_MATRIX_LIST", std::string(""));
    
    if (!hdf5)
    {
        if (aFilename.length() == 0)
        {
            logged << "No orientation matrix list provided (and no HDF5 file either). Exiting now." << std::endl;
            sendLogAndExit();
        }
    }
    
    if (aFilename.length() && !FileReader::exists(aFilename))
    {
        std::ostringstream logged;
        logged << "File specified in ORIENTATION_MATRIX_LIST " << aFilename << " doesn't exist." << std::endl;
        
        if (hdf5)
        {
            logged << "You have specified HDF5 image files, so you don't necessarily need to specify your ORIENTATION_MATRIX_LIST." << std::endl;
        }
        else
        {
            logged << "You must either specify some HDF5 image files using HDF5_SOURCE_FILES, or you must specify ORIENTATION_MATRIX_LIST to load .img files from the file system instead." << std::endl;
        }
        
        LoggableObject::staticLogAndExit(logged);
    }
    
    if (aFilename.length())
    {
        readDataFromOrientationMatrixList(&aFilename, areImages, targetImages);
    }
    
    if (hdf5)
    {
        readFromHdf5(&images);
    }
    
    if (hdf5)
    {
        Hdf5ImagePtr maskImage = Hdf5ImagePtr(new Hdf5Image("tag-mask"));
        maskImage->setMask();
    }

	std::vector<MtzPtr> mtzManagers = getAllMtzs();

	logged << "Summary\nImages: " << images.size() << "\nCrystals: " << mtzManagers.size() << std::endl;
    sendLog();


	bool shouldApplyUnrefinedPartiality = FileParser::getKey("APPLY_UNREFINED_PARTIALITY", false);

	if (shouldApplyUnrefinedPartiality)
	{
		logged << "Applying unrefined partiality to images." << std::endl;
		sendLog();

		for (int i = 0; i < mtzManagers.size(); i++)
		{
			mtzManagers[i]->applyUnrefinedPartiality();
		}
	}
}

void MtzRefiner::readFromHdf5(std::vector<ImagePtr> *newImages)
{
    // in the case where HDF5_OUTPUT_FILE is set.
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    std::vector<double> beamCentre = FileParser::getKey("BEAM_CENTRE", std::vector<double>(0, 2));
    vector<double> givenUnitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    
    int count = 0;
    
    std::string orientationList = FileParser::getKey("ORIENTATION_MATRIX_LIST", std::string(""));
    bool hasList = (orientationList.length() > 0);
    
    std::string hdf5OutputFile = FileParser::getKey("HDF5_OUTPUT_FILE", std::string(""));

    if (hdf5OutputFile.length())
    {
        hdf5ProcessingPtr = Hdf5ManagerProcessingPtr(new Hdf5ManagerProcessing(hdf5OutputFile));
    }
    
    int startImage = FileParser::getKey("IMAGE_SKIP", 0);
    int endImage = startImage + FileParser::getKey("IMAGE_LIMIT", 0);
    
    int loadAll = (startImage == 0 && endImage == 0);
    
    if (!loadAll)
    {
        logged << "Starting from image " << startImage << " and ending on image " << endImage << std::endl;
    }
    else
    {
        logged << "Loading all images from HDF5 source files." << std::endl;
    }

    sendLog();
    
    int inputHdf5Count = Hdf5ManagerCheetah::cheetahManagerCount();
    
    for (int i = 0; i < inputHdf5Count; i++)
    {
        Hdf5ManagerCheetahPtr manager = Hdf5ManagerCheetah::cheetahManager(i);
        
        for (int j = 0; j < manager->imageAddressCount(); j++)
        {
            if (!loadAll)
            {
                if (count < startImage || count >= endImage)
                {
                    count++;
                    continue;
                }
            }
            
            std::string imageAddress = manager->imageAddress(j);
            
            std::string imgName = manager->lastComponent(imageAddress) + ".img";
            
            if (!hasList)
            {
                Hdf5ImagePtr hdf5Image = Hdf5ImagePtr(new Hdf5Image(imgName, wavelength,
                                                                    detectorDistance));
                ImagePtr newImage = boost::static_pointer_cast<Image>(hdf5Image);
                hdf5Image->loadCrystals();
                newImages->push_back(newImage);
            }

            count++;
        }
    }

    logged << "N: Images loaded from HDF5: " << newImages->size() << std::endl;;
    
    sendLog();
}

void MtzRefiner::readMatricesAndMtzs()
{
    readMatricesAndImages(NULL, false);
}

void MtzRefiner::readMatrix(double (&matrix)[9], std::string line)
{
    vector<std::string> components = FileReader::split(line, ' ');
    
    for (int j = 1; j <= 9; j++)
    {
        std::string component = components[j];
        
        /* Locate the substring to replace. */
        int index = (int)component.find(std::string("*^"), 0);
        if (index != std::string::npos)
        {
            component.replace(index, 2, std::string("e"));
        }
        
        double matVar = atof(component.c_str());
        matrix[j - 1] = matVar;
    }
}

// MARK: Merging

void MtzRefiner::correlationAndInverse(bool shouldFlip)
{
    if (MtzManager::getReferenceManager() == NULL)
    {
        MtzManager::setReference(&*reference);
    }

	std::vector<MtzPtr> mtzManagers = getAllMtzs();

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        double correl = mtzManagers[i]->correlation(true);
        double invCorrel = correl;
        
        if (MtzManager::getReferenceManager()->ambiguityCount() > 1)
        {
            mtzManagers[i]->setActiveAmbiguity(1);
            invCorrel = mtzManagers[i]->correlation(true);
            mtzManagers[i]->setActiveAmbiguity(0);
            
            if (invCorrel > correl && shouldFlip)
                mtzManagers[i]->setActiveAmbiguity(1);
        }
        double newCorrel = mtzManagers[i]->correlation(true);
        mtzManagers[i]->setRefCorrelation(newCorrel);
    }
}

void MtzRefiner::linearScaling()
{
    loadInitialMtz();
    readMatricesAndMtzs();
	std::vector<MtzPtr> mtzManagers = getAllMtzs();

    MtzMerger merger;
    merger.setSomeMtzs(mtzManagers);
    merger.setCycle(-1);
    merger.setScalingType(ScalingTypeReference);
    merger.setFilename("scaling.mtz");
    merger.scale();
    
    referencePtr = merger.getMergedMtz();
}

void MtzRefiner::merge(int cycle)
{
	loadImageFiles();
	loadInitialMtz();

	std::vector<MtzPtr> mtzManagers = getAllMtzs();

	if (cycle < -1)
	{
		correlationAndInverse();
	}

    bool anomalousMerge = FileParser::getKey("MERGE_ANOMALOUS", false);
	int scalingInt = FileParser::getKey("SCALING_STRATEGY",
										(int) SCALING_STRATEGY);
	ScalingType scaling = (ScalingType) scalingInt;
	bool replaceReference = FileParser::getKey("REPLACE_REFERENCE", true);

	MtzPtr mergedMtz;
	std::string filename = "allMerge" + i_to_str(cycle) + ".mtz";

	MtzManager::setReference(&*reference);

	MtzMerger merger;
	merger.setAllMtzs(mtzManagers);
	merger.setCycle(cycleNum);
	merger.setScalingType(scaling);
	merger.mergeFull();
	mergedMtz = merger.getMergedMtz();

	referencePtr = mergedMtz;

	if (anomalousMerge)
	{
		merger.mergeFull(true);
	}

	// *************************
	// ***** BINNED MERGE ******
	// *************************

	redumpBins();

	for (int i = 0; i < binList.size(); i++)
	{
		if (binList.size() <= 1)
		{
			break;
		}

		MtzMerger merger;
		merger.setAllMtzs(binList[i]);
		merger.setCycle(cycleNum);
		merger.setScalingType(scaling);
		merger.setFilename("allBin" + i_to_str(i) +
						   "_" + i_to_str(cycleNum) + ".mtz");
		merger.mergeFull();
		mergedMtz = merger.getMergedMtz();
	}

	if (replaceReference && mergedMtz)
	{
		reference = mergedMtz;
		MtzManager::setReference(&*reference);
	}
}

// MARK: Integrating images

void MtzRefiner::integrateImagesWrapper(MtzRefiner *object,
                                        vector<MtzPtr> *&mtzSubset, int offset)
{
    object->integrateImages(mtzSubset, offset);
}

void MtzRefiner::integrateImages(vector<MtzPtr> *&mtzSubset,
                                 int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < images.size(); i += maxThreads)
    {
        std::ostringstream logged;
        logged << "Integrating image " << i << std::endl;
        Logger::mainLogger->addStream(&logged);

        images[i]->refineOrientations();
        
        vector<MtzPtr> mtzs = images[i]->currentMtzs();

        mtzSubset->insert(mtzSubset->end(), mtzs.begin(), mtzs.end());
        images[i]->dropImage();
    }
}

void MtzRefiner::maximumImageThread(MtzRefiner *me, ImagePtr maxImage, int offset)
{
    std::vector<ImagePtr> selection;
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < me->images.size(); i += maxThreads)
    {
        selection.push_back(me->images[i]);
    }
    
    maxImage->makeMaximumFromImages(selection);
}

void MtzRefiner::maximumImage()
{
    if (images.size() == 0)
    {
        loadImageFiles();
    }
    
    ImagePtr maximum = ImagePtr(new Image("maxImage.img", 0, 0));
    int maxThreads = FileParser::getMaxThreads();
    
    std::vector<ImagePtr> maxSelections;
    
    for (int i = 0; i < maxThreads; i++)
    {
        ImagePtr anImage = ImagePtr(new Image("maxImage_" + i_to_str(i) + ".img", 0, 0));
        maxSelections.push_back(anImage);
    }

    boost::thread_group threads;
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(maximumImageThread, this, maxSelections[i], i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    maximum->makeMaximumFromImages(maxSelections, true);
    Detector::setDrawImage(maximum);
}

void MtzRefiner::loadImageFiles()
{
    loadPanels();
    
    std::string filename = FileParser::getKey("ORIENTATION_MATRIX_LIST",
                                              std::string(""));
    std::string hdf5 = FileParser::getKey("HDF5_OUTPUT_FILE",
                                          std::string(""));
    
    readMatricesAndImages(&filename);
    
    if (images.size() > 0)
    {
        Detector::setDrawImage(images[0]);
        // force loading of everything
        images[0]->valueAt(0, 0);
        SpotPtr spot = SpotPtr(new Spot(images[0]));
    }
}

void MtzRefiner::integrate()
{
    loadImageFiles();
    
    int crystals = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        crystals += images[i]->mtzCount();
    }
    
    logged << images.size() << " images with " << crystals << " crystal orientations." << std::endl;
    sendLog();

    int maxThreads = FileParser::getMaxThreads();
    
    boost::thread_group threads;
    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(integrateImagesWrapper, this,
                                               &managerSubsets[i], i);
        threads.add_thread(thr);
    }
    
    threads.join_all();

    writeNewOrientations(false, true);
    
    integrationSummary();
}


void MtzRefiner::integrationSummary()
{
    std::ostringstream refineSummary;
    
    int imageCount = (int)images.size();
    int mtzCount = 0;
	int goodSpots = 0;

    refineSummary << MtzManager::parameterHeaders() << std::endl;
    
    for (int i = 0; i < images.size(); i++)
    {
		if (images[i]->acceptableSpotCount())
		{
			goodSpots++;
		}

        for (int j = 0; j < images[i]->mtzCount(); j++)
        {
            refineSummary << images[i]->mtz(j)->writeParameterSummary() << std::endl;

            mtzCount++;
        }
    }
    
    std::string summaryString = refineSummary.str();
    Logger::mainLogger->addStream(&refineSummary);
    
    std::replace(summaryString.begin(), summaryString.end(), '\t', ',');
    
    std::ofstream summaryCSV;
    summaryCSV.open("integration.csv");
    summaryCSV << summaryString;
    summaryCSV.close();
    
    float pctIndexed = (float)mtzCount / (float)imageCount * 100;
	float pctPlausible = (float)goodSpots / (float)imageCount * 100;
	logged << "N: Spot count plausible for " << goodSpots << " images out of " << imageCount << " images (" << pctPlausible << "%)." << std::endl;
	logged << "N: Indexed " << mtzCount << " crystals out of " << imageCount << " images (" << pctIndexed << "%)." << std::endl;
    sendLog();
    
    Logger::mainLogger->addString("Written integration summary to integration.csv");
    
}

void MtzRefiner::loadPanels(bool mustFail)
{
    bool useNewDetectorFormat = Detector::isActive();
    
    if (!useNewDetectorFormat && mustFail)
    {
        logged << "DETECTOR_LIST has not been provided." << std::endl;
        FileParser::printCommandInfo("DETECTOR_LIST");
        sendLogAndExit();
        return;
    }
    
    if (useNewDetectorFormat)
    {
        if (Detector::getMaster())
        {
            return;
        }
        
        std::string detectorList = FileParser::getKey("DETECTOR_LIST", std::string(""));
        
        logged << "Loading new detector format from " << detectorList << "." << std::endl;
        sendLog();
        int formatInt = FileParser::getKey("GEOMETRY_FORMAT", 0);
        
        GeometryFormat format = (GeometryFormat)formatInt;
        
        GeometryParser geomParser = GeometryParser(detectorList, format);
        geomParser.parse();
        
        return;
    }
}

// MARK: indexing

void MtzRefiner::writeAllNewOrientations()
{
    std::string filename = FileParser::getKey("NEW_MATRIX_LIST", std::string("new_orientations.dat"));
    bool includeRots = false;
    
    std::ofstream allMats;
    allMats.open(FileReader::addOutputDirectory("all-" + filename));
    
    allMats << "version 3.0" << std::endl;
    
    for (int i = 0; i < images.size(); i++)
    {
        std::string imageName = images[i]->getBasename();
        allMats << "image " << imageName << std::endl;
        
		allMats << "distance_offset " << Image::getDistanceOffset(&*images[i]) << std::endl;

		if (images[i]->getSpotsFile().length() > 0)
        {
            allMats << "spots " << images[i]->getSpotsFile() << std::endl;
        }
        
        for (int j = 0; j < images[i]->mtzCount(); j++)
        {
            MtzPtr mtz = images[i]->mtz(j);
            int crystalNum = j;
            
			allMats << "crystal " << crystalNum << std::endl;
			allMats << mtz->getOrientationMatrixListScript();
        }
    }
    
    Logger::mainLogger->addString("Written to matrix list: all-" + filename);
    allMats.close();
}

void MtzRefiner::writeNewOrientations(bool includeRots, bool detailed)
{
    std::ofstream mergeMats;
    std::ofstream refineMats;
    std::ofstream integrateMats;

    std::string filename = FileParser::getKey("NEW_MATRIX_LIST", std::string("new_orientations.dat"));
    
    mergeMats.open(FileReader::addOutputDirectory("merge-" + filename));
    refineMats.open(FileReader::addOutputDirectory("refine-" + filename));
    integrateMats.open(FileReader::addOutputDirectory("integrate-" + filename));
	std::vector<MtzPtr> mtzManagers = getAllMtzs();

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        MtzPtr manager = mtzManagers[i];
        
        // write out matrices etc.
        std::string imgFilename = manager->getBasename();
        
        ImagePtr image = manager->getImagePtr();
        std::string parentName;
        
        if (image)
        {
            parentName = image->getBasename();
        }
        
        MatrixPtr matrix = manager->getMatrix()->copy();
        
        if (includeRots)
        {
			double hRad = MtzManager::getHRot(&*manager) * M_PI / 180;
			double kRad = MtzManager::getKRot(&*manager) * M_PI / 180;

            matrix->rotate(hRad, kRad, 0);
        }
        
        std::string prefix = "img-";
        
        if (imgFilename.substr(0, 3) == "img")
            prefix = "";
        
        
        if (!detailed)
        {
            refineMats << prefix << imgFilename << " ";
            std::string description = matrix->description();
            refineMats << description << std::endl;
        }
        else
        {
            refineMats << "image " << prefix << imgFilename << std::endl;
            mergeMats << "image ref-" << prefix << imgFilename << std::endl;
            mergeMats << "parent " << imgFilename << std::endl;
            std::string description = matrix->description(true);
            refineMats << description << std::endl;
            mergeMats << description << std::endl;
            refineMats << "wavelength " << image->getWavelength() << std::endl;
            
            if (hasRefined)
            {
                refineMats << manager->getParamLine() << std::endl;
            }
        }
    }
    
    for (int i = 0; i < images.size(); i++)
    {
        ImagePtr image = images[i];
        
        // write out matrices etc.
        std::string imgFilename = image->getBasename();
        
        integrateMats << "image " << imgFilename << std::endl;
        
        for (int j = 0; j < image->mtzCount(); j++)
        {
            MatrixPtr matrix = image->mtz(j)->getMatrix()->copy();
            
            std::string desc90 = matrix->description(true);
            integrateMats << desc90 << std::endl;
        }
        
        if (image->getSpotsFile().length() > 0)
        {
            integrateMats << "spots " << image->getSpotsFile() << std::endl;
        }
        
        integrateMats << "wavelength " << image->getWavelength() << std::endl;
    }
    
    if (!images.size())
    {
		std::vector<MtzPtr> mtzManagers = getAllMtzs();

		for (int i = 0; i < mtzManagers.size(); i++)
        {
            std::string filename = mtzManagers[i]->getFilename();
            filename.erase(filename.end() - 6, filename.end());
            filename.erase(filename.begin(), filename.begin() + 4);
            
            integrateMats << "image " << filename << std::endl;
            
            MatrixPtr matrix = mtzManagers[i]->getMatrix()->copy();
            
            std::string desc90 = matrix->description(true);
            integrateMats << desc90 << std::endl;

            integrateMats << "wavelength " << mtzManagers[i]->getWavelength() << std::endl;
        }
    }
    
    Logger::mainLogger->addString("Written to filename set: " + filename);
    
    refineMats.close();
    mergeMats.close();
    integrateMats.close();
    
    writeAllNewOrientations();
}

void MtzRefiner::findSpotsThread(MtzRefiner *me, int offset)
{
    double maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < me->images.size(); i += maxThreads)
    {
        me->images[i]->processSpotList();
        me->images[i]->dropImage();
    }
}

void MtzRefiner::index()
{
    loadImageFiles();
    
    logged << "N: Total images loaded: " << images.size() << std::endl;
    sendLog();
    
    if (!indexManager)
        indexManager = new IndexManager(images);
    
    indexManager->index();

	takeTwoPNG();
    
    writeNewOrientations(false, true);
    integrationSummary();
}

void MtzRefiner::powderPattern()
{
    loadImageFiles();
    
    if (!indexManager)
        indexManager = new IndexManager(images);
    
    double maxThreads = FileParser::getMaxThreads();
    boost::thread_group threads;
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(findSpotsThread, this, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    indexManager->powderPattern();
}

// MARK: SACLA stuff for April 2016


void MtzRefiner::integrateSpotsThread(MtzRefiner *me, int offset)
{
    double maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < me->images.size(); i += maxThreads)
    {
        me->images[i]->integrateSpots();
    }
}

void MtzRefiner::integrateSpots()
{
    loadPanels();
    this->readMatricesAndImages();
    std::cout << "N: Total images loaded: " << images.size() << std::endl;
    double maxThreads = FileParser::getMaxThreads();
    
    boost::thread_group threads;
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(integrateSpotsThread, this, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    std::cout << "Finished integrating spots." << std::endl;
}

// MARK: Miscellaneous

void MtzRefiner::plotPixelValueVsFiducial()
{
	loadImageFiles();

	CSVPtr csv = CSVPtr(new CSV(2, "fiducial", "value"));

	std::vector<int> pix = FileParser::getKey("SPECIAL_PIXEL", std::vector<int>(50, 2));

	int x = pix[0];
	int y = pix[1];

	for (int i = 0; i < images.size(); i++)
	{
		ImagePtr image = images[i];

		double value = image->valueAt(x, y);
		double fiducial = image->getFrameNumber();

		if (fiducial < 0)
		{
			logged << "Fiducial not found for image " << image->getFilename() << std::endl;
			sendLog();
			continue;
		}

		csv->addEntry(2, fiducial, value);

		image->dropImage();
	}

	std::string filename = "pixel_vs_fiducial_" + i_to_str(x) +  "_" +  i_to_str(y);
	csv->writeToFile(filename + ".csv");

	std::map<std::string, std::string> plotMap;
	plotMap["filename"] = filename;
	plotMap["xHeader0"] = "fiducial";
	plotMap["yHeader0"] = "value";
	plotMap["xTitle0"] = "Fiducial for frame";
	plotMap["yTitle0"] = "Pixel value ADU";
	plotMap["style0"] = "line";

	csv->plotPNG(plotMap);
}

void MtzRefiner::plotIntegrationWindows()
{
    int h, k, l;
    std::vector<int> hkl = FileParser::getKey("MILLER_INDEX", std::vector<int>());
    GraphDrawer drawer = GraphDrawer(NULL);
    
    if (hkl.size() < 3)
    {
        logged << "Please set MILLER_INDEX" << std::endl;
        sendLog();
    }
    else
    {
        h = hkl[0];
        k = hkl[1];
        l = hkl[2];
        
		std::vector<MtzPtr> mtzManagers = getAllMtzs();
        logged << "Drawing integration windows for " << h << " " << k << " " << l << " for " << mtzManagers.size() << " mtzs." << std::endl;
        sendLog();

		drawer.cutoutIntegrationAreas(mtzManagers, h, k, l);
    }

}

void MtzRefiner::differenceCorrelation()
{
	std::vector<MtzPtr> mtzs = getAllMtzs();
	AmbiguityBreaker breaker = AmbiguityBreaker(mtzs);
	breaker.plotDifferences();

}

void MtzRefiner::plotIntensities()
{
    int h, k, l;
    std::vector<int> hkl = FileParser::getKey("MILLER_INDEX", std::vector<int>());
    GraphDrawer drawer = GraphDrawer(NULL);
    
    if (hkl.size() < 3)
    {
        logged << "To plot a specific Miller index please use MILLER_INDEX; this will dump a huge number" << std::endl;
        sendLog();
        drawer.plotReflectionFromMtzs(getAllMtzs());
    }
    else
    {
        h = hkl[0];
        k = hkl[1];
        l = hkl[2];
        
        drawer.plotReflectionFromMtzs(getAllMtzs(), h, k, l);
    }
    
    
}

void MtzRefiner::refineUnitCell()
{
	std::vector<MtzPtr> mtzManagers = getAllMtzs();
	UnitCellLattice::getMainLattice()->refineMtzs(mtzManagers);

	FileParser::setKey("FIX_UNIT_CELL", false);
}

// MARK: geometry

void MtzRefiner::flattenDetector()
{
	loadImageFiles();

	Detector::getMaster()->flattenDetector();
}

void MtzRefiner::refineMetrology(bool global)
{
    loadImageFiles();
	std::vector<MtzPtr> mtzManagers = getAllMtzs();
	Detector::getMaster()->clearMillers();

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->millersToDetector();
    }
    
    GeometryRefiner refiner;
    refiner.setImages(images);
    refiner.refineGeometry();
}

void MtzRefiner::reportMetrology()
{
    GeometryRefiner refiner;
    refiner.setImages(images);
    refiner.startingGraphs();
    refiner.reportProgress();
    
    Detector::getMaster()->reportMillerScores();
}

void MtzRefiner::fakeSpotsThread(std::vector<ImagePtr> *images, int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < images->size(); i += maxThreads)
    {
        ImagePtr image = images->at(i);
        
        image->fakeSpots();
    }
}

void MtzRefiner::fakeSpots()
{
    loadImageFiles();
    
    boost::thread_group threads;
    
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(fakeSpotsThread, &images, i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    writeAllNewOrientations();
}

MtzRefiner::~MtzRefiner()
{
    //    std::cout << "Deallocating MtzRefiner." << std::endl;
    
    images.clear();
    
    // TODO Auto-generated destructor stub
}

void MtzRefiner::displayIndexingHands()
{
	std::vector<MtzPtr> mtzManagers = getAllMtzs();
	std::ostringstream idxLogged, invLogged;
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        std::ostringstream &which = mtzManagers[i]->getActiveAmbiguity() == 0 ? idxLogged : invLogged;
        
        which << mtzManagers[i]->getFilename() << std::endl;
    }
    
    Logger::mainLogger->addString("**** Indexing hands ****\n - First indexing hand");
    Logger::mainLogger->addStream(&idxLogged);
    Logger::mainLogger->addString("- Second indexing hand");
    Logger::mainLogger->addStream(&invLogged);
    
}

void MtzRefiner::applyUnrefinedPartiality()
{
	std::vector<MtzPtr> mtzManagers = getAllMtzs();

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->applyUnrefinedPartiality();
    }
}

void MtzRefiner::orientationPlot()
{
    GraphDrawer drawer = GraphDrawer(NULL);
    drawer.plotOrientationStats(getAllMtzs());
}

void MtzRefiner::writePNGs(int total)
{
    loadPanels();
    
    if (!images.size())
    {
        readMatricesAndImages();
    }
    
    if (total == 0)
    {
        total = FileParser::getKey("PNG_TOTAL", 10);
    }
    
    int totalImages = (int)images.size();
    int skip = totalImages / total;
    
    if (skip == 0)
        skip = 1;
    
    bool allLattices = FileParser::getKey("PNG_ALL_LATTICES", false);
    
    for (int i = 0; i < totalImages; i += skip)
    {
        images[i]->drawSpotsOnPNG();
        
        if (allLattices && images[i]->mtzCount() > 0)
        {
            images[i]->drawCrystalsOnPNG(-1);
        }
        else
        {
            for (int j = 0; j < images[i]->mtzCount(); j++)
            {
                images[i]->drawCrystalsOnPNG(j);
            }
        }
    }
}

void MtzRefiner::takeTwoPNG()
{
    ImagePtr image = Detector::getSpecialImage();
    
    if (!image)
    {
        logged << "No images/special image loaded, cannot do TakeTwo vector plot." << std::endl;
        sendLog();
        return;
    }
    
    image->plotTakeTwoVectors(images);
}

void MtzRefiner::imageToDetectorMap()
{
    if (images.size())
    {
        double min, max;
        Detector::getMaster()->zLimits(&min, &max);
    }
}