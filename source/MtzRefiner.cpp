/*
 * MtzRefiner.cpp
 *
 *  Created on: 4 Jan 2015
 *      Author: helenginn
 */

#include "MtzRefiner.h"
#include <vector>
#include "lbfgs_cluster.h"
#include "misc.h"
#include "Vector.h"
#include <boost/thread/thread.hpp>
#include "GraphDrawer.h"
#include "Image.h"
#include <fstream>
#include "AmbiguityBreaker.h"

#include "FileParser.h"
#include "Index.h"
#include "parameters.h"
#include "FileReader.h"
#include "PanelParser.h"
#include "Panel.h"
#include "Logger.h"
#include "XManager.h"

bool MtzRefiner::hasPanelParser;
int MtzRefiner::imageLimit;
int MtzRefiner::cycleNum;

MtzRefiner::MtzRefiner()
{
    // TODO Auto-generated constructor stub
    reference = NULL;
    panelParser = NULL;
    
    hasPanelParser = false;
    int logInt = FileParser::getKey("VERBOSITY_LEVEL", 0);
    Logger::mainLogger->changePriorityLevel((LogLevel)logInt);
    
    imageLimit = FileParser::getKey("IMAGE_LIMIT", 0);
}

void print_parameters()
{
    std::cout << "N: ===== Default parameters =====" << std::endl;
    
    std::cout << "N: Initial parameters for MTZ grid search:" << std::endl;
    std::cout << "N: Bandwidth: " << INITIAL_BANDWIDTH << std::endl;
    std::cout << "N: Mosaicity: " << INITIAL_MOSAICITY << std::endl;
    std::cout << "N: Spot size: " << INITIAL_SPOT_SIZE << std::endl;
    std::cout << "N: Exponent: " << INITIAL_EXPONENT << std::endl;
    std::cout << std::endl;
    
    std::cout << "N: Parameters refined during grid search: ";
    std::cout << (OPTIMISED_ROT ? "rotations; " : "");
    std::cout << (OPTIMISED_BANDWIDTH ? "bandwidth; " : "");
    std::cout << (OPTIMISED_WAVELENGTH ? "wavelength; " : "");
    std::cout << (OPTIMISED_EXPONENT ? "exponent; " : "");
    std::cout << (OPTIMISED_MOSAICITY ? "mosaicity; " : "");
    std::cout << (OPTIMISED_SPOT_SIZE ? "spot size; " : "");
    std::cout << std::endl << std::endl;
    
    std::cout << "N: Grid search steps: ";
    std::cout << "N: Bandwidth: " << BANDWIDTH_STEP << std::endl;
    std::cout << "N: Rotation: " << ROT_STEP << std::endl;
    std::cout << "N: Spot size: " << SPOT_STEP << std::endl;
    std::cout << "N: Exponent: " << EXPONENT_STEP << std::endl;
    std::cout << "N: Wavelength: " << MEAN_STEP << std::endl;
    std::cout << std::endl;
    
    std::cout << "N: Grid search tolerances: ";
    std::cout << "N: Bandwidth: " << BANDWIDTH_TOLERANCE << std::endl;
    std::cout << "N: Rotation: " << ROT_TOLERANCE << std::endl;
    std::cout << "N: Spot size: " << SPOT_STEP << std::endl;
    std::cout << "N: Exponent: " << SPOT_SIZE_TOLERANCE << std::endl;
    std::cout << "N: Wavelength: " << MEAN_TOLERANCE << std::endl;
    std::cout << std::endl;
    
    std::cout << "N: Maximum resolution used for most grid searches: "
    << MAX_OPTIMISATION_RESOLUTION << std::endl;
    std::cout << "N: Maximum resolution used for spot size: "
    << MAX_SPOT_SIZE_OPT_RESOLUTION << std::endl;
    std::cout << "N: Further optimisation is "
    << (FURTHER_OPTIMISATION ? "" : "not ") << "used." << std::endl;
    
    std::cout << "N: Default detector distance: " << DEFAULT_DETECTOR_DISTANCE
    << std::endl;
    std::cout << "N: Default wavelength: " << DEFAULT_WAVELENGTH << std::endl;
    std::cout << "N: Intensity threshold: " << INTENSITY_THRESHOLD << std::endl;
    
    std::cout << std::endl << "N: Polarisation correction: "
    << (POLARISATION_CORRECTION ? "on" : "off") << std::endl;
    std::cout << "N: Polarisation factor (horizontal): "
    << HORIZONTAL_POLARISATION_FACTOR << std::endl;
    std::cout << std::endl << "N: Minimum miller count for rejection: "
    << MIN_MILLER_COUNT << std::endl;
    std::cout << "N: Rejecting miller indices: "
    << (REJECTING_MILLERS ? "on" : "off") << std::endl;
}

// MARK: Refinement

void MtzRefiner::cycleThreadWrapper(MtzRefiner *object, int offset)
{
    object->cycleThread(offset);
}

void MtzRefiner::cycleThread(int offset)
{
    int img_num = (int)mtzManagers.size();
    int j = 0;
    
    bool initialGridSearch = FileParser::getKey("INITIAL_GRID_SEARCH", false);
    
    int maxThreads = FileParser::getMaxThreads();
    
    for (int i = offset; i < img_num; i += maxThreads)
    {
        j++;
        
        ostringstream logged;
        
        MtzPtr image = mtzManagers[i];
        
        if (!image->isRejected())
        {
            logged << "Refining image " << i << " " << image->getFilename() << std::endl;
            Logger::mainLogger->addStream(&logged);
            if (initialGridSearch && cycleNum == 0)
                image->findSteps();
            image->gridSearch(true);
            
            image->writeToDat();
        }
    }
}

void MtzRefiner::cycle()
{
    MtzManager::setReference(reference);
    
    time_t startcputime;
    time(&startcputime);
    
    boost::thread_group threads;
    
    int maxThreads = FileParser::getMaxThreads();
    
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
    
  /*  Lbfgs_Cluster *lbfgs = new Lbfgs_Cluster();
    lbfgs->initialise_cluster_lbfgs(mtzManagers, &originalMerge);
    reference = originalMerge;
    delete lbfgs;
    
    reference->writeToFile("initialMerge.mtz");*/
    
     AmbiguityBreaker breaker = AmbiguityBreaker(mtzManagers);
     breaker.run();
     originalMerge = breaker.getMergedMtz();
     reference = originalMerge;
}

void MtzRefiner::refine()
{
    MtzManager *originalMerge = NULL;
    
    readMatricesAndMtzs();
    
    if (!loadInitialMtz())
    {
        initialMerge();
    }
    
    bool denormalise = FileParser::getKey("DENORMALISE_PARTIALITY", false);
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        if (denormalise)
            mtzManagers[i]->denormaliseMillers();
    }
    
    originalMerge = reference;
    
    std::cout << "N: Total images loaded: " << mtzManagers.size() << std::endl;
    
    originalMerge->writeToFile("originalMerge.mtz");
    
    MtzManager::currentManager = originalMerge;
    MtzManager::setReference(reference);
    double correl = originalMerge->correlation(true);
    std::cout << "Merged correlation = " << correl << std::endl;
    
    MtzManager::setLowRes(0);
    MtzManager::setHighRes(0);
    
    refineCycle();
}

void MtzRefiner::refineCycle(bool once)
{
    int i = 0;
    bool finished = false;
    
    int minimumCycles = FileParser::getKey("MINIMUM_CYCLES", 6);
    
    bool stop = FileParser::getKey("STOP_REFINEMENT", true);
    double correlationThreshold = FileParser::getKey("CORRELATION_THRESHOLD",
                                                     CORRELATION_THRESHOLD);
    double initialCorrelationThreshold = FileParser::getKey(
                                                            "INITIAL_CORRELATION_THRESHOLD", INITIAL_CORRELATION_THRESHOLD);
    int thresholdSwap = FileParser::getKey("THRESHOLD_SWAP",
                                           THRESHOLD_SWAP);
    bool exclusion = FileParser::getKey("EXCLUDE_OWN_REFLECTIONS", false);
    
    double resolution = FileParser::getKey("MAX_RESOLUTION_ALL",
                                           MAX_OPTIMISATION_RESOLUTION);
    int scalingInt = FileParser::getKey("SCALING_STRATEGY",
                                        (int) SCALING_STRATEGY);
    ScalingType scaling = (ScalingType) scalingInt;
    
    while (!finished)
    {
        cycleNum = i;
        cycle();
        
        std::cout << "Grouping final MTZs" << std::endl;
        MtzGrouper *grouper = new MtzGrouper();
        MtzManager::setReference(reference);
        if (i >= 0)
            grouper->setScalingType(scaling);
        if (once)
            grouper->setScalingType(ScalingTypeAverage);
        grouper->setWeighting(WeightTypePartialitySigma);
        grouper->setExpectedResolution(resolution);
        grouper->setMtzManagers(mtzManagers);
        grouper->setExcludeWorst(true);
        if (i < thresholdSwap)
            grouper->setCorrelationThreshold(initialCorrelationThreshold);
        else
            grouper->setCorrelationThreshold(correlationThreshold);
        
        MtzManager *mergedMtz = NULL;
        MtzManager *unmergedMtz = NULL;
        grouper->merge(&mergedMtz, &unmergedMtz, i);
        
        MtzManager::currentManager = mergedMtz;
        MtzManager::setReference(reference);
        
        std::cout << "Holders: " << mergedMtz->holderCount() << std::endl;
        if (!once)
        {
            double correl = mergedMtz->correlation(true);
            std::cout << "N: Merged correlation = " << correl << std::endl;
        }
        mergedMtz->description();
        
        double scale = 1000 / mergedMtz->averageIntensity();
        mergedMtz->applyScaleFactor(scale);
        
        std::cout << "Here" << std::endl;
        
        string filename = "allMerge" + i_to_str(i) + ".mtz";
        mergedMtz->writeToFile(filename.c_str(), true);
        
        MtzManager *reloadMerged = new MtzManager();
        
        reloadMerged->setFilename(filename.c_str());
        reloadMerged->loadReflections(1);
        reloadMerged->description();
        
        MtzManager::currentManager = reloadMerged;
        double reloaded = reloadMerged->correlation(true);
        std::cout << "Reloaded correlation = " << reloaded << std::endl;
        
        if (reloaded > 0.999 && i >= minimumCycles && stop)
            finished = true;
        
        if (once)
            finished = true;
        
        delete grouper;
        if (!once)
            delete reference;
        reference = exclusion ? unmergedMtz : mergedMtz;
        i++;
    }
}

// MARK: Symmetry-related reflection refinement

void MtzRefiner::refineSymmetry()
{
    readMatricesAndMtzs();
    
    std::cout << "N: Total images loaded: " << mtzManagers.size() << std::endl;
    
    std::cout << "Refining images using symmetry: " << std::endl;
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        int symmCount = mtzManagers[i]->symmetryRelatedReflectionCount();
        int paramCount = mtzManagers[i]->refinedParameterCount();
        
        if (symmCount < paramCount * 8 && paramCount > 30)
        {
            mtzManagers[i]->setRejected(true);
            continue;
        }
        
        std::cout << mtzManagers[i]->getFilename() << "\t" << symmCount << "\t" << paramCount << std::endl;
        
        mtzManagers[i]->setDefaultScoreType(ScoreTypeSymmetry);
    }
    
    refineCycle(true);
    
    MtzManager::setReference(reference);
    
    int defaultScoreInt = FileParser::getKey("DEFAULT_TARGET_FUNCTION",
                                             (int) DEFAULT_SCORE_TYPE);
    ScoreType desiredScoreType = (ScoreType) defaultScoreInt;
    
    double spotSum = 0;
    double mosaicitySum = 0;
    int count = 0;
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        double spotSize = abs(mtzManagers[i]->getSpotSize());
        double mosaicity = abs(mtzManagers[i]->getMosaicity());
        
        if (mtzManagers[i]->isRejected() == false)
        {
            spotSum += spotSize;
            mosaicitySum += mosaicity;
            count++;
        }
    }
    
    spotSum /= count;
    mosaicitySum /= count;
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        if (mtzManagers[i]->isRejected())
        {
            mtzManagers[i]->setSpotSize(spotSum);
            mtzManagers[i]->setMosaicity(mosaicitySum);
            mtzManagers[i]->setRejected(false);
        }
        
        mtzManagers[i]->setOptimisingOrientation(true);
        mtzManagers[i]->setOptimisingWavelength(true);
        
        mtzManagers[i]->setDefaultScoreType(desiredScoreType);
    }
}

// MARK: Loading data

bool MtzRefiner::loadInitialMtz(bool force)
{
    bool hasInitialMtz = FileParser::hasKey("INITIAL_MTZ");
    
    if (reference != NULL && !force)
        return true;
    
    ostringstream logged;
    
    logged << "Initial orientation matrix has "
    << (hasInitialMtz ? "" : "not ") << "been provided." << std::endl;
    
    Logger::mainLogger->addStream(&logged);
    
    if (hasInitialMtz)
    {
        
        std::string referenceFile = FileParser::getKey("INITIAL_MTZ",
                                                       string(""));
        
        reference = new MtzManager();
        
        reference->setFilename(referenceFile.c_str());
        reference->loadReflections(1);
        reference->setSigmaToUnity();
    }
    
    MtzManager::setReference(reference);
    
    return hasInitialMtz;
}


int MtzRefiner::imageMax(size_t lineCount)
{
    int end = (int)lineCount;
    
    if (imageLimit != 0)
        end = imageLimit < lineCount ? imageLimit : (int)lineCount;
    
    return end;
}

void MtzRefiner::applyParametersToImages()
{
    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    
    if (unitCell.size() == 0)
    {
        std::cout
        << "Please provide unit cell dimensions under keyword UNIT_CELL"
        << std::endl;
        exit(1);
    }
    
    int spg_num = FileParser::getKey("SPACE_GROUP", -1);
    
    if (spg_num == -1)
    {
        std::cout << "Please set space group number under keyword SPACE_GROUP"
        << std::endl;
    }
    
    double overPredSpotSize = FileParser::getKey("OVER_PRED_RLP_SIZE",
                                                 OVER_PRED_SPOT_SIZE);
    double overPredBandwidth = FileParser::getKey("OVER_PRED_BANDWIDTH",
                                                  OVER_PRED_BANDWIDTH);
    overPredBandwidth /= 2;
    
    // @TODO add orientation tolerance as flexible
    double orientationTolerance = FileParser::getKey(
                                                     "INDEXING_ORIENTATION_TOLERANCE", INDEXING_ORIENTATION_TOLERANCE);
    
    double orientationStep = FileParser::getKey("INITIAL_ORIENTATION_STEP", INITIAL_ORIENTATION_STEP);
    
    // @TODO intensityThreshold should be flexible
    double intensityThreshold = FileParser::getKey("INTENSITY_THRESHOLD",
                                                   INTENSITY_THRESHOLD);
    
    double metrologySearchSize = FileParser::getKey("METROLOGY_SEARCH_SIZE",
                                                    METROLOGY_SEARCH_SIZE);
    
    double maxIntegratedResolution = FileParser::getKey(
                                                        "MAX_INTEGRATED_RESOLUTION", MAX_INTEGRATED_RESOLUTION);
    
    bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);
    
    for (int i = 0; i < images.size(); i++)
    {
        Image *newImage = images[i];
        
        newImage->setPinPoint(true);
        
        CCP4SPG *spg = ccp4spg_load_by_standard_num(spg_num);
        
        newImage->setSpaceGroup(spg);
        newImage->setMaxResolution(maxIntegratedResolution);
        newImage->setSearchSize(metrologySearchSize);
        newImage->setIntensityThreshold(intensityThreshold);
        newImage->setUnitCell(unitCell);
        
        if (fixUnitCell)
            for (int j = 0; j < newImage->indexerCount(); j++)
                newImage->getIndexer(j)->getMatrix()->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
        
        newImage->setInitialStep(orientationStep);
      
        newImage->setTestSpotSize(overPredSpotSize);
        newImage->setTestBandwidth(overPredBandwidth);
    }
}

void MtzRefiner::readSingleImageV2(string *filename, std::vector<Image *> *newImages, int offset)
{
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    
    
    if (wavelength == 0)
    {
        std::cout
        << "Please provide initial wavelength for integration under keyword INTEGRATION_WAVELENGTH"
        << std::endl;
        exit(1);
    }
    
    if (detectorDistance == 0)
    {
        std::cout
        << "Please provide detector distance for integration under keyword DETECTOR_DISTANCE"
        << std::endl;
        exit(1);
    }
    
    const std::string contents = FileReader::get_file_contents(
                                                               filename->c_str());
    
    std::vector<std::string> imageList = FileReader::split(contents, "\nimage ");
    
    int maxThreads = FileParser::getMaxThreads();
    int end = imageMax(imageList.size());
    
    for (int i = offset; i < end; i += maxThreads)
    {
        if (imageList[i].length() == 0)
            continue;
        
        std::vector<std::string> lines = FileReader::split(imageList[i], '\n');
        
        std::vector<std::string> components = FileReader::split(lines[0], ' ');
        
        if (components.size() <= 1)
            continue;
        
        string imgName = components[1] + ".img";
        
        ostringstream logged;
        
        if (!FileReader::exists(imgName))
        {
            logged << "Skipping image " << imgName << std::endl;
            Logger::mainLogger->addStream(&logged);
            continue;
        }
        
        logged << "Loading image " << i << " (" << imgName << ")"
        << std::endl;
        Logger::mainLogger->addStream(&logged);

        double usedWavelength = wavelength;
        double usedDistance = detectorDistance;
        
        if (components.size() >= 3)
        {
            usedWavelength = atof(components[1].c_str());
            usedDistance = atof(components[2].c_str());
        }
        
        Image *newImage = new Image(imgName, wavelength,
                                    detectorDistance);
        
        for (int i = 1; i < lines.size(); i++)
        {
            std::vector<std::string> components = FileReader::split(lines[i], ' ');
            
            if (components[0] == "matrix")
            {
                // individual matrices
                double matrix[9];
                readMatrix(matrix, lines[i]);
                
                MatrixPtr newMatrix = MatrixPtr(new Matrix(matrix));
                
                double upRot = 0.75;
                newMatrix->rotate(0, 0, M_PI / 2);
                newMatrix->rotateHK(0, upRot);
                newImage->setUpIndexer(newMatrix);
            }
        }
        
        newImages->push_back(newImage);
    }
    
}

void MtzRefiner::readMatricesAndImages(string *filename)
{
    if (images.size() > 0)
        return;
    
    if (filename == NULL)
    {
        std::string aFilename = FileParser::getKey("ORIENTATION_MATRIX_LIST", string(""));
        filename = &aFilename;
        
        if (filename->length() == 0)
        {
            std::cout << "No orientation matrix list provided. Exiting now." << std::endl;
            exit(1);
        }
    }

    double version = FileParser::getKey("MATRIX_LIST_VERSION", 1.0);
    
    // thought: turn the vector concatenation into a templated function
    
    boost::thread_group threads;
    
    int maxThreads = FileParser::getMaxThreads();
    
    std::vector<std::vector<Image *> > imageSubsets;
    imageSubsets.resize(maxThreads);
    
    for (int i = 0; i < maxThreads; i++)
    {
        if (version == 1.0)
        {
            boost::thread *thr = new boost::thread(singleLoadImages, filename,
                                               &imageSubsets[i], i);
            threads.add_thread(thr);
        }
        else if (version == 2.0)
        {
            boost::thread *thr = new boost::thread(readSingleImageV2, filename,
                                                   &imageSubsets[i], i);
            threads.add_thread(thr);
        }
    }
    
    threads.join_all();
    
    
    int total = 0;
   
    for (int i = 0; i < maxThreads; i++)
    {
        total += imageSubsets[i].size();
    }
    
    mtzManagers.reserve(total);
    int lastPos = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        images.insert(images.begin() + lastPos,
                           imageSubsets[i].begin(), imageSubsets[i].end());
        lastPos += imageSubsets[i].size();
    }
    
    if (version == 2.0)
        applyParametersToImages();
}

void MtzRefiner::singleLoadImages(string *filename, std::vector<Image *> *newImages, int offset)
{
    const std::string contents = FileReader::get_file_contents(
                                                               filename->c_str());
    
    std::vector<std::string> lines = FileReader::split(contents, '\n');
    
    int spg_num = FileParser::getKey("SPACE_GROUP", -1);
    
    if (spg_num == -1)
    {
        std::cout << "Please set space group number under keyword SPACE_GROUP"
        << std::endl;
    }
    
    double overPredSpotSize = FileParser::getKey("OVER_PRED_RLP_SIZE",
                                                 OVER_PRED_SPOT_SIZE);
    double overPredBandwidth = FileParser::getKey("OVER_PRED_BANDWIDTH",
                                                  OVER_PRED_BANDWIDTH);
    overPredBandwidth /= 2;
    
    // @TODO add orientation tolerance as flexible
    double orientationTolerance = FileParser::getKey(
                                                     "INDEXING_ORIENTATION_TOLERANCE", INDEXING_ORIENTATION_TOLERANCE);
    
    double orientationStep = FileParser::getKey("INITIAL_ORIENTATION_STEP", INITIAL_ORIENTATION_STEP);
    
    // @TODO intensityThreshold should be flexible
    double intensityThreshold = FileParser::getKey("INTENSITY_THRESHOLD",
                                                   INTENSITY_THRESHOLD);
    
    double metrologySearchSize = FileParser::getKey("METROLOGY_SEARCH_SIZE",
                                                    METROLOGY_SEARCH_SIZE);
    
    double maxIntegratedResolution = FileParser::getKey(
                                                        "MAX_INTEGRATED_RESOLUTION", MAX_INTEGRATED_RESOLUTION);
    
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    
    bool fixUnitCell = FileParser::getKey("FIX_UNIT_CELL", true);
    
    if (wavelength == 0)
    {
        std::cout
        << "Please provide initial wavelength for integration under keyword INTEGRATION_WAVELENGTH"
        << std::endl;
        exit(1);
    }
    
    if (detectorDistance == 0)
    {
        std::cout
        << "Please provide detector distance for integration under keyword DETECTOR_DISTANCE"
        << std::endl;
        exit(1);
    }
    
    vector<double> unitCell = FileParser::getKey("UNIT_CELL", vector<double>());
    
    if (unitCell.size() == 0)
    {
        std::cout
        << "Please provide unit cell dimensions under keyword UNIT_CELL"
        << std::endl;
        exit(1);
    }
    
    int maxThreads = FileParser::getMaxThreads();
    
    int end = imageMax(lines.size());
    
    int skip = FileParser::getKey("IMAGE_SKIP", 0);
    
    if (skip > lines.size())
    {
        std::cout << "Image skip beyond image count" << std::endl;
        exit(1);
    }
    
    if (skip > 0)
    {
        std::cout << "Skipping " << skip << " lines" << std::endl;
    }
    
    for (int i = offset + skip; i < end + skip; i += maxThreads)
    {
        std::vector<std::string> components = FileReader::split(lines[i], ' ');
        
        string imgName = components[0] + ".img";
        if (!FileReader::exists(imgName))
        {
            continue;
        }
        
        std::cout << "Loading image " << i << " (" << imgName << ")"
        << std::endl;
        
        if (components.size() > 11)
        {
            wavelength = atof(components[10].c_str());
            detectorDistance = atof(components[11].c_str());
        }
        
        
        Image *newImage = new Image(imgName, wavelength,
                                    detectorDistance);
        
        MatrixPtr newMat = MatrixPtr();
        
        if (components.size() >= 10)
        {
            double matrix[9];
            readMatrix(matrix, lines[i]);
            newMat = MatrixPtr(new Matrix(matrix));
            if (fixUnitCell)
                newMat->changeOrientationMatrixDimensions(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
        }
        else
        {
            newMat = Matrix::matrixFromUnitCell(unitCell[0], unitCell[1], unitCell[2], unitCell[3], unitCell[4], unitCell[5]);
        }
        
        newImage->setUpIndexer(newMat);
        newImage->setPinPoint(true);
        
        CCP4SPG *spg = ccp4spg_load_by_standard_num(spg_num);
        
        newImage->setSpaceGroup(spg);
        newImage->setMaxResolution(maxIntegratedResolution);
        newImage->setSearchSize(metrologySearchSize);
        newImage->setIntensityThreshold(intensityThreshold);
        newImage->setUnitCell(unitCell);
        newImage->setInitialStep(orientationStep);
        
        // @TODO: masks should be addable through input file
        /*
         newImage->addMask(0, 840, 1765, 920);
         newImage->addMask(446, 1046, 1324, 1765);
         newImage->addMask(820, 450, 839, 473);
         */
        
        newImage->setTestSpotSize(overPredSpotSize);
        newImage->setTestBandwidth(overPredBandwidth);
        
        newImages->push_back(newImage);
    }
}

void MtzRefiner::readMatricesAndMtzs()
{
    std::string filename = FileParser::getKey("ORIENTATION_MATRIX_LIST",
                                              string(""));
    
    if (filename == "")
    {
        std::cout
        << "Orientation matrix list has not been provided for refinement. Please provide under keyword ORIENTATION_MATRIX_LIST."
        << std::endl;
        exit(1);
    }
    
    ostringstream logged;
    
    if (mtzManagers.size() > 0)
    {
        logged << "Mtzs already present; not reloading from list" << std::endl;
        Logger::mainLogger->addStream(&logged);
        return;
    }
    
    logged << "Loading MTZs from list" << std::endl;
    
    const std::string contents = FileReader::get_file_contents(
                                                               filename.c_str());
    
    std::vector<std::string> lines = FileReader::split(contents, '\n');
    
    int maxThreads = FileParser::getMaxThreads();
    
    boost::thread_group threads;
    
    std::vector<std::vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(singleThreadRead, lines,
                                               &managerSubsets[i], i);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    int total = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }
    
    mtzManagers.reserve(total);
    int lastPos = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        mtzManagers.insert(mtzManagers.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }
    
    logged << "Mtz count: " << mtzManagers.size() << std::endl;
    Logger::mainLogger->addStream(&logged);
}

void MtzRefiner::singleThreadRead(std::vector<std::string> lines,
                                  vector<MtzPtr> *mtzManagers, int offset)
{
    int maxThreads = FileParser::getMaxThreads();
    
    int end = imageMax(lines.size());
    
    int skip = FileParser::getKey("IMAGE_SKIP", 0);
    
    if (skip > lines.size())
    {
        std::cout << "Mtz skip beyond Mtz count" << std::endl;
        exit(1);
    }
    
    if (skip > 0)
    {
        std::cout << "Skipping " << skip << " lines" << std::endl;
    }
    
    for (int i = skip + offset; i < end + skip; i += maxThreads)
    {
        std::ostringstream log;
        
        std::vector<std::string> components = FileReader::split(lines[i], ' ');
        
        if (components.size() < 10)
            continue;
        
        string mtzName = components[0] + ".mtz";
        if (!FileReader::exists(mtzName))
        {
            log << "Skipping file " << mtzName << endl;
            continue;
        }
        
        log << "Loading file " << mtzName << endl;
        
        double matrix[9];
        
        readMatrix(matrix, lines[i]);
        
        MtzPtr newManager = MtzPtr(new MtzManager());
        
        newManager->setFilename(mtzName.c_str());
        newManager->setMatrix(matrix);
        newManager->loadReflections(1);
        newManager->setSigmaToUnity();
        newManager->loadParametersMap();
        
        mtzManagers->push_back(newManager);
        
        Logger::mainLogger->addStream(&log);
    }
}

void MtzRefiner::readMatrix(double (&matrix)[9], std::string line)
{
    std::vector<std::string> components = FileReader::split(line, ' ');
    
    for (int j = 1; j <= 9; j++)
    {
        string component = components[j];
        
        /* Locate the substring to replace. */
        int index = (int)component.find(string("*^"), 0);
        if (index != string::npos)
        {
            component.replace(index, 2, string("e"));
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
        MtzManager::setReference(reference);
    }
    
    
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
        
        std::cout << mtzManagers[i]->getFilename() << "\t" << correl << "\t"
        << invCorrel << std::endl;
    }
}

void MtzRefiner::merge()
{
    loadInitialMtz();
    MtzManager::setReference(this->reference);
    readMatricesAndMtzs();
    
    bool partialityRejection = FileParser::getKey("PARTIALITY_REJECTION", false);
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->excludeFromLogCorrelation();
        if (partialityRejection)
            mtzManagers[i]->excludePartialityOutliers();
    }
    
    bool denormalise = FileParser::getKey("DENORMALISE_PARTIALITY", false);
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        if (denormalise)
            mtzManagers[i]->denormaliseMillers();
    }

    correlationAndInverse(true);
    
    double correlationThreshold = FileParser::getKey("CORRELATION_THRESHOLD",
                                                     CORRELATION_THRESHOLD);
    
    vector<MtzPtr> idxOnly;
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        if ((mtzManagers[i]->getActiveAmbiguity() == 0))
            idxOnly.push_back(mtzManagers[i]);
    }
   
    
    bool mergeAnomalous = FileParser::getKey("MERGE_ANOMALOUS", false);
    int scalingInt = FileParser::getKey("SCALING_STRATEGY",
                                        (int) SCALING_STRATEGY);
    ScalingType scaling = (ScalingType) scalingInt;
    
    
    MtzGrouper *grouper = new MtzGrouper();
    grouper->setScalingType(scaling);
    grouper->setWeighting(WeightTypePartialitySigma);
    grouper->setMtzManagers(mtzManagers);
    
    grouper->setCutResolution(false);
    
    grouper->setCorrelationThreshold(correlationThreshold);
    
    MtzManager *mergedMtz = NULL;
    MtzManager *unmergedMtz = NULL;
    grouper->merge(&mergedMtz, NULL, -1, mergeAnomalous);
    mergedMtz->writeToFile("remerged.mtz");
    
    delete unmergedMtz;
    delete mergedMtz;
    delete grouper;
}

// MARK: Integrating images

void MtzRefiner::integrateImagesWrapper(MtzRefiner *object,
                                        std::vector<MtzPtr> *&mtzSubset, int offset, bool orientation)
{
    object->integrateImages(mtzSubset, offset, orientation);
}

void MtzRefiner::integrateImages(std::vector<MtzPtr> *&mtzSubset,
                                 int offset, bool orientation)
{
    int maxThreads = FileParser::getMaxThreads();
    
    bool refineDistances = FileParser::getKey("REFINE_DISTANCES", false);
    
    for (int i = offset; i < images.size(); i += maxThreads)
    {
        ostringstream logged;
        logged << "Integrating image " << i << std::endl;
        Logger::mainLogger->addStream(&logged);
        
        if (refineDistances)
            images[i]->refineDistances();

        if (orientation)
            images[i]->refineOrientations();
        
        std::vector<MtzPtr> mtzs = images[i]->currentMtzs();
        mtzSubset->insert(mtzSubset->end(), mtzs.begin(), mtzs.end());
        images[i]->dropImage();
    }
}

void MtzRefiner::integrate(bool orientation)
{
    orientation = FileParser::getKey("REFINE_ORIENTATIONS", false);
    
    loadPanels();
    
    std::string filename = FileParser::getKey("ORIENTATION_MATRIX_LIST",
                                              string(""));
    if (filename == "")
    {
        std::cout
        << "Filename list has not been provided for integration. Please provide under keyword ORIENTATION_MATRIX_LIST."
        << std::endl;
        exit(1);
    }
    
    readMatricesAndImages(&filename);
    
    int crystals = 0;
    
    for (int i = 0; i < images.size(); i++)
    {
        crystals += images[i]->indexerCount();
    }
    
    std::cout << images.size() << " images with " << crystals << " crystal orientations." << std::endl;
    
    mtzManagers.clear();
    
    int maxThreads = FileParser::getMaxThreads();
    
    boost::thread_group threads;
    std::vector<std::vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    
    for (int i = 0; i < maxThreads; i++)
    {
        boost::thread *thr = new boost::thread(integrateImagesWrapper, this,
                                               &managerSubsets[i], i, orientation);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    Panel::finaliseMillerArrays();
    
    int total = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }
    
    std::cout << "N: Total images loaded: " << total << std::endl;
    
    mtzManagers.reserve(total);
    int lastPos = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        mtzManagers.insert(mtzManagers.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }
    
    writeNewOrientations();
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->indexerCount(); j++)
        {
            std::string filename = images[i]->getFilename();
            int totalReflections = images[i]->getIndexer(j)->getTotalReflectionsWithinBandwidth();
            double lastScore = images[i]->getIndexer(j)->getLastScore();
            double bestHRot = 0;
            double bestKRot = 0;
            images[i]->getIndexer(j)->getBestRots(&bestHRot, &bestKRot);
            
            std::cout << filename << "\t" << totalReflections << "\t" << lastScore << "\t"
                      << bestHRot << "\t" << bestKRot << std::endl;
        }
    }
}


void MtzRefiner::loadPanels()
{
    std::string panelList = FileParser::getKey("PANEL_LIST", string(""));
    
    panelParser = new PanelParser(panelList);
    
    if (panelList != "" && Panel::panelCount() == 0)
    {
        std::cout << "Loading panels" << std::endl;
        panelParser->parse();
        hasPanelParser = true;
    }
    else if (panelList != "")
    {
        std::cout << "No panel list provided. Continuing regardless..."
        << std::endl;
    }
    else if (Panel::panelCount() > 0)
    {
        std::cout << "Panels already present" << std::endl;
    }
}

void MtzRefiner::findSteps()
{
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->findSteps();
    }
}

void MtzRefiner::refineDetectorGeometry()
{
    if (hasPanelParser)
    {
        for (int i = 0; i < 5; i++)
        {
            Logger::mainLogger->addString("Refining beam center and detector distance.");
            Panel::expectedBeamCentre();
            Panel::refineDetectorDistance();
            
            Panel::clearAllMillers();
            integrate(false);
        }
    }
    else
    {
        Logger::mainLogger->addString("Individual panel information has not been set up.");
    }
    
    Panel::plotAll(PlotTypeAbsolute);
}

void MtzRefiner::loadMillersIntoPanels()
{
    loadPanels();
    
    loadInitialMtz();
    MtzManager::setReference(reference);
    
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.0);
    double mmPerPixel = FileParser::getKey("MM_PER_PIXEL", MM_PER_PIXEL);

    vector<double> beam = FileParser::getKey("BEAM_CENTRE", vector<double>());
    
    if (beam.size() == 0)
    {
        beam.push_back(BEAM_CENTRE_X);
        beam.push_back(BEAM_CENTRE_Y);
    }

    for (int i = 0; i < mtzManagers.size(); i++)
    {
        MtzPtr mtz = mtzManagers[i];
        
        double hRot = mtz->getHRot();
        double kRot = mtz->getKRot();
        
        for (int j = 0; j < mtz->holderCount(); j++)
        {
            for (int k = 0; k < mtz->holder(j)->millerCount(); k++)
            {
                MillerPtr miller = mtz->holder(j)->miller(k);
                
                double x, y;
                
                miller->calculatePosition(detectorDistance, wavelength, beam[0], beam[1], mmPerPixel, hRot, kRot, &x, &y);
                
                if (miller->accepted())
                {
                    Panel::addMillerToPanelArray(miller);
                }
            }
        }
    }
    
    Panel::finaliseMillerArrays();
}

// MARK: indexing

void MtzRefiner::writeNewOrientations()
{
    ofstream newMats;
    newMats.open("new_orientations.dat");
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        MtzPtr manager = mtzManagers[i];
        
        // write out matrices etc.
        string imgFilename = manager->filenameRoot();
        newMats << "img-" << imgFilename << " ";
        
        MatrixPtr matrix = manager->getMatrix();
        
        string description = matrix->description();
        newMats << description << std::endl;
    }
    
    Logger::mainLogger->addString("Written to new_orientations.dat");
    
    newMats.close();
}

void MtzRefiner::indexImage(int offset, vector<MtzPtr> *mtzSubset)
{
    for (int i = offset; i < images.size(); i+= MAX_THREADS)
    {
        Image *newImage = images[i];
        IndexerPtr indexer = newImage->getIndexer(0);
        
        newImage->addMask(0, 840, 1765, 920);
        newImage->addMask(446, 1046, 1324, 1765);
        newImage->addMask(820, 450, 839, 473);
        newImage->addMask(1472, 789, 1655, 829);
        
        indexer->findSpots();
        
        indexer->matchMatrixToSpots();
        
        indexer->refineRoundBeamAxis();
        std::vector<MtzPtr> managers = newImage->currentMtzs();
        MtzPtr manager = managers[0];
        
        manager->description();
        manager->writeToFile(manager->getFilename());
        manager->writeToDat();
        
        mtzSubset->push_back(manager);
        
        newImage->dropImage();
    }
}

void MtzRefiner::indexImageWrapper(MtzRefiner *object, int offset, vector<MtzPtr> *mtzSubset)
{
    object->indexImage(offset, mtzSubset);
}

void MtzRefiner::index()
{
    this->readMatricesAndImages();
    loadPanels();
    std::vector<std::vector<MtzPtr>> mtzSubsets;
    mtzSubsets.resize(MAX_THREADS);
    
    boost::thread_group threads;
    
    for (int i = 0; i < MAX_THREADS; i++)
    {
        boost::thread *thr = new boost::thread(indexImageWrapper, this, i, &mtzSubsets[i]);
        threads.add_thread(thr);
    }
    
    threads.join_all();
    
    std::cout << "Images size: " << images.size() << std::endl;
    
    mtzManagers.reserve(images.size());
    
    unsigned int position = 0;
    
    for (int i = 0; i < MAX_THREADS; i++)
    {
        mtzManagers.insert(mtzManagers.begin() + position, mtzSubsets[i].begin(), mtzSubsets[i].end());
        position += mtzSubsets[i].size();
    }
    
    std::cout << "Mtzs: " << mtzManagers.size() << std::endl;
    
    writeNewOrientations();
}

// MARK: Miscellaneous

void MtzRefiner::refineMetrology()
{
    loadInitialMtz();
    
    if (!Panel::hasMillers())
    {
        loadMillersIntoPanels();
    }
    
    Panel::printToFile("new_panels.txt");
    Panel::plotAll(PlotTypeAbsolute);
}

void MtzRefiner::plotDetectorGains()
{
    bool hasMillers = Panel::hasMillers();
    
    if (!hasMillers)
    {
        loadMillersIntoPanels();
    }
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        MtzManager *ref = MtzManager::getReferenceManager();
        
        double scale = mtzManagers[i]->gradientAgainstManager(*ref);
        mtzManagers[i]->applyScaleFactor(scale);
    }
    
    Panel::checkDetectorGains();
}

void MtzRefiner::xFiles()
{
    std::string filename = FileParser::getKey("ORIENTATION_MATRIX_LIST",
                                              string(""));
    
    if (filename == "")
    {
        std::cout << "Orientation matrix list has not been provided.";
        exit(1);
    }
    
    readXFiles(filename);
    
    loadInitialMtz();
}

void MtzRefiner::readXFiles(string filename)
{
    if (mtzManagers.size() > 0)
    {
        return;
    }
    
    const std::string contents = FileReader::get_file_contents(filename.c_str());
    
    std::vector<std::string> lines = FileReader::split(contents, '\n');
    
    for (int i = 0; i < lines.size(); i++)
    {
        std::vector<std::string> filenames = FileReader::split(lines[i], ' ');
        
        XManager *xManager = new XManager();
        xManager->setFilenames(filenames);
        xManager->loadReflections(PartialityModelFixed);
        xManager->loadParametersMap();
        
        MtzPtr mtz = MtzPtr();
        mtz.reset(xManager);
        
        mtzManagers.push_back(mtz);
    }
    
    ostringstream logged;
    logged << "Mtz count: " << mtzManagers.size() << std::endl;
    Logger::mainLogger->addStream(&logged);
}

void MtzRefiner::polarisationGraph()
{
    GraphDrawer drawer = GraphDrawer(reference);
    
#ifdef MAC
    drawer.plotPolarisation(mtzManagers);
#endif
}

MtzRefiner::~MtzRefiner()
{
    delete reference;
    delete panelParser;
    mtzManagers.clear();
    
    for (int i = 0; i < images.size(); i++)
    {
        delete images[i];
    }
    
    images.clear();

    // TODO Auto-generated destructor stub
}

void MtzRefiner::removeSigmaValues()
{
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->setSigmaToUnity();
    }
    
    if (MtzManager::getReferenceManager() != NULL)
        MtzManager::getReferenceManager()->setSigmaToUnity();
}


void MtzRefiner::displayIndexingHands()
{
    ostringstream idxLogged, invLogged;
    
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        ostringstream &which = mtzManagers[i]->getActiveAmbiguity() == 0 ? idxLogged : invLogged;
        
        which << mtzManagers[i]->getFilename() << std::endl;
    }
    
    Logger::mainLogger->addString("**** Indexing hands ****\n - First indexing hand");
    Logger::mainLogger->addStream(&idxLogged);
    Logger::mainLogger->addString("- Second indexing hand");
    Logger::mainLogger->addStream(&invLogged);
    
}

void MtzRefiner::applyUnrefinedPartiality()
{
    for (int i = 0; i < mtzManagers.size(); i++)
    {
        mtzManagers[i]->applyUnrefinedPartiality();
    }
}

void MtzRefiner::refineDistances()
{
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->indexerCount(); j++)
        {
            images[i]->getIndexer(j)->refineDetectorAndWavelength(reference);
        }
    }
}

void MtzRefiner::orientationPlot()
{
#ifdef MAC
    GraphDrawer drawer = GraphDrawer(reference);
    drawer.plotOrientationStats(mtzManagers);
#endif
}