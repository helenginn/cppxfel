/*
 * MtzRefiner.h
 *
 *  Created on: 4 Jan 2015
 *      Author: helenginn
 */

#ifndef MTZREFINER_H_
#define MTZREFINER_H_

#include "parameters.h"
#include "MtzManager.h"
#include "PanelParser.h"
#include "LoggableObject.h"
#include <scitbx/mat3.h>

class IndexManager;

class MtzRefiner : public LoggableObject
{
private:
    vector<MtzPtr> mtzManagers;
	MtzPtr reference;
    MtzPtr referencePtr;
	vector<ImagePtr> images;
    static bool hasPanelParser;
    PanelParser *panelParser;
    static int imageLimit;
    static int imageMax(size_t lineCount);
    static void singleLoadImages(std::string *filename, vector<ImagePtr> *newImages, int offset);
    static void readSingleImageV2(std::string *filename, vector<ImagePtr> *newImages, vector<MtzPtr> *newMtzs, int offset, bool v3 = false);
    static void findSpotsThread(MtzRefiner *me, int offset);
    void readFromHdf5(std::vector<ImagePtr> *newImages);
    bool readRefinedMtzs;
IndexManager *indexManager;
    void applyParametersToImages();
    static int cycleNum;
    bool hasRefined;
    int maxThreads;
    bool isPython;
    static int imageSkip(size_t totalCount);
    static void radialAverageThread(MtzRefiner *me, int offset);
    static void integrateSpotsThread(MtzRefiner *me, int offset);
    Hdf5ManagerProcessingPtr hdf5ProcessingPtr;
    void readDataFromOrientationMatrixList(std::string *filename, bool areImages, std::vector<ImagePtr> *targetImages);
public:
	MtzRefiner();
	virtual ~MtzRefiner();

    void index();
    void indexFromScratch();
    void powderPattern();
    void hitAnalysis();
	bool loadInitialMtz(bool force = false);
    void indexingParameterAnalysis();
    
	void cycle();
	void cycleThread(int offset);
	static void cycleThreadWrapper(MtzRefiner *object, int offset);

    void refineSymmetry();
	void refine();
	void refineCycle(bool once = false);
	void readMatricesAndMtzs();
    void refineMetrology();
    void initialMerge();
    void orientationPlot();
    void applyUnrefinedPartiality();
    void loadImageFiles();
    void findSpots();
    
    void loadPanels(bool mustFail = true);
	void integrate();
    void integrationSummary();
	static void integrateImagesWrapper(MtzRefiner *object,
			vector<MtzPtr> *&mtzSubset, int offset, bool orientation);
	void integrateImages(vector<MtzPtr> *&mtzSubset, int offset, bool orientation);
	void readMatricesAndImages(std::string *filename = NULL, bool areImages = true, std::vector<ImagePtr> *targetImages = NULL);
    void combineLists();
    
	static void readMatrix(double (&matrix)[9], std::string line);
	static void singleThreadRead(vector<std::string> lines,
			vector<MtzPtr> *mtzManagers, int offset);
	void merge(bool mergeOnly = false);
    void correlationAndInverse(bool shouldFlip = false);
    void refreshCurrentPartialities();
    void maximumImage();
    static void maximumImageThread(MtzRefiner *me, ImagePtr maxImage, int offset);
    
    static int getCycleNum()
    {
        return cycleNum;
    }
    
    vector<MtzPtr> getMtzManagers()
    {
        return mtzManagers;
    }
    
    bool isFromPython()
    {
        return isPython;
    }
    
    void setFromPython(bool newValue)
    {
        isPython = newValue;
    }
    
    void setupFreeMillers();
    void refineDistances();
    void polarisationGraph();
    void displayIndexingHands();
    void findSteps();
    
    void writeAllNewOrientations();
    void writeNewOrientations(bool includeRots = false, bool detailed = false);
    void removeSigmaValues();
    void radialAverage();
    void integrateSpots();
    void linearScaling();
    
    void plotIntensities();
    void plotIntegrationWindows();
    
    void writePNGs(int total = 0);
};

#endif /* MTZREFINER_H_ */
