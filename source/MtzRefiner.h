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
	MtzManager *reference;
	vector<ImagePtr> images;
    static bool hasPanelParser;
    PanelParser *panelParser;
    static int imageLimit;
    static int imageMax(size_t lineCount);
    static void singleLoadImages(std::string *filename, vector<ImagePtr> *newImages, int offset);
    static void readSingleImageV2(std::string *filename, vector<ImagePtr> *newImages, vector<MtzPtr> *newMtzs, int offset);
    IndexManager *indexManager;
    void applyParametersToImages();
    static int cycleNum;
    bool hasRefined;
    bool isPython;
    static int imageSkip(size_t totalCount);
    static void radialAverageThread(MtzRefiner *me, int offset);
    static void integrateSpotsThread(MtzRefiner *me, int offset);
public:
	MtzRefiner();
	virtual ~MtzRefiner();

    void index();
    void indexFromScratch();
    void powderPattern();
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
    
    void loadPanels();
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
    
    void writeNewOrientations(bool includeRots = false, bool detailed = false);
    void removeSigmaValues();
    void radialAverage();
    void integrateSpots();
    
};

#endif /* MTZREFINER_H_ */
