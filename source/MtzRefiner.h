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

class MtzRefiner
{
private:
	vector<MtzPtr> mtzManagers;
	MtzManager *reference;
	vector<Image *> images;
    static bool hasPanelParser;
    PanelParser *panelParser;
    static int imageLimit;
    static int imageMax(size_t lineCount);
    static void singleLoadImages(string *filename, std::vector<Image *> *newImages, int offset);
    static void readSingleImageV2(string *filename, std::vector<Image *> *newImages, int offset);
    void indexImage(int offset, vector<MtzPtr> *mtzSubset);
    static void indexImageWrapper(MtzRefiner *object, int offset, vector<MtzPtr> *mtzSubset);
    void applyParametersToImages();
    static int cycleNum;
public:
	MtzRefiner();
	virtual ~MtzRefiner();

    void index();
	bool loadInitialMtz(bool force = false);

	void cycle();
	void cycleThread(int offset);
	static void cycleThreadWrapper(MtzRefiner *object, int offset);

    void refineSymmetry();
	void refine();
	void refineCycle(bool once = false);
	void readMatricesAndMtzs();
    void refineDetectorGeometry();
    void refineMetrology();
    void loadMillersIntoPanels();
    void plotDetectorGains();
    void initialMerge();
    void orientationPlot();
    void applyUnrefinedPartiality();

    void loadPanels();
	void integrate(bool orientation);
	static void integrateImagesWrapper(MtzRefiner *object,
			std::vector<MtzPtr> *&mtzSubset, int offset, bool orientation);
	void integrateImages(std::vector<MtzPtr> *&mtzSubset, int offset, bool orientation);
	void readMatricesAndImages(string *filename = NULL);

	static void readMatrix(double (&matrix)[9], std::string line);
	static void singleThreadRead(std::vector<std::string> lines,
			std::vector<MtzPtr> *mtzManagers, int offset);
	void merge();
    void correlationAndInverse(bool shouldFlip = false);
    
    void refineDistances();
    void polarisationGraph();
    void displayIndexingHands();
    void findSteps();
    
    void writeNewOrientations(bool includeRots = false);
    void removeSigmaValues();
    void readXFiles(string filename);
    void xFiles();
};

#endif /* MTZREFINER_H_ */
