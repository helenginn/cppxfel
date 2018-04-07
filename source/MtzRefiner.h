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
#include "LoggableObject.h"
#include <map>

typedef std::map<int, std::vector<MtzPtr> > BinList;

class IndexManager;
class ProfileFit;

class MtzRefiner : public LoggableObject
{
private:
        MtzPtr reference;
    MtzPtr referencePtr;
        vector<ImagePtr> images;
    static int imageLimit;
    static int imageMax(size_t lineCount);
    static void readSingleImageV2(std::string *filename, vector<ImagePtr> *newImages, vector<MtzPtr> *newMtzs, int offset, bool v3 = false, MtzRefiner *me = NULL);
    static void findSpotsThread(MtzRefiner *me, int offset);
    void readFromHdf5(std::vector<ImagePtr> *newImages);
    bool readRefinedMtzs;
        std::vector<MtzPtr> getAllMtzs();
        IndexManager *indexManager;
    ProfileFit *profileFit;
    static int cycleNum;
    bool hasRefined;
    int maxThreads;
    bool isPython;
    static int imageSkip(size_t totalCount);
    static void radialAverageThread(MtzRefiner *me, int offset);
    static void integrateSpotsThread(MtzRefiner *me, int offset);
    Hdf5ManagerProcessingPtr hdf5ProcessingPtr;
    void readDataFromOrientationMatrixList(std::string *filename, bool areImages, std::vector<ImagePtr> *targetImages);
        void redumpBins();

        BinList binList;
public:
        MtzRefiner();
        virtual ~MtzRefiner();

    
    bool scaleBFactors;
    void profileFitter();
    void refineAllBFactors();
    void smoothenRefinedAllBFactors();
    
    //bool fileSort(MtzPtr a, MtzPtr b);
    //void smoothenScale();
    void index();
    void powderPattern();
        bool loadInitialMtz(bool force = false);

        void cycle();
        void cycleThread(int offset);
        static void cycleThreadWrapper(MtzRefiner *object, int offset);

    void refine();
        void refineCycle(bool once = false);
        void readMatricesAndMtzs();

        void refineMetrology(bool global);
    void reportMetrology();
        void flattenDetector();
    inline long imageNumber(const std::string& inFile);
    void initialMerge();
    void orientationPlot();
    void applyUnrefinedPartiality();
    void loadImageFiles();
    void findSpots();
        void differenceCorrelation();

    void loadPanels(bool mustFail = true);
        void integrate();
    static void fakeSpotsThread(std::vector<ImagePtr> *images, int offset);
    void fakeSpots();
    void integrationSummary();
        static void integrateImagesWrapper(MtzRefiner *object, vector<MtzPtr> *&mtzSubset, int offset);
        void integrateImages(vector<MtzPtr> *&mtzSubset, int offset);
        void readMatricesAndImages(std::string *filename = NULL, bool areImages = true, std::vector<ImagePtr> *targetImages = NULL);
        void refineUnitCell();

        static void readMatrix(double (&matrix)[9], std::string line);
        void merge(int cycle = -2);
    void correlationAndInverse(bool shouldFlip = false);
    void refreshCurrentPartialities();
    void maximumImage();
    static void maximumImageThread(MtzRefiner *me, ImagePtr maxImage, int offset);

    static int getCycleNum()
    {
        return cycleNum;
    }

    bool isFromPython()
    {
        return isPython;
    }

    void setFromPython(bool newValue)
    {
        isPython = newValue;
    }

    void polarisationGraph();
    void displayIndexingHands();
        void plotPixelValueVsFiducial();

    void writeAllNewOrientations();
    void writeNewOrientations(bool includeRots = false, bool detailed = false);
    void integrateSpots();
    void linearScaling();
    //void plotProfile();
    void plotIntensities();
    void plotIntegrationWindows();

    void imageToDetectorMap();
    void writePNGs(int total = 0);
    void takeTwoPNG();
};

#endif /* MTZREFINER_H_ */

