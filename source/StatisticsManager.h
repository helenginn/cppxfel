#ifndef statistics
#define statistics


#include <string>
#include <iostream>
#include "MtzManager.h"

struct Partial
{
    MillerPtr miller;
        double partiality;
        double percentage;
        double resolution;
        double wavelength;
};

class StatisticsManager
{
private:

public:
        StatisticsManager(void);
        ~StatisticsManager(void);

    vector<MtzPtr> mtzs;
        void loadFiles(char **filenames, int filenum, int partiality);

        void setMtzs(vector<MtzPtr> mtzs);

        static double midPointBetweenResolutions(double minD, double maxD);
        static void generateResolutionBins(double minD, double maxD, int binCount,
                        vector<double> *bins);
        static void convertResolutions(double lowAngstroms, double highAngstroms,
                        double *lowReciprocal, double *highReciprocal);

        double cc_through_origin(int num1, int num2, int silent, int inverted,
                        int *hits);
        double cc_through_origin(int num1, int num2, int silent, int inverted,
                        int *hits, double lowResolution, double highResolution, bool log);

        static double cc_pearson(MtzManager *shot1, MtzManager *shot2, int silent = 1,
            int *hits = NULL, double *multiplicity = NULL, double lowResolution = 0,
                        double highResolution = 0, bool log = false, bool freeOnly = false);
        double cc_pearson(int num1, int num2, int silent, int *hits,
                        double *multiplicity, double lowResolution = 0, double highResolution =
                                        0, bool log = false, bool freeOnly = false);

        static double r_factor(RFactorType rFactor, MtzManager *shot1, int *hits,
                        double *multiplicity, double lowResolution, double highResolution, bool freeOnly = false);

        static double r_split(MtzManager *shot1, MtzManager *shot2, int silent,
                        int *hits, double *multiplicity, double lowResolution,
                        double highResolution, bool log, bool freeOnly = false);

        int mtz_num;
};

#endif // statistics
