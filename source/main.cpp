//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <string>

#include "InputFileParser.h"
#include "MtzManager.h"
#include "StatisticsManager.h"
#include "misc.h"
#include "GraphDrawer.h"
#include "Logger.h"
#include <fstream>
#include <unistd.h>
#include "Hdf5ManagerCheetahSacla.h"
#include <execinfo.h>
#include <signal.h>

void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, (int)size, STDERR_FILENO);
    exit(1);
}


void finishJobNotification(int argc, char *argv[], int minutes)
{
    const char *jobNotificationFileStr = std::getenv("JOB_NOTIFICATION_FILE");

    if (jobNotificationFileStr == NULL)
    {
        return;
    }

    std::ostringstream command;
    command << "cppxfel.run ";
    for (int i = 1; i < argc; i++)
    {
        command << argv[i] << " ";
    }

    std::string workingDirectory;

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
    {
        workingDirectory = cwd;
    }

    std::ostringstream notification;
    notification << "osascript -e 'display notification \"" << workingDirectory << "\" with title \"" << command.str() << "\" subtitle \"" << minutes << " minutes to complete\" sound name \"Glass\"'" << std::endl;

    std::ofstream jobNotificationFile;
    jobNotificationFile.open(jobNotificationFileStr, std::ofstream::out | std::ofstream::app);
    jobNotificationFile << notification.str();
    jobNotificationFile.close();

    std::cout << "Job notification posted." << std::endl;
}

void new_main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    new_main(argc, argv);
}

void new_main(int argc, char *argv[])
{
        std::cout.sync_with_stdio(false);

    signal(SIGSEGV, handler);
    signal(SIGABRT, handler);

    time_t startcputime;
    time(&startcputime);

        if (argc == 1)
        {
        std::cout << "Welcome to cppxfel version 2.0!" << std::endl;
        std::cout << "Please refer to & cite paper (Ginn et al., J. Appl. Cryst. (2016). 49, 1065-1072)" << std::endl << std::endl;

                std::cout << "Specify input file, e.g.:" << std::endl;
        std::cout << "\tcppxfel.run -i index.txt" << std::endl;
        std::cout << "\tcppxfel.run -i integrate.txt" << std::endl;
        std::cout << "\tcppxfel.run -i refine.txt" << std::endl;
        std::cout << "\tcppxfel.run -i merge.txt" << std::endl << std::endl;;
        std::cout << "Other functions for assessing data quality:" << std::endl << std::endl;

        std::cout << "Correlation between two MTZ files:" << std::endl << std::endl;
        std::cout << "\tcppxfel.run -cc firstFile.mtz secondFile.mtz [ambiguity] [lowRes] [highRes] [bins]" << std::endl << std::endl;
        std::cout << "ambiguity: 0, 1, 2 or 3 - will compare different indexing solutions where the Bravais lattice symmetry is higher than that of the point group for certain space groups. Default 0" << std::endl;
        std::cout << "lowRes and highRes: set to resolution in Angstroms to bound the results, or set to 0 to take lowest/highest resolution data. Default 0, 0" << std::endl;
        std::cout << "bins: number of bins to report correlation statistics. Default 20." << std::endl << std::endl;

        std::cout << "Partiality CSV files:" << std::endl << std::endl;
        std::cout << "\tcppxfel.run -partiality reference.mtz ref-img-shot-number.mtz [highRes]" << std::endl << std::endl;
        std::cout << "highRes: - highest resolution reflection to report results on. Defaults to the edge of the image." << std::endl;
                std::cout << "This outputs ref-img-shot-number.mtz_partiality_[r].png where r denotes the resolution range, and also the relevant CSV with the .csv file extension." << std::endl << std::endl;;


        std::cout << "Merging statistics:" << std::endl << std::endl;
        std::cout << "\tcppxfel.run -rpim unmerged_file.mtz [lowRes] [highRes] [bins]" << std::endl;
        std::cout << "\tcppxfel.run -rmeas unmerged_file.mtz [lowRes] [highRes] [bins]" << std::endl;
        std::cout << "\tcppxfel.run -rmerge unmerged_file.mtz [lowRes] [highRes] [bins]" << std::endl << std::endl;
        std::cout << "lowRes and highRes: set to resolution in Angstroms to bound the results, or set to 0 to take lowest/highest resolution data. Default 0, 0" << std::endl;
        std::cout << "bins: number of bins to report correlation statistics. Default 20." << std::endl << std::endl;

        std::cout << "General help for most of the parameters which are specified in input files can be found using --help or -h, for example:" << std::endl;
        std::cout << "\tcppxfel.run --help                       (full list)" << std::endl;
        std::cout << "\tcppxfel.run --help verbosity_level" << std::endl;
        std::cout << "\tcppxfel.run --help intensity_threshold" << std::endl;

        std::cout << "\nNow testing your CCP4 library installation..." << std::endl << std::endl;

        CSym::CCP4SPG *spg = ccp4spg_load_by_ccp4_num(1);
        if (!spg)
        {
            std::cout << std::endl << "*****************************************************" << std::endl;
            std::cout <<              "************* WARNING! VERY IMPORTANT!! *************" << std::endl;
            std::cout <<              "*****************************************************" << std::endl;
            std::cout << std::endl << "I can't seem to load space group P1!" << std::endl;
            std::cout << "If CCP4 is installed, it may need sourcing." << std::endl;
            std::cout << "Please make sure the value of the environment variable $SYMINFO is set correctly." << std::endl;
            std::cout << "Something like (bash): export SYMINFO=/path/to/ccp4-vX.X.X/lib/data/syminfo.lib" << std::endl;
            std::cout << "Something like (csh): setenv SYMINFO /path/to/ccp4-vX.X.X/lib/data/syminfo.lib" << std::endl;
            std::cout << std::endl << "(... look carefully, those paths are almost certainly wrong!)" << std::endl;
        }
        else
        {
            std::cout << "I found space group P1! Everything is fine." << std::endl;
        }


        exit(1);
        }

    Logger::mainLogger = LoggerPtr(new Logger());
    boost::thread thr = boost::thread(Logger::awaitPrintingWrapper, Logger::mainLogger);

    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0 )
    {
        FileParser parser;

        if (argc == 2)
        {
            FileParser::printAllCommands();
        }
        else
        {
            std::string whichCommand = argv[2];
            FileParser::printCommandInfo(whichCommand);
        }
    }

        std::cout << "Welcome to cppxfel!" << std::endl;

        if (strcmp(argv[1], "-i") == 0 || strcmp(argv[1], "-dry") == 0)
        {
        bool dry = (strcmp(argv[1], "-dry") == 0);

                if (argc < 3)
                {
                        std::cout << "arguments: -i <input_script>" << std::endl;
                        exit(1);
                }

        std::vector<std::string> extras;

        for (int i = 3; i < argc; i++)
        {
            extras.push_back(std::string(argv[i]));
        }

        InputFileParser *parser = new InputFileParser(std::string(argv[2]), extras);
        parser->setDry(dry);
                parser->parse(false);

        delete parser;
        }

    if (strcmp(argv[1], "-b") == 0)
    {
        float bFactor = atof(argv[2]);

        MtzManager *mtz1 = new MtzManager();
        mtz1->setFilename(std::string(argv[3]));
        mtz1->loadReflections();

        mtz1->applyBFactor(bFactor);

        mtz1->writeToFile("b-" + std::string(argv[3]));
    }

        if (strcmp(argv[1], "-rmerge") == 0 || strcmp(argv[1], "-rpim") == 0
                        || strcmp(argv[1], "-rmeas") == 0)
        {
                if (argc < 3)
                {
                        std::cout << "arguments: -r{merge} <file1> [lowRes] [highRes] [bins]."
                                        << std::endl;
                        exit(1);
                }

                RFactorType rFactor = RFactorTypeMerge;

                if (strcmp(argv[1], "-rpim") == 0)
                {
                        rFactor = RFactorTypePim;
                }
                else if (strcmp(argv[1], "-rmeas") == 0)
                {
                        rFactor = RFactorTypeMeas;
                }

                double highRes = 0;
                double lowRes = 0;
                int bins = 20;

                if (argc > 5)
                {
                        lowRes = atof(argv[4]);
                        highRes = atof(argv[5]);
                }
                if (argc > 6)
                {
                        bins = atoi(argv[6]);
                }

        MtzManager *mtz = new MtzManager();
                mtz->setFilename(std::string(argv[2]));
                mtz->loadReflections();

                mtz->rFactorWithManager(rFactor, false, false, lowRes, highRes, bins);
        }

        if (strcmp(argv[1], "-cc") == 0 || strcmp(argv[1], "-rsplit") == 0)
        {
                if (argc < 4)
                {
                        std::cout
                                        << "arguments: -cc <file1> <file2> (0 / 1) [lowRes] [highRes] [bins]."
                                        << std::endl;
                        exit(1);
                }

                bool rsplit = (strcmp(argv[1], "-rsplit") == 0);

                int inverted = 0;
                double highRes = 0;
                double lowRes = 0;
                int bins = 20;

                if (argc > 4)
                {
                        inverted = atoi(argv[4]);
                }
                if (argc > 6)
                {
                        lowRes = atof(argv[5]);
                        highRes = atof(argv[6]);
                }
                if (argc > 7)
                {
                        bins = atoi(argv[7]);
                }

                MtzManager *mtz1 = new MtzManager();
                mtz1->setFilename(std::string(argv[2]));
                mtz1->loadReflections();

                MtzManager *mtz2 = new MtzManager();
                mtz2->setFilename(std::string(argv[3]));
                mtz2->loadReflections();

                if (inverted)
            mtz1->setActiveAmbiguity(1);

                if (rsplit)
                {
                        mtz1->rSplitWithManager(mtz2, 1, 0, lowRes, highRes, bins);
                }
                else
                {
                        mtz1->correlationWithManager(mtz2, 1, 0, lowRes, highRes, bins, NULL, true);

                }

                delete mtz1;
                delete mtz2;

        std::ostringstream logged;
        logged << "Done" << std::endl;
        LoggableObject::staticLogAndExit(logged, "DONE");
        }

        if (strcmp(argv[1], "-inv") == 0)
        {
                if (argc <= 2)
                {
                        std::cout << "arguments: -inv <file>." << std::endl;
                        exit(1);
                }

        MtzManager *mtz = new MtzManager();
                mtz->setFilename(std::string(argv[2]));
                mtz->loadReflections();
        mtz->setActiveAmbiguity(1);

                mtz->writeToFile(std::string("inv-") + argv[2], true, true);
        }

    if (strcmp(argv[1], "-intensities") == 0)
    {
        if (argc <= 5)
        {
            std::cout << "arguments: -intensities h k l <file1> {<file2> ...}." << std::endl;
            exit(1);
        }

        int h = atoi(argv[2]);
        int k = atoi(argv[3]);
        int l = atoi(argv[4]);

        if (h == 0 && k == 0 && l == 0)
        {
            h = (double)(rand() / RAND_MAX) * 40;
            k = (double)rand() / RAND_MAX * (40 - h) + h;
            l = (double)rand() / RAND_MAX * (40 - k) + k;
        }

        std::cout << "Plotting intensities for (" << h << ", " << k << ", " << l << ")" << std::endl;


        std::vector<MtzPtr> mtzs;

        for (int i = 5; i < argc; i++)
        {
            MtzPtr mtz = MtzPtr(new MtzManager());
            mtz->setFilename(argv[i]);
            mtz->loadReflections();

            mtzs.push_back(mtz);
        }

        GraphDrawer drawer = GraphDrawer(&*mtzs[0]);
        drawer.plotReflectionFromMtzs(mtzs, h, k, l);
        }

    if (strcmp(argv[1], "-partiality") == 0)
        {
                if (argc <= 3)
                {
                        std::cout << "arguments: -partiality <ref> <filein> {<maxres>}." << std::endl;
                        exit(1);
                }

                MtzManager *reference = new MtzManager();
                reference->setFilename(argv[2]);
                reference->loadReflections();
                MtzManager::setReference(reference);

        double maxRes = 0;

        if (argc > 4)
        {
            maxRes = atof(argv[4]);
        }


        MtzPtr mtz = MtzPtr(new MtzManager());
        mtz->setFilename(argv[3]);
        std::cout << "Partiality plot for " << argv[3] << std::endl;
        mtz->loadReflections();

        vector<Reflection *>refReflections, imageReflections;

        GraphDrawer graph = GraphDrawer(&*mtz);

        graph.partialityPNG(mtz, maxRes);

                delete reference;
        }

        if (strcmp(argv[1], "-shellscale") == 0)
        {
                MtzManager *reference = new MtzManager();
                reference->setFilename(argv[2]);
                reference->loadReflections();
                MtzManager::setReference(reference);

                MtzManager *change = new MtzManager();
                change->setFilename(argv[3]);
                change->loadReflections();
                change->description();

                change->applyScaleFactorsForBins(40);
                change->description();

                change->writeToFile("shsc-" + change->getFilename(), true, false, true);

        }

    if (strcmp(argv[1], "-bfac") == 0)
    {
        MtzManager *reference = new MtzManager();
        reference->setFilename(argv[2]);
        reference->loadReflections();

        double bFactor = atof(argv[3]);
		std::cout << "Calculated B factor: " << bFactor << std::endl;

        reference->applyBFactor(bFactor);
        reference->writeToFile("bfac-" + reference->getFilename(), true);

    }

        if (strcmp(argv[1], "-bfactor") == 0)
        {
                if (argc <= 2)
                {
                        std::cout << "arguments: -bfactor <ref> <file1> {<file2> ...}." << std::endl;
                        exit(1);
                }

                vector<MtzManager *> managers = vector<MtzManager *>();

                MtzManager *reference = new MtzManager();

                reference->setFilename(argv[2]);
                reference->loadReflections();
                MtzManager::setReference(reference);
  //      FileParser::setKey("REFINE_B_FACTOR", true);

                for (int i = 3; i < argc; i++)
                {
                        MtzManager *mtz = new MtzManager();
                        mtz->setFilename(argv[i]);
                        mtz->loadReflections();
            managers.push_back(mtz);
                }

                GraphDrawer graph = GraphDrawer(reference);
        graph.resolutionStatsCSV(managers);


                for (int i = 0; i < managers.size(); i++)
                        delete managers[i];

        }

    time_t endcputime;
    time(&endcputime);

    clock_t difference = endcputime - startcputime;
    double seconds = difference;

    int finalSeconds = (int) seconds % 60;
    int minutes = seconds / 60;

    std::ostringstream logged;
    logged << "N: Total time: " << minutes << " minutes, "
    << finalSeconds << " seconds (" << seconds << " seconds)." << std::endl;

    if (strcmp(argv[1], "-i") == 0)
        finishJobNotification(argc, argv, minutes);

    Hdf5ManagerCheetahSacla::closeHdf5Files();

    logged << "Done" << std::endl;
    LoggableObject::staticLogAndExit(logged, "DONE");
}
