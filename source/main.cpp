#include <iostream>
#include <string>

#include "InputFileParser.h"
#include "MtzManager.h"
#include "StatisticsManager.h"
#include "misc.h"
#include "lbfgs_cluster.h"
#include "GraphDrawer.h"
#include "Wiki.h"
#include "Index.h"
#include "Logger.h"



int main(int argc, char *argv[])
{
    time_t startcputime;
    time(&startcputime);
    
	if (argc == 1)
	{
		std::cout << "No tasks selected. Try -i {input_file}." << std::endl;
		exit(1);
	}

	if (strcmp(argv[1], "-wiki") == 0)
	{
		if (argc <= 2)
		{
			std::cout << "arguments: -wiki <logfile>" << std::endl;
			exit(1);
		}

        Wiki wiki = Wiki(std::string(argv[2]));
		wiki.process();

		exit(1);
	}
    
#ifdef MAC
//	std::cout << "Original: " << getenv("SYMINFO") << std::endl;

	setenv("SYMINFO", "/Applications/ccp4-6.4.0/lib/data/syminfo.lib", 1);

	std::cout << getenv("SYMINFO") << std::endl;
#endif
    
	std::cout << "Welcome to Helen's XFEL tasks" << std::endl;
    
	if (strcmp(argv[1], "-i") == 0)
	{
		if (argc < 3)
		{
			std::cout << "arguments: -i <input_script>" << std::endl;
			exit(1);
		}

        InputFileParser *parser = new InputFileParser(std::string(argv[2]));
		parser->parse();
        
        delete parser;
	}
    else
    {
        Logger::mainLogger = LoggerPtr(new Logger());
        boost::thread thr = boost::thread(Logger::awaitPrintingWrapper, Logger::mainLogger);
    }
    
	if (strcmp(argv[1], "-index") == 0)
	{
		index(argv, argc);
	}

	if (strcmp(argv[1], "-allcc") == 0)
	{
		if (argc <= 2)
		{
			std::cout << "arguments: -allcc <file1> <file2> ... <filen>." << std::endl;
			exit(1);
		}

		Lbfgs_Cluster *lbfgs = new Lbfgs_Cluster();

		(*lbfgs).initialise_cluster_lbfgs(&argv[2], argc - 2, NULL);

		delete lbfgs;
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
		mtz->loadReflections(1);

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
		mtz1->loadReflections(1);

		MtzManager *mtz2 = new MtzManager();
		mtz2->setFilename(std::string(argv[3]));
		mtz2->loadReflections(1);

		if (inverted)
            mtz1->setActiveAmbiguity(1);

		if (rsplit)
		{
			mtz1->rSplitWithManager(mtz2, 1, 0, lowRes, highRes, bins);
		}
		else
		{
			mtz1->correlationWithManager(mtz2, 1, 0, lowRes, highRes, bins);

		}

		delete mtz1;
		delete mtz2;

		exit(1);
	}

	if (strcmp(argv[1], "-gradscaling") == 0)
	{
		if (argc <= 3)
		{
			std::cout << "arguments: -gradscaling <ref> <file2> ... <filen>."
					<< std::endl;
			exit(1);
		}

		MtzManager *reference = new MtzManager();
		reference->setFilename(argv[2]);
		reference->loadReflections(1);

		for (int i = 3; i < argc; i++)
		{
			MtzManager *image = new MtzManager();
			image->setFilename(argv[i]);
			image->loadReflections(1);

			double gradientOld = image->gradientAgainstManager(*reference);
			double gradientBefore = image->minimizeRFactor(reference);

			std::cout << image->getFilename() << "\t" << gradientOld << "\t"
					<< gradientBefore << std::endl;
		}
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
		mtz->loadReflections(1);
        mtz->setActiveAmbiguity(1);

		mtz->writeToFile(std::string("inv-") + argv[2], false, false, true);
	}

	if (strcmp(argv[1], "-stats") == 0)
	{
		if (argc < 3)
		{
			std::cout << "arguments: -stats <filein> <threshold>." << std::endl;
			exit(1);
		}

		if (argc >= 3)
		{
			StatisticsManager stats;
			stats.loadFiles(&argv[2], 1, 0);
            
            double threshold = -100;
            int h = 0; int k = 0; int l = 0;
            
            if (argc >= 4)
            {
                threshold = atof(argv[3]);
            }
            
            if (argc >= 7)
            {
                h = atoi(argv[4]);
                k = atoi(argv[5]);
                l = atoi(argv[6]);
            }

   //         stats.partialityStats(0, threshold, h, k, l);

#ifdef MAC
            GraphDrawer drawer = GraphDrawer(&*stats.mtzs[0]);
            drawer.plotPartialityStats();
#endif
		}
	}

#ifdef MAC
    if (strcmp(argv[1], "-partimg") == 0)
	{
		if (argc <= 3)
		{
			std::cout << "arguments: -partimg <ref> <filein>." << std::endl;
			exit(1);
		}

		MtzManager *reference = new MtzManager();
		reference->setFilename(argv[2]);
		reference->loadReflections(1);
		MtzManager::setReference(reference);

		for (int i = 3; i < argc; i++)
		{
			MtzManager *mtz = new MtzManager();
			mtz->setFilename(argv[i]);
            std::cout << argv[i] << std::endl;
			mtz->loadReflections(PartialityModelNone);

            vector<Reflection *>refReflections, imageReflections;
            
			GraphDrawer graph = GraphDrawer(mtz);

			graph.partialityPlot("partiality", GraphMap());

			delete mtz;
		}

		delete reference;
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
		reference->loadReflections(1);
		MtzManager::setReference(reference);

		for (int i = 3; i < argc; i++)
		{
			MtzManager *mtz = new MtzManager();
			mtz->setFilename(argv[i]);
			mtz->loadReflections(1);

			managers.push_back(mtz);
		}
        
        
        for (int i = 0; i < managers.size(); i++)
        {
            managers[i]->denormaliseMillers();
        }

		GraphDrawer graph = GraphDrawer(reference);

		graph.resolutionStatsPlot(managers, "intensity_bins", GraphMap(), true,
				false);
		graph.resolutionStatsPlot(managers, "intensity_bins_2", GraphMap(),
				true, true);
		graph.resolutionStatsPlot(managers);
        graph.bFactorPlot(managers);
        
		for (int i = 0; i < managers.size(); i++)
			delete managers[i];

	}

	if (strcmp(argv[1], "-ccplot") == 0)
	{
		MtzManager *reference = new MtzManager();
		reference->setFilename(argv[2]);
		reference->loadReflections(1);
		MtzManager::setReference(reference);

		for (int i = 3; i < argc; i++)
		{
			MtzManager *image = new MtzManager();
			image->setFilename(argv[i]);
			image->loadReflections(1);

			GraphDrawer graph = GraphDrawer(image);
			graph.correlationPlot("correl");

			delete image;
		}

		delete reference;
	}
#endif
    
    time_t endcputime;
    time(&endcputime);
    
    clock_t difference = endcputime - startcputime;
    double seconds = difference;
    
    int finalSeconds = (int) seconds % 60;
    int minutes = seconds / 60;
    
    std::ostringstream logged;
    logged << "N: Total time: " << minutes << " minutes, "
    << finalSeconds << " seconds (" << seconds << " seconds)." << std::endl;

	logged << "Done" << std::endl;
    Logger::mainLogger->addStream(&logged);
    
    sleep(2);
}
