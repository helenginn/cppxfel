/*
 * Index.cpp
 *
 *  Created on: 26 Dec 2014
 *      Author: helenginn
 */

#include "Index.h"
#include "definitions.h"
#include "Matrix.h"
#include <string>
#include "Image.h"
#include "Vector.h"
#include "parameters.h"
#include "headers/csymlib.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include <boost/thread/thread.hpp>
#include "FileReader.h"
#include <vector>
#include "FileParser.h"
#include "MtzManager.h"

Index::Index()
{
	// TODO Auto-generated constructor stub

}

Index::~Index()
{
	// TODO Auto-generated destructor stub
}



void findSpots(char **argv, int argc, int offset, vector<Image *>**images)
{
    int maxThreads = FileParser::getMaxThreads();
    
	for (int i = 2 + offset; i < argc; i += maxThreads)
	{
		double wavelength = 1.304;
		double detectorDistance = 2520;

		MatrixPtr newMat = Matrix::matrixFromUnitCell(295.14, 295.14, 678.23, 90,
				90, 120);
		std::string imgName = std::string(argv[i]);
		Image *newImage = new Image(imgName, wavelength,
				detectorDistance);
		(*images)->push_back(newImage);
		newImage->setUpIndexer(newMat);
		newImage->setPinPoint(false);

	//	newImage->newShoebox(8, 10, 16);
		
		newImage->setMaxResolution(40);
		newImage->setInitialStep(10);
		newImage->setSearchSize(20);
		CCP4SPG *spg = ccp4spg_load_by_standard_num(146);
		newImage->setSpaceGroup(spg);

		newImage->addMask(0, 840, 1765, 920);
		newImage->addMask(446, 1046, 1324, 1765);
		newImage->addMask(820, 450, 839, 473);

        for (int j = 0; j < newImage->indexerCount(); j++)
        {
            newImage->getIndexer(j)->findSpots();
        }
		newImage->dropImage();
	}
}

void indexThread(int offset, vector<MtzPtr>**mtzList,
		vector<Image *> *images)
{
    int maxThreads = FileParser::getMaxThreads();
    
	for (int i = offset; i < (*images).size(); i += maxThreads)
	{
        Image *newImage = (*images)[i];
        for (int j = 0; j < newImage->indexerCount(); j++)
        {
            IndexerPtr indexer = newImage->getIndexer(j);
            
            indexer->setTestSpotSize(0.002);
            indexer->setTestBandwidth(0.2);
            
            indexer->matchMatrixToSpots();
            indexer->setTestBandwidth(0.2);
            
            indexer->refineRoundBeamAxis();
            
            vector<MtzPtr> managers = newImage->currentMtzs();
            
            for (int i = 0; i < managers.size(); i++)
            {
                MtzPtr manager = managers[i];
                manager->description();
                manager->writeToFile(manager->getFilename());
                manager->writeToDat();
                (*mtzList)->push_back(manager);
            }
            
            newImage->dropImage();
        }
	}

}

void index(char **argv, int argc)
{
	boost::thread_group threads;

	vector<vector<Image *> *> allImages;
	vector<Image *> images;
	vector<vector<MtzPtr> *> mtzs;
    
    int maxThreads = FileParser::getMaxThreads();

	for (int i = 0; i < maxThreads; i++)
	{
		vector<Image *> *imageList = new vector<Image *>();
		allImages.push_back(imageList);

		vector<MtzPtr> *mtzList = new vector<MtzPtr>();
		mtzs.push_back(mtzList);
	}

	for (int i = 0; i < maxThreads; i++)
	{
		boost::thread *thr = new boost::thread(findSpots, argv, argc, i,
				&allImages[i]);
		threads.add_thread(thr);
	}

	threads.join_all();

	for (int i = 0; i < maxThreads; i++)
	{
		for (int j = 0; j < allImages[i]->size(); j++)
		{
			images.push_back((*allImages[i])[j]);
		}
	}

    std::ostringstream logged;
    logged << "Image count: " << images.size() << std::endl;
    Logger::mainLogger->addStream(&logged);

	for (int i = 0; i < images.size(); i++)
	{
		std::string name = images[i]->getFilename();
		int lastindex = (int)name.find_last_of(".");
		std::string rootName = name.substr(0, lastindex);
		std::string datName = rootName + ".dat";

        for (int j = 0; j < images[i]->indexerCount(); j++)
        {
            images[i]->getIndexer(j)->writeDatFromSpots(datName);
        }
	}

	boost::thread_group secondThread;

	for (int i = 0; i < maxThreads; i++)
	{
		boost::thread *thr = new boost::thread(indexThread, i, &mtzs[i],
				&images);
		secondThread.add_thread(thr);
	}

	secondThread.join_all();

	std::ofstream newMats;
	newMats.open("new_orientations.dat");

	for (int i = 0; i < mtzs.size(); i++)
	{
		for (int j = 0; j < mtzs[i]->size(); j++)
		{
			MtzManager *manager = &(*(*mtzs[i])[j]);

			// write out matrices etc.
            std::string imgFilename = manager->filenameRoot();
			newMats << imgFilename << " ";

			std::string description = manager->getMatrix()->description();
			newMats << description;
		}
	}

    Logger::mainLogger->addString("Written to new_orientations.dat");

	newMats.close();
}
