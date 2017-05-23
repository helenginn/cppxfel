//
//  IndexManager.cpp
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "IndexManager.h"
#include "FileParser.h"
#include "Matrix.h"
#include "MtzManager.h"
#include "Logger.h"
#include <algorithm>
#include <iomanip>
#include "parameters.h"
#include <fstream>
#include "SpotVector.h"
#include "IndexingSolution.h"
#include "FreeLattice.h"
#include "CSV.h"
#include "Reflection.h"
#include "Miller.h"
#include "misc.h"
#include "Detector.h"
#include <set>

int IndexManager::_cycleNum = 0;


std::string targetString(PseudoScoreType type)
{
	switch (type) {
		case PseudoScoreTypeAllInterPanel:
			return "all_interpanel";
			break;
		case PseudoScoreTypeInterPanel:
			return "interpanel";
			break;
		case PseudoScoreTypeIntraPanel:
			return "intrapanel";
			break;

		default:
			return "";
			break;
	}

}

IndexManager::IndexManager(std::vector<ImagePtr> newImages)
{
	writtenPDB = false;
    images = newImages;
    scoreType = PseudoScoreTypeIntraPanel;
    proportionDistance = FileParser::getKey("DISTANCE_VS_ANGLE_FRACTION", 1.0);
    _canLockVectors = false;
    _axisWeighting = PseudoScoreWeightingAxisNone;
    _maxFrequency = -1;

	interPanelDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);
	intraPanelDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);

    MatrixPtr rotationMat = MatrixPtr(new Matrix());
    
    maxDistance = 0;
    lattice = UnitCellLattice::getMainLattice();
 
    sendLog(LogLevelDebug);
}


double IndexManager::pseudoScore(void *object)
{
	return pseudoAngleScore(object);
}

void IndexManager::plotGoodVectors()
{
    if (!goodVectors.size())
    {
        return;
    }
    
    ImagePtr special = Detector::getSpecialImage();
    std::string filename = getActiveDetector()->getTag() + "_" + targetString(getPseudoScoreType()) + "_vec.png";
    
    special->plotVectorsOnPNG(goodVectors, filename);
}

CSVPtr IndexManager::pseudoAnglePDB(bool writePost)
{
	double angleDistance = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);

	if (angleDistance <= 0)
	{
		return CSVPtr();
	}

	CSVPtr angleCSV = CSVPtr(new CSV(3, "angle", "aLength", "bLength"));
	goodVectorPairs.clear();

	CSVPtr predictedAngles = getLattice()->getAngleCSV();

	SpotVectorPair nullPair;
	goodVectorPairs.push_back(nullPair);
	DetectorPtr myDet = getActiveDetector();

	logged << "Preparing vector pairs for " << myDet->getTag() << std::endl;
	sendLog();

	for (int i = 0; i < images.size() && angleCSV->entryCount() < 100000; i++)
	{
		for (int j = 0; j < images[i]->spotVectorCount(); j++)
		{
			SpotVectorPtr vec1 = images[i]->spotVector(j);
			vec1->calculateDistance();

			if (!vec1->originalDistanceLessThan(angleDistance))
			{
				continue;
			}

			if (scoreType != PseudoScoreTypeBeamCentre && !vec1->isOnlyFromDetector(myDet))
			{
				continue;
			}

			for (int k = 0; k < j; k++)
			{
				SpotVectorPtr vec2 = images[i]->spotVector(k);

				if (checkVectors(vec1, vec2) != scoreType)
				{
					continue;
				}

				vec2->calculateDistance();

				SpotVectorPtr vec3 = vec1->completeTriangleWith(vec2);
				if (vec3->originalDistanceLessThan(angleDistance))
				{
					continue;
				}

				SpotVectorPair pair;
				pair.vecs[0] = vec1;
				pair.vecs[1] = vec2;
				pair.vecs[2] = vec3;

				double angle = vec1->angleWithVector(vec2);
				angle *= 180 / M_PI;
				angle = (angle > 90) ? 180 - angle : angle;

				double dist1 = vec1->distance();
				double dist2 = vec2->distance();

				SpotVectorPtr minVec = (dist1 > dist2) ? vec2 : vec1;
				SpotVectorPtr maxVec = (dist1 <= dist2) ? vec2 : vec1;

				angleCSV->addEntry(3, angle, dist1, dist2);

				goodVectorPairs.push_back(pair);
				goodVectors.push_back(maxVec);
				goodVectors.push_back(minVec);
				goodVectors.push_back(vec3);

				goodSpots.insert(vec1->getFirstSpot());
				goodSpots.insert(vec1->getSecondSpot());
				goodSpots.insert(vec2->getFirstSpot());
				goodSpots.insert(vec2->getSecondSpot());
			}
		}
	}

	if (getActiveDetector() == Detector::getMaster())
	{
		return angleCSV;
	}

	std::string prePost = writePost ? "_post_" : "_pre_";
	std::string filename = "angles_" + getActiveDetector()->getTag() +
	"_" + targetString(scoreType) + prePost + i_to_str(getCycleNum()) + ".pdb";
	angleCSV->plotPDB(filename, "aLength", "bLength", "angle");

	writtenPDB = true;

	return angleCSV;
}

PseudoScoreType IndexManager::checkVectors(SpotVectorPtr vec1, SpotVectorPtr vec2)
{
    DetectorPtr activeDetector = getActiveDetector();
    bool isInterPanel = scoreType == PseudoScoreTypeAllInterPanel || scoreType == PseudoScoreTypeInterPanel;

	if (scoreType == PseudoScoreTypeBeamCentre)
	{
		if (vec1->usesBeamCentre() && vec2->usesBeamCentre())
		{
			return PseudoScoreTypeBeamCentre;
		}

		return PseudoScoreTypeInvalid;
	}

	if ((!vec1->originalDistanceLessThan(intraPanelDistance) || !vec2->originalDistanceLessThan(intraPanelDistance)))
	{
		return PseudoScoreTypeInvalid;
	}

	// intra panel dealings

	if (!isInterPanel && !vec1->hasCommonSpotWithVector(vec2))
	{
		return PseudoScoreTypeInvalid;
	}

	if (!isInterPanel && vec1->isIntraPanelVector() && vec2->isIntraPanelVector())
	{
		if (vec1->isOnlyFromDetector(activeDetector) && vec2->isOnlyFromDetector(activeDetector))
		{
			return PseudoScoreTypeIntraPanel;
		}

		return PseudoScoreTypeInvalid;
	}

	if (!isInterPanel)
	{
		return PseudoScoreTypeInvalid;
	}

	// now definitely interpanel (specific) or all interpanel

	if (scoreType == PseudoScoreTypeAllInterPanel)
	{
		if (vec1->isIntraPanelVector() || vec2->isIntraPanelVector())
		{
			return PseudoScoreTypeInvalid;
		}

		if (vec1->hasCommonSpotWithVector(vec2))
		{
			return PseudoScoreTypeAllInterPanel;
		}

		return PseudoScoreTypeInvalid;
	}

	// now definitely specific interpanel (not all)

	if (vec1->isIntraPanelVector() && vec2->isIntraPanelVector())
	{
		return PseudoScoreTypeInvalid;
	}

	if (!vec1->hasCommonSpotWithVector(vec2))
	{
		return PseudoScoreTypeInvalid;
	}

	if ((vec1->spansChildrenOfDetector(activeDetector) && (vec2->spansChildrenOfDetector(activeDetector))))
	{
		return PseudoScoreTypeInterPanel;
	}

	if ((vec1->isOnlyFromDetector(activeDetector) && (vec2->isOnlyFromDetector(activeDetector))
		 && ((vec1->spansChildrenOfDetector(activeDetector) || (vec1->spansChildrenOfDetector(activeDetector))))
		 && vec1->hasCommonSpotWithVector(vec2)))
	{
		return PseudoScoreTypeInterPanel;
	}

    return PseudoScoreTypeInvalid;
}

double IndexManager::pseudoAngleScore(void *object)
{
	double maxAngleDist = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);

	if (maxAngleDist <= 0)
	{
		return 0;
	}

	IndexManager *me = static_cast<IndexManager *>(object);

	if (me->goodVectorPairs.size() == 0)
	{
		me->pseudoAnglePDB();
	}

	double totalScore = 0;
	double sumDistances = 0;

	for (std::set<SpotPtr>::iterator it = me->goodSpots.begin(); it != me->goodSpots.end(); it++)
	{
		SpotPtr spot = *it;
		spot->storeEstimatedVector();
	}

	for (int i = 0; i < me->goodVectors.size(); i++)
	{
		me->goodVectors[i]->quickDistance();
	}

	for (int i = 1; i < me->goodVectorPairs.size(); i++)
	{
		SpotVectorPair pair = me->goodVectorPairs[i];

		double score = 1;

		for (int i = 0; i < 3; i++)
		{
			int j = (i + 1) % 3;

			double aLength = pair.vecs[i]->distance();
			double bLength = pair.vecs[j]->distance();
			double angle = pair.vecs[i]->angleWithVector(pair.vecs[j]);

			angle *= 180 / M_PI;
			angle = (angle > 90) ? 180 - angle : angle;
			sumDistances += aLength;
			sumDistances += bLength;

			score *= me->getLattice()->weightForPair(aLength, bLength, angle);
			totalScore += score;

		}

	}

	totalScore /= me->goodVectorPairs.size();

	return -totalScore;
}

// MARK: Actual indexing

void IndexManager::indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset)
{
	std::ostringstream logged;

	while (true)
	{
		ImagePtr image = indexer->getNextImage();

		if (!image)
		{
			logged << "Finishing thread " << offset << std::endl;
			Logger::mainLogger->addStream(&logged); logged.str("");
			return;
		}

		logged << "Starting image " << image->getFilename() << " on thread " << offset << std::endl;
		Logger::mainLogger->addStream(&logged); logged.str("");

		image->findIndexingSolutions();

		std::vector<MtzPtr> mtzs = image->currentMtzs();

		image->dropImage();

		mtzSubset->reserve(mtzSubset->size() + mtzs.size());
		mtzSubset->insert(mtzSubset->begin(), mtzs.begin(), mtzs.end());
	}
}

ImagePtr IndexManager::getNextImage()
{
    std::lock_guard<std::mutex> lock(indexMutex);

    nextImage++;

    if (nextImage >= images.size())
    {
        return ImagePtr();
    }
    else
    {
        return images[nextImage];
    }
}

void IndexManager::index()
{
    int maxThreads = FileParser::getMaxThreads();
    IndexingSolution::setupStandardVectors();
    
    boost::thread_group threads;
    vector<vector<MtzPtr> > managerSubsets;
    managerSubsets.resize(maxThreads);
    nextImage = -1;
    
    for (int num = 0; num < 1; num++)
    {
        time_t startcputime;
        time(&startcputime);

        for (int i = 0; i < maxThreads; i++)
        {
            boost::thread *thr = new boost::thread(indexThread, this, &managerSubsets[i], i);
            threads.add_thread(thr);
        }

        threads.join_all();
    
        time_t endcputime;
        time(&endcputime);
        
        nextImage = -1;
    }
    
    
    int total = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        total += managerSubsets[i].size();
    }
    
    mtzs.reserve(total);
    int lastPos = 0;
    
    for (int i = 0; i < maxThreads; i++)
    {
        mtzs.insert(mtzs.begin() + lastPos,
                           managerSubsets[i].begin(), managerSubsets[i].end());
        lastPos += managerSubsets[i].size();
    }
}

// MARK: Normal powder pattern

PowderHistogram IndexManager::generatePowderHistogram(int intraPanel, int perfectPadding)
{
    PowderHistogram frequencies;
    double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);
    
    double maxCatDistance = FileParser::getKey("MAX_RECIPROCAL_DISTANCE", 0.15);
    int maxCategory = maxCatDistance / step;
    
    for (int i = 0; i < maxCategory; i++)
    {
        frequencies[i] = std::make_pair(0, 0);
    }
    
    for (int i = 0; i < lattice->standardVectorCount(); i++)
    {
        double distance = lattice->standardVector(i)->distance();
        int categoryNum = distance / step;
        
        if (categoryNum > maxCategory)
            continue;
        
        for (int i = categoryNum - perfectPadding; i <= categoryNum + perfectPadding; i++)
        {
            if (i < 0 || i > maxCategory)
                continue;
            
            frequencies[i].second++;
        }
    }
    
    for (int i = 0; i < images.size(); i++)
    {
        for (int j = 0; j < images[i]->spotVectorCount(); j++)
        {
            SpotVectorPtr spotVec = images[i]->spotVector(j);

            double distance = spotVec->distance();
            int categoryNum = distance / step;
            
            if (categoryNum > maxCategory || categoryNum < 0)
                continue;

            bool isIntraPanel = spotVec->isIntraPanelVector();
            
            if (intraPanel == -1)
            {
                frequencies[categoryNum].first++;
            }
            else if (intraPanel == 0 && !isIntraPanel)
            {
                frequencies[categoryNum].first++;
            }
            else if (intraPanel == 1 && isIntraPanel)
            {
                frequencies[categoryNum].first++;
            }
        }
    }
    
    return frequencies;
}


void IndexManager::powderPattern(std::string csvName, bool force)
{
	bool alwaysFilterSpots = FileParser::getKey("ALWAYS_FILTER_SPOTS", false);

	std::ostringstream pdbLog;


	for (int i = 0; i < images.size(); i++)
	{
		if (force)
		{
			images[i]->compileDistancesFromSpots(maxDistance, smallestDistance, alwaysFilterSpots);
		}

		for (int j = 0; j < images[i]->spotVectorCount(); j++)
		{
			SpotVectorPtr spotVec = images[i]->spotVector(j);
			spotVec->calculateDistance();
		}
	}


	PowderHistogram allFrequencies = generatePowderHistogram();
	PowderHistogram intraFrequencies = generatePowderHistogram(1);
	PowderHistogram interFrequencies = generatePowderHistogram(0);

	if (images.size() == 0)
	{
		logged << "No images specified." << std::endl;
		sendLog();
		return;
	}

	logged << "******* DISTANCE FREQUENCY *******" << std::endl;

	double step = FileParser::getKey("POWDER_PATTERN_STEP", 0.00005);

	CSV powder(5, "Distance", "Frequency", "Perfect frequency", "Intra-panel", "Inter-panel");

	for (PowderHistogram::iterator it = allFrequencies.begin(); it != allFrequencies.end(); it++)
	{
		double distance = it->first * step;
		double freq = it->second.first;
		double perfect = it->second.second;
		double intrapanel = intraFrequencies[it->first].first;
		double interpanel = interFrequencies[it->first].first;

		powder.addEntry(0, distance, freq, perfect, intrapanel, interpanel);
	}

	if (angleCSV)
	{
		angleCSV->writeToFile("angle_" + csvName);
	}
	powder.writeToFile(csvName);

	std::map<std::string, std::string> plotMap;
	plotMap["xHeader0"] = "Distance";
	plotMap["xTitle0"] = "Reciprocal distance between vectors (Ang)";
	plotMap["yHeader0"] = "Intra-panel";
	plotMap["yTitle0"] = "Frequency of distance";
	plotMap["style0"] = "line";
	plotMap["round0"] = "true";
	plotMap["xHeader1"] = "Distance";
	plotMap["xTitle1"] = "Reciprocal distance between vectors (Ang)";
	plotMap["yHeader1"] = "Inter-panel";
	plotMap["style1"] = "line";
	plotMap["colour1"] = "blue";
	plotMap["filename"] = csvName;
	plotMap["xHeader2"] = "Distance";
	plotMap["xTitle2"] = "Reciprocal distance between vectors (Ang)";
	plotMap["yHeader2"] = "Perfect frequency";
	plotMap["style2"] = "line";
	plotMap["colour2"] = "blue";

	powder.plotPNG(plotMap);
	logged << "Written to " << csvName << std::endl;
	sendLog();

	if (force)
	{
		powder.plotColumns(0, 1);
	}

	lattice->powderPattern(false, "unitCellLatticePowder.csv");
}
