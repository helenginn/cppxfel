//
//  GeometryParser.cpp
//  cppxfel
//
//  Created by Helen Ginn on 27/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "FileParser.h"
#include "FileReader.h"
#include "GeometryParser.h"
#include "misc.h"
#include "Detector.h"
#include "Vector.h"
#include <fstream>


GeometryParser::GeometryParser(std::string aFilename, GeometryFormat aFormat)
{
    filename = aFilename;
    
    format = aFormat;
}

DetectorPtr GeometryParser::makeDetector(DetectorPtr parent, int min_fs, int min_ss, int max_fs, int max_ss,
                                         double ss_x, double ss_y, double ss_z,
                                         double fs_x, double fs_y, double fs_z,
                                         double midpoint_x, double midpoint_y, double midpoint_z,
                                         double alpha, double beta, double gamma, bool ghost, bool refinable,
										 double gain)
{
    Coord topLeft = std::make_pair(min_fs, min_ss);
    Coord bottomRight = std::make_pair(max_fs, max_ss);
    
    vec slowAxis = new_vector(ss_x, ss_y, ss_z);
    vec fastAxis = new_vector(fs_x, fs_y, fs_z);
    vec middle = new_vector(midpoint_x, midpoint_y, midpoint_z);
    
    DetectorPtr segment = DetectorPtr(new Detector(parent));
    
    if (parent)
    {
        parent->addChild(segment);
    }
    segment->initialise(topLeft, bottomRight, slowAxis, fastAxis, middle, true, ghost);
    
    segment->prepareRotationAngles(alpha, beta, gamma);
	segment->setRefinable(refinable);
	segment->setGain(gain);
    
    return segment;
}

void GeometryParser::parsePanelListLines(std::vector<std::string> lines)
{
    logged << "Warning: 'panels list' definition of prior days of cppxfel is woefully incomplete." << std::endl;
    logged << "If you can parse in a CrystFEL or cppxfel (version 2) geometry file instead, this" << std::endl;
    logged << "is highly recommended." << std::endl;
    
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    std::vector<double> beamCentre = FileParser::getKey("BEAM_CENTRE", std::vector<double>(2, 0));
    
    DetectorPtr master = DetectorPtr(new Detector(detectorDistance, 0, 0));
    Detector::setMaster(master);
    
    for (int i = 0; i < lines.size(); i++)
    {
        std::string line = lines[i];
        
        double x1 = 0;
        double x2 = 0;
        double y1 = 0;
        double y2 = 0;
        double o1 = 0;
        double o2 = 0;
        double sw = 0;
        double gain = 1;
        
        trim(line);
        
        if (line[0] == '#')
        {
            continue;
        }

        std::vector<std::string> components = FileReader::split(line, ' ');
        
        if (!components.size())
            continue;
        
        if (components[0] == "PANEL")
        {
            std::string panelNum = "panel" + i_to_str(i);
            Coord arrangedTopLeft, arrangedBottomRight;
            
            
            if (components.size() >= 5)
            {
                x1 = atof(components[1].c_str());
                y1 = atof(components[2].c_str());
                x2 = atof(components[3].c_str());
                y2 = atof(components[4].c_str());
                
                arrangedTopLeft = std::make_pair(x1, y1);
                arrangedBottomRight = std::make_pair(x2, y2);
            }
            
            if (components.size() >= 7)
            {
                o1 = atof(components[5].c_str());
                o2 = atof(components[6].c_str());
            }
            
            if (components.size() >= 10)
            {
                sw = atof(components[9].c_str());
            }
            
            if (components.size() >= 11)
            {
                gain = atof(components[10].c_str());
            }
            
            DetectorPtr segment = DetectorPtr(new Detector(master, arrangedTopLeft, arrangedBottomRight, sw, o1, o2, gain));
            segment->setTag(panelNum);
            master->addChild(segment);
        }
    }
    
    Detector::setArrangedMidPointX(&*master, beamCentre[0]);
    Detector::setArrangedMidPointY(&*master, beamCentre[1]);

    for (int i = 0; i < master->childrenCount(); i++)
    {
        master->getChild(i)->removeMidPointRelativeToParent();
    }
    
    Detector::setArrangedMidPointX(&*master, 0);
    Detector::setArrangedMidPointY(&*master, 0);

}

void GeometryParser::parseCppxfelLines(std::vector<std::string> lines)
{
    std::vector<std::string> panelStack;
    std::vector<DetectorPtr> detectorStack;
    int bracesIn = 0;
    
    bool nextLineShouldBeBrace = false;
    
    int    min_fs = 0;
    int    min_ss = 0;
    int    max_fs = 0;
    int    max_ss = 0;
    double fs_x = 0;
    double fs_y = 0;
    double fs_z = 0;
    double ss_x = 0;
    double ss_y = 0;
    double ss_z = 0;
    double midpoint_x = 0;
    double midpoint_y = 0;
    double midpoint_z = 0;
    double alpha = 0;
    double beta = 0;
    double gamma = 0;
    bool ghost = false;
    bool refinable = false;
	double gain = 1;
    
    for (int i = 0; i < lines.size(); i++)
    {
        std::string line = lines[i];
        
        trim(line);
        
        if (line[0] == '#')
        {
            continue;
        }
        
        std::vector<std::string> components = FileReader::split(line, ' ');
        
        if (components.size() < 1)
            continue;
        
        if (nextLineShouldBeBrace && components[0] != "{")
        {
            logged << "Warning: expecting { after start of panel " << panelStack.at(panelStack.size() - 1) << "." << std::endl;
            sendLog();
        }
        else if (components[0] == "{")
        {
            bracesIn++;
        }
        
        nextLineShouldBeBrace = false;

        if (components[0] == "panel" || components[0] == "}")
        {
            if (panelStack.size())
            {
                if (components[0] == "}" && detectorStack.size() == panelStack.size())
                {
                    bracesIn--;
                    detectorStack.pop_back();
                    panelStack.pop_back();
                    continue;
                }
                
                if (!(components[0] == "panel" && detectorStack.size() == panelStack.size()))
                {
                    DetectorPtr parent = DetectorPtr();
                    
                    if (Detector::getMaster())
                    {
                        parent = detectorStack[detectorStack.size() - 1];
                        //parent->setRefinable(true);
                    }
                    
                    std::string tag = panelStack[panelStack.size() - 1];
                    DetectorPtr segment = makeDetector(parent, min_fs, min_ss, max_fs, max_ss,
                                                       ss_x, ss_y, ss_z, fs_x, fs_y, fs_z,
                                                       midpoint_x, midpoint_y, midpoint_z,
                                                       alpha, beta, gamma, ghost, refinable, gain);
                    segment->setTag(tag);
                    
                    if (!Detector::getMaster())
                    {
                        Detector::setMaster(segment);
                    }
                    
                    detectorStack.push_back(segment);
                    logged << "Added detector " << segment->getTag() << " at hierarchy level " << detectorStack.size() << std::endl;
                    sendLog();
                }
            }

            if (components[0] == "panel")
            {
                nextLineShouldBeBrace = true;
                
                if (components.size() < 2)
                {
                    continue;
                }
                
                panelStack.push_back(components[1]);
            }
            else
            {
                bracesIn--;
                panelStack.pop_back();
                detectorStack.pop_back();
            }
        }
        
        if (components.size() < 3)
            continue;
        
        if (components[0] == "min_fs")
            min_fs = atoi(components[2].c_str());
        if (components[0] == "min_ss")
            min_ss = atoi(components[2].c_str());

        if (components[0] == "max_fs")
            max_fs = atoi(components[2].c_str());
        if (components[0] == "max_ss")
            max_ss = atoi(components[2].c_str());

        if (components[0] == "fs_x")
            fs_x = atof(components[2].c_str());
        if (components[0] == "fs_y")
            fs_y = atof(components[2].c_str());
        if (components[0] == "fs_z")
            fs_z = atof(components[2].c_str());

        if (components[0] == "ss_x")
            ss_x = atof(components[2].c_str());
        if (components[0] == "ss_y")
            ss_y = atof(components[2].c_str());
        if (components[0] == "ss_z")
            ss_z = atof(components[2].c_str());

        if (components[0] == "midpoint_x")
            midpoint_x = atof(components[2].c_str());
        if (components[0] == "midpoint_y")
            midpoint_y = atof(components[2].c_str());
        if (components[0] == "midpoint_z")
            midpoint_z = atof(components[2].c_str());

        if (components[0] == "alpha")
            alpha = atof(components[2].c_str());
        if (components[0] == "beta")
            beta = atof(components[2].c_str());
        if (components[0] == "gamma")
            gamma = atof(components[2].c_str());
        
        if (components[0] == "ghost")
            ghost = (strcmp(components[2].c_str(), "true") == 0);
        if (components[0] == "refine")
            refinable = (strcmp(components[2].c_str(), "true") == 0);
		if (components[0] == "gain")
			gain = atof(components[2].c_str());
    }

    if (!Detector::getMaster())
    {
        throw std::exception();
    }
    
    Detector::getMaster()->updateCurrentRotation();
}

DetectorPtr GeometryParser::makeDetectorFromPanelMap(std::vector<PanelMap> panelMaps, std::string tag, DetectorPtr parent)
{
    for (int i = 0; i < panelMaps.size(); i++)
    {
        if (panelMaps[i]["name"] == tag)
        {
            return makeDetectorFromPanelMap(panelMaps[i], parent);
        }
    }
    
    logged << "Couldn't find tag " << tag << std::endl;
    sendLog();
    
    return DetectorPtr();
}

bool detectorFurtherThan(DetectorPtr a, DetectorPtr b)
{
	vec aVec = a->getAverageMidpoint();
	vec bVec = b->getAverageMidpoint();

	return (length_of_vector(aVec) > length_of_vector(bVec));
}

DetectorPtr GeometryParser::makeDetectorFromPanelMap(PanelMap panelMap, DetectorPtr parent)
{
    int    min_fs;
    int    min_ss;
    int    max_fs;
    int    max_ss;
    double fs_x;
    double fs_y;
    double fs_z;
    double ss_x;
    double ss_y;
    double ss_z;
    double corner_x;
    double corner_y;
    std::string name;
    
    try
    {
        min_fs = atoi(panelMap["min_fs"].c_str());
        min_ss = atoi(panelMap["min_ss"].c_str());
        max_fs = atoi(panelMap["max_fs"].c_str());
        max_ss = atoi(panelMap["max_ss"].c_str());
        fs_x = atof(panelMap["fs_x"].c_str());
        fs_y = atof(panelMap["fs_y"].c_str());
        fs_z = atof(panelMap["fs_z"].c_str());
        ss_x = atof(panelMap["ss_x"].c_str());
        ss_y = atof(panelMap["ss_y"].c_str());
        ss_z = atof(panelMap["ss_z"].c_str());
        corner_x = atof(panelMap["corner_x"].c_str());
        corner_y = atof(panelMap["corner_y"].c_str());
        name = panelMap["name"];
    }
    catch (std::exception e)
    {
        logged << "Missing a vital parameter in panel" << std::endl;
        sendLogAndExit();
    }
    
    Coord topLeft = std::make_pair(min_fs, min_ss);
    Coord bottomRight = std::make_pair(max_fs, max_ss);
    vec fastAxis = new_vector(fs_x, fs_y, fs_z);
    vec slowAxis = new_vector(ss_x, ss_y, ss_z);
    vec arrangedTopLeft = new_vector(corner_x, corner_y, 0);
    
    DetectorPtr segment = DetectorPtr(new Detector(parent));
    segment->initialise(topLeft, bottomRight, slowAxis, fastAxis, arrangedTopLeft);
    
    segment->setTag(name);
    parent->addChild(segment);

    return segment;
}

void GeometryParser::parseCrystFELLines(std::vector<std::string> lines)
{
    std::vector<PanelMap> panelMaps;
    
    std::string lastPanel = "";
    std::map<std::string, std::string> panelMap;
    bool isSacla = false;

	if (FileParser::getKey("DETECTOR_DISTANCE", 0.) <= 0)
	{
		logged << "Detector distance not set or set to zero/negative." << std::endl;
		logged << "You must specify a detector distance for CrystFEL geometry format!" << std::endl;
		sendLogAndExit();
	}

    for (int i = 0; i < lines.size(); i++)
    {
        std::string line = lines[i];
        
        if (line[0] == ';')
        {
            continue;
        }
        
        std::vector<std::string> forwardSlashes = FileReader::split(line, '/');
        
        if (forwardSlashes.size() < 2)
        {
            /* this is not a panel definition */
            continue;
        }
        
        unsigned long nameLength = forwardSlashes[0].length();
        
        if (!(nameLength == 2 || (nameLength >= 4 && nameLength <= 5)))
        {
            continue;
        }
        
        if (forwardSlashes[0] != lastPanel)
        {
            if (lastPanel.length())
            {
                panelMaps.push_back(panelMap);
            }
            
            panelMap = std::map<std::string, std::string>();
            panelMap["name"] = forwardSlashes[0];
            lastPanel = forwardSlashes[0];
            
            if (!isSacla)
                isSacla = (forwardSlashes[0].length() == 2 && forwardSlashes[0][0] == 'q');
        }
        
        
        /* Hacky as f.*/
        
        std::vector<std::string> equals = FileReader::split(forwardSlashes[1], '=');
        
        if (equals.size() < 2)
        {
            /* this is a weird one, ignore... */
            continue;
        }
        
        trim(equals[0]);
        trim(equals[1]);
        
        if (equals[0] == "fs" || equals[0] == "ss")
        {
            // need to split into _x + _y
            std::string *str = &equals[1];
            
            str->erase(std::remove(str->begin(), str->end(), 'x'), str->end());
            str->erase(std::remove(str->begin(), str->end(), 'y'), str->end());
            str->erase(std::remove(str->begin(), str->end(), 'z'), str->end());
            
            std::vector<std::string> spaces = FileReader::split(equals[1], ' ');
            
            panelMap[equals[0] + "_x"] = spaces[0];
            
            if (spaces.size() > 1)
            {
                panelMap[equals[0] + "_y"] = spaces[1];
            }
            else
            {
                panelMap[equals[0] + "_y"] = "0";
            }
            
            if (spaces.size() > 2)
            {
                panelMap[equals[0] + "_z"] = spaces[2];
            }
            else
            {
                panelMap[equals[0] + "_z"] = "0";
            }
        }
        else
        {
            panelMap[equals[0]] = equals[1];
        }
    }
    
    panelMaps.push_back(panelMap);
    
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    std::vector<double> beamCentre = FileParser::getKey("BEAM_CENTRE", std::vector<double>(2, 0));
    
    DetectorPtr master = DetectorPtr(new Detector(detectorDistance, beamCentre[0], beamCentre[1]));
    Detector::setMaster(master);
    master->setRefinable(refineQuadrantGeometry);
    
    logged << "Making the assumption that this detector is a " << (isSacla ? "MPCCD" : "CSPAD") << " detector." << std::endl;
    sendLog();
    
    if (isSacla)
    {
        Detector::setDetectorType(DetectorTypeMPCCD);
        
        DetectorPtr q1256 = DetectorPtr(new Detector(master, new_vector(0, 0, 0), "q1256"));
        master->addChild(q1256);
        q1256->setRefinable(refineLocalGeometry);
        
        DetectorPtr q12 = DetectorPtr(new Detector(q1256, new_vector(0, 0, 0), "q12"));
        q1256->addChild(q12);
        q12->setRefinable(refineLocalGeometry);

        DetectorPtr q1 = makeDetectorFromPanelMap(panelMaps, "q1", q12);
        q1->setRefinable(refineLocalGeometry);
        
        DetectorPtr q2 = makeDetectorFromPanelMap(panelMaps, "q2", q12);
        q2->setRefinable(refineLocalGeometry);
        
        DetectorPtr q56 = DetectorPtr(new Detector(q1256, new_vector(0, 0, 0), "q56"));
        q1256->addChild(q56);
        q56->setRefinable(refineLocalGeometry);
        
        DetectorPtr q5 = makeDetectorFromPanelMap(panelMaps, "q5", q56);
        q5->setRefinable(refineLocalGeometry);
        
        DetectorPtr q6 = makeDetectorFromPanelMap(panelMaps, "q6", q56);
        q6->setRefinable(refineLocalGeometry);

        DetectorPtr q3478 = DetectorPtr(new Detector(master, new_vector(0, 0, 0), "q3478"));
        master->addChild(q3478);
        q3478->setRefinable(refineLocalGeometry);
        
        DetectorPtr q34 = DetectorPtr(new Detector(q3478, new_vector(0, 0, 0), "q34"));
        q3478->addChild(q34);
        q34->setRefinable(refineLocalGeometry);
        
        DetectorPtr q3 = makeDetectorFromPanelMap(panelMaps, "q3", q34);
        q3->setRefinable(refineLocalGeometry);
        
        DetectorPtr q4 = makeDetectorFromPanelMap(panelMaps, "q4", q34);
        q4->setRefinable(refineLocalGeometry);
        
        DetectorPtr q78 = DetectorPtr(new Detector(q3478, new_vector(0, 0, 0), "q78"));
        q3478->addChild(q78);
        q34->setRefinable(refineLocalGeometry);
        
        DetectorPtr q7 = makeDetectorFromPanelMap(panelMaps, "q7", q78);
        q7->setRefinable(refineLocalGeometry);
        
        DetectorPtr q8 = makeDetectorFromPanelMap(panelMaps, "q8", q78);
        q8->setRefinable(refineLocalGeometry);

    }
	else if (!isSacla)
	{
		Detector::setDetectorType(DetectorTypeCSPAD);

		std::vector<DetectorPtr> quarters;
		DetectorPtr q0_1 = DetectorPtr(new Detector(master, new_vector(0, 0, 0), "q0_1"));
		q0_1->setRefinable(refineQuadrantGeometry);
		master->addChild(q0_1);

		DetectorPtr q2_3 = DetectorPtr(new Detector(master, new_vector(0, 0, 0), "q2_3"));
		q2_3->setRefinable(refineQuadrantGeometry);
		master->addChild(q2_3);

		// hard coded for now - convert to some kind of input

		std::vector<std::vector<DetectorPtr> > pairs;
		/* Top right */

		DetectorPtr q0 = DetectorPtr(new Detector(q0_1, new_vector(0, 0, 0), "q0"));
		q0_1->addChild(q0);
		q0->setRefinable(refineQuadrantGeometry);
		quarters.push_back(q0);

		DetectorPtr q1 = DetectorPtr(new Detector(q0_1, new_vector(0, 0, 0), "q1"));
		q0_1->addChild(q1);
		q1->setRefinable(refineQuadrantGeometry);
		quarters.push_back(q1);

		DetectorPtr q2 = DetectorPtr(new Detector(q2_3, new_vector(0, 0, 0), "q2"));
		q2_3->addChild(q2);
		q2->setRefinable(refineQuadrantGeometry);
		quarters.push_back(q2);

		DetectorPtr q3 = DetectorPtr(new Detector(q2_3, new_vector(0, 0, 0), "q3"));
		q2_3->addChild(q3);
		q3->setRefinable(refineQuadrantGeometry);
		quarters.push_back(q3);

		for (int i = 0; i < quarters.size(); i++)
		{
			pairs.push_back(std::vector<DetectorPtr>());
			DetectorPtr qX = quarters[i];
			std::string qXName = qX->getTag();

			DetectorPtr qXa0_8 = DetectorPtr(new Detector(qX, new_vector(0, 0, 0), qXName + "a0_8"));
			qXa0_8->setRefinable(refineLocalGeometry);
			qX->addChild(qXa0_8);

			DetectorPtr qXa0_3 = DetectorPtr(new Detector(qXa0_8, new_vector(0, 0, 0), qXName + "a0_3"));
			qXa0_3->setRefinable(refineLocalGeometry);
			qXa0_8->addChild(qXa0_3);

			DetectorPtr qXa0_1 = DetectorPtr(new Detector(qXa0_3, new_vector(0, 0, 0), qXName + "a0_1")); //
			qXa0_1->setRefinable(refineLocalGeometry);
			qXa0_3->addChild(qXa0_1);
			pairs[i].push_back(qXa0_1);

			DetectorPtr qXa2_3 = DetectorPtr(new Detector(qXa0_3, new_vector(0, 0, 0), qXName + "a2_3")); //
			qXa2_3->setRefinable(refineLocalGeometry);
			qXa0_3->addChild(qXa2_3);
			pairs[i].push_back(qXa2_3);

			DetectorPtr qXa4_7 = DetectorPtr(new Detector(qXa0_8, new_vector(0, 0, 0), qXName + "a4_7"));
			qXa4_7->setRefinable(refineLocalGeometry);
			qXa0_8->addChild(qXa4_7);

			DetectorPtr qXa4_5 = DetectorPtr(new Detector(qXa4_7, new_vector(0, 0, 0), qXName + "a4_5")); //
			qXa4_5->setRefinable(refineLocalGeometry);
			qXa4_7->addChild(qXa4_5);
			pairs[i].push_back(qXa4_5);

			DetectorPtr qXa6_7 = DetectorPtr(new Detector(qXa4_7, new_vector(0, 0, 0), qXName + "a6_7")); //
			qXa6_7->setRefinable(refineLocalGeometry);
			qXa4_7->addChild(qXa6_7);
			pairs[i].push_back(qXa6_7);

			DetectorPtr qXa8_15 = DetectorPtr(new Detector(qX, new_vector(0, 0, 0), qXName + "a8_15"));
			qXa8_15->setRefinable(refineLocalGeometry);
			qX->addChild(qXa8_15);

			DetectorPtr qXa8_11 = DetectorPtr(new Detector(qXa8_15, new_vector(0, 0, 0), qXName + "a8_11"));
			qXa8_11->setRefinable(refineLocalGeometry);
			qXa8_15->addChild(qXa8_11);

			DetectorPtr qXa8_9 = DetectorPtr(new Detector(qXa8_11, new_vector(0, 0, 0), qXName + "a8_9")); //
			qXa8_9->setRefinable(refineLocalGeometry);
			qXa8_11->addChild(qXa8_9);
			pairs[i].push_back(qXa8_9);

			DetectorPtr qXa10_11 = DetectorPtr(new Detector(qXa8_11, new_vector(0, 0, 0), qXName + "a10_11")); //
			qXa10_11->setRefinable(refineLocalGeometry);
			qXa8_11->addChild(qXa10_11);
			pairs[i].push_back(qXa10_11);

			DetectorPtr qXa12_15 = DetectorPtr(new Detector(qXa8_15, new_vector(0, 0, 0), qXName + "a12_15"));
			qXa12_15->setRefinable(refineLocalGeometry);
			qXa8_15->addChild(qXa12_15);

			DetectorPtr qXa12_13 = DetectorPtr(new Detector(qXa12_15, new_vector(0, 0, 0), qXName + "a12_13")); //
			qXa12_13->setRefinable(refineLocalGeometry);
			qXa12_15->addChild(qXa12_13);
			pairs[i].push_back(qXa12_13);

			DetectorPtr qXa14_15 = DetectorPtr(new Detector(qXa12_15, new_vector(0, 0, 0), qXName + "a14_15")); //
			qXa14_15->setRefinable(refineLocalGeometry);
			qXa12_15->addChild(qXa14_15);
			pairs[i].push_back(qXa14_15);

		}

		for (int i = 0; i < panelMaps.size(); i++)
		{
			PanelMap map = panelMaps[i];
			std::string name = map["name"];
			if (name.length() < 4)
			{
				logged << "I'm struggling..." << std::endl;
				sendLog();
			}

			int aLength = (int)name.length() - 3;
			int qNum = atoi(name.substr(1, 1).c_str());
			int aNum = atoi(name.substr(3, aLength).c_str());
			int pairNum = aNum / 2;

			DetectorPtr myParent;
			DetectorPtr segment = makeDetectorFromPanelMap(panelMaps[i], pairs[qNum][pairNum]);
			segment->setRefinable(false);
		}
		
		/*
		std::vector<DetectorPtr> pairs;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 16; j += 2)
			{
				std::string qaCombo;
				qaCombo = "q" + i_to_str(i) + "a" + i_to_str(j) + "_" + i_to_str(j + 1);
				DetectorPtr asicPair = DetectorPtr(new Detector(master, new_vector(0, 0, 0), qaCombo));
				asicPair->setRefinable(true);
				pairs.push_back(asicPair);
			}
		}

		for (int i = 0; i < panelMaps.size(); i++)
		{
			PanelMap map = panelMaps[i];
			std::string name = map["name"];
			if (name.length() < 4)
			{
				logged << "I'm struggling..." << std::endl;
				sendLog();
			}

			int aLength = (int)name.length() - 3;
			int qNum = atoi(name.substr(1, 1).c_str());
			int aNum = atoi(name.substr(3, aLength).c_str());
			int pairNum = qNum * 8 + aNum / 2;

			DetectorPtr myParent;
			DetectorPtr segment = makeDetectorFromPanelMap(panelMaps[i], pairs[pairNum]);
			segment->setRefinable(false);
		}

		std::sort(pairs.begin(), pairs.end(), detectorFurtherThan);
		DetectorPtr lastParent = Detector::getMaster();

		for (int i = 0; i < pairs.size(); i++)
		{
			std::string setName = "panel_set_" + i_to_str(i);
			DetectorPtr theRest = DetectorPtr(new Detector(DetectorPtr(), new_vector(0, 0, 0), setName));
			theRest->setRefinable(true);

			lastParent->addChild(pairs[i]);

			if (i < pairs.size() - 1)
			{
				lastParent->addChild(theRest);
			}

			lastParent = theRest;
		}*/
	}

    Detector::getMaster()->updateCurrentRotation();
    Detector::getMaster()->fixMidpoints();
}

void GeometryParser::parse()
{
    if (!FileReader::exists(filename))
    {
        logged << "File " << filename << " does not exist." << std::endl;
        sendLogAndExit();
    }
    
    refineGlobalGeometry = !FileParser::getKey("TRUST_GLOBAL_GEOMETRY", false);
    refineQuadrantGeometry = !FileParser::getKey("TRUST_QUADRANT_GEOMETRY", false);
    refineLocalGeometry = !FileParser::getKey("TRUST_LOCAL_GEOMETRY", false);
    
    if (refineQuadrantGeometry)
    {
        refineGlobalGeometry = true;
    }
    
    if (refineLocalGeometry)
    {
        refineQuadrantGeometry = true;
        refineGlobalGeometry = true;
    }
    
    std::string contents = FileReader::get_file_contents(filename.c_str());
    std::vector<std::string> lines = FileReader::split(contents, '\n');
    
    if (format == GeometryFormatCrystFEL)
    {
        try
        {
            parseCrystFELLines(lines);
        }
        catch (std::exception)
        {
            logged << "I had trouble trying to parse that geometry format and crashed." << std::endl;
            logged << "Sorry. Are you sure the GEOMETRY_FORMAT should be CrystFEL?" << std::endl;
            sendLogAndExit();
        }
    }
    
    if (format == GeometryFormatCppxfel)
    {
        try
        {
            parseCppxfelLines(lines);
        }
        catch (std::exception)
        {
            logged << "I had trouble trying to parse that geometry format and crashed." << std::endl;
            logged << "Sorry. Are you sure the GEOMETRY_FORMAT should be cppxfel?" << std::endl;
            sendLogAndExit();
        }
    }
    else if (format == GeometryFormatPanelList)
    {
        parsePanelListLines(lines);
    }
    
	Detector::getMaster()->postInit();
	filename = "converted.cgeom";
    writeToFile(filename, 0);

    Detector::fullDescription();
}

void GeometryParser::writeToFile(std::string newName, int fileCount)
{
    Detector::getMaster()->prepareInterNudges();
    std::string geometry = Detector::getMaster()->writeGeometryFile(fileCount);
    
    std::ofstream file;
    file.open(FileReader::addOutputDirectory(newName.c_str()));
    
    file << "# Geometry file (cppxfel format v2)" << std::endl << std::endl;

    file << "version 2" << std::endl;
    file << geometry;
    
    file.close();

	std::string crystFELName = newName + ".crystfel";
	std::string crystFELContents = Detector::getMaster()->writeCrystFELFile();

	std::ofstream crystFile;
	crystFile.open(FileReader::addOutputDirectory(crystFELName.c_str()));
	crystFile << crystFELContents;
	crystFile.close();

}