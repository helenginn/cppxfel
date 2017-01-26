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

DetectorPtr GeometryParser::makeDetector(int min_fs, int min_ss, int max_fs, int max_ss,
                                         double ss_x, double ss_y, double ss_z,
                                         double fs_x, double fs_y, double fs_z,
                                         double midpoint_x, double midpoint_y, double midpoint_z,
                                         double alpha, double beta, double gamma)
{
    Coord topLeft = std::make_pair(min_fs, min_ss);
    Coord bottomRight = std::make_pair(max_fs, max_ss);
    
    vec slowAxis = new_vector(ss_x, ss_y, ss_z);
    vec fastAxis = new_vector(fs_x, fs_y, fs_z);
    vec middle = new_vector(midpoint_x, midpoint_y, midpoint_z);
    
    DetectorPtr segment = DetectorPtr(new Detector(DetectorPtr(), topLeft, bottomRight,
                                                   slowAxis, fastAxis, middle, true));
    
    segment->prepareRotationAngles(alpha, beta, gamma);
    
    return segment;
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
                    DetectorPtr segment = makeDetector(min_fs, min_ss, max_fs, max_ss,
                                                       ss_x, ss_y, ss_z, fs_x, fs_y, fs_z,
                                                       midpoint_x, midpoint_y, midpoint_z,
                                                       alpha, beta, gamma);
                    segment->setTag(panelStack[panelStack.size() - 1]);
                    
                    if (!Detector::getMaster())
                    {
                        Detector::setMaster(segment);
                    }
                    else
                    {
                        DetectorPtr lastDetector = detectorStack[detectorStack.size() - 1];
                        lastDetector->addChild(segment);
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
    }
    
    Detector::getMaster()->applyRotations();
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
        Logger::mainLogger->addStream(&logged, LogLevelNormal, true);
    }
    
    Coord topLeft = std::make_pair(min_fs, min_ss);
    Coord bottomRight = std::make_pair(max_fs, max_ss);
    vec fastAxis = new_vector(fs_x, fs_y, fs_z);
    vec slowAxis = new_vector(ss_x, ss_y, ss_z);
    vec arrangedTopLeft = new_vector(corner_x, corner_y, 0);
    
    DetectorPtr segment = DetectorPtr(new Detector(parent, topLeft, bottomRight,
                                                   slowAxis, fastAxis, arrangedTopLeft));
    
    segment->setRotationAngles(0, 0, 0);
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
    
    for (int i = 0; i < lines.size(); i++)
    {
        std::string line = lines[i];
        
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
    
    logged << "Making the assumption that this detector is a " << (isSacla ? "MPCCD" : "CSPAD") << " detector." << std::endl;
    sendLog();
    
    if (isSacla)
    {
        Detector::setDetectorType(DetectorTypeMPCCD);
        
        for (int i = 0; i < panelMaps.size(); i++)
        {
            DetectorPtr segment = makeDetectorFromPanelMap(panelMaps[i], master);
        }
        
    }
    else if (!isSacla) // for now, CSPAD
    {
        Detector::setDetectorType(DetectorTypeCSPAD);

        std::vector<DetectorPtr> quarters;
        /* Top right */
        DetectorPtr q0 = DetectorPtr(new Detector(master, new_vector(+441, +441, 0), "q0"));
        master->addChild(q0);
        quarters.push_back(q0);
        
        /* Bottom left */
        DetectorPtr q1 = DetectorPtr(new Detector(master, new_vector(-441, +441, 0), "q1"));
        master->addChild(q1);
        quarters.push_back(q1);
        
        /* Top left */
        DetectorPtr q2 = DetectorPtr(new Detector(master, new_vector(-441, -441, 0), "q2"));
        master->addChild(q2);
        quarters.push_back(q2);
        
        /* Bottom right */
        DetectorPtr q3 = DetectorPtr(new Detector(master, new_vector(+441, -441, 0), "q3"));
        master->addChild(q3);
        quarters.push_back(q3);

        for (int i = 0; i < panelMaps.size(); i++)
        {
            PanelMap map = panelMaps[i];
            std::string name = map["name"];
            if (name.length() < 4)
            {
                logged << "Huh?" << std::endl;
                sendLog();
            }
            
            int qNum = atoi(name.substr(1, 1).c_str());
            
            DetectorPtr myParent;
            makeDetectorFromPanelMap(panelMaps[i], quarters[qNum]);
        }
    }

    master->setRotationAngles(0, 0, 0);
}

void GeometryParser::parse()
{
    if (!FileReader::exists(filename))
    {
        logged << "File " << filename << " does not exist." << std::endl;
        Logger::mainLogger->addStream(&logged, LogLevelNormal, true);
    }
    
    std::string contents = FileReader::get_file_contents(filename.c_str());
    std::vector<std::string> lines = FileReader::split(contents, '\n');
    
    if (format == GeometryFormatCrystFEL)
    {
        parseCrystFELLines(lines);
    }
    else if (format == GeometryFormatCppxfel)
    {
        parseCppxfelLines(lines);
    }
    
    Detector::fullDescription();
}

void GeometryParser::writeToFile(std::string newName)
{
    std::string geometry = Detector::getMaster()->writeGeometryFile();
    
    std::ofstream file;
    file.open(newName.c_str());
    
    file << "# Geometry file (cppxfel format v2)" << std::endl << std::endl;

    file << "version 2" << std::endl;
    file << geometry;
    
    file.close();
    
}