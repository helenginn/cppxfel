//
//  GeometryParser.h
//  cppxfel
//
//  Created by Helen Ginn on 27/12/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__GeometryParser__
#define __cppxfel__GeometryParser__

#include <stdio.h>
#include "parameters.h"
#include "LoggableObject.h"

class GeometryParser : public LoggableObject
{
private:
    std::string filename;
    GeometryFormat format;
    
    void parseCrystFELLines(std::vector<std::string> lines);
    void parseCppxfelLines(std::vector<std::string> lines);
    
    DetectorPtr makeDetector(int min_fs, int min_ss, int max_fs, int max_ss,
                             double ss_x, double ss_y, double ss_z,
                             double fs_x, double fs_y, double fs_z,
                             double midpoint_x, double midpoint_y, double midpoint_z,
                             double alpha, double beta, double gamma);
public:
    GeometryParser(std::string aFilename, GeometryFormat aFormat = GeometryFormatCrystFEL);
    
    void parse();
    void writeToFile(std::string newName);
};

#endif /* defined(__cppxfel__GeometryParser__) */
