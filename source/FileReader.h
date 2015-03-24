//
//  FileReader.h
//  GameDriver
//
//  Created by Helen Ginn on 21/05/2014.
//  Copyright (c) 2014 Helen Ginn. All rights reserved.
//

#ifndef __GameDriver__FileReader__
#define __GameDriver__FileReader__

#include <iostream>
#include <string>
#include <vector>

namespace FileReader
{
    std::string get_file_contents(const char *filename);
    
    std::vector<std::string> split(const std::string s, const std::string &delim);
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    std::vector<std::string> split(const std::string &s, char delim);
    bool exists(const std::string& name);
    
    int splitAtIndices(const std::string &s, std::vector<int> &positions, std::vector<std::string> &elems);
};

#endif /* defined(__GameDriver__FileReader__) */
