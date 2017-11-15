//
//  FileReader.h
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __GameDriver__FileReader__
#define __GameDriver__FileReader__

#include <iostream>
#include <string>
#include "parameters.h"
#include <vector>

namespace FileReader
{
    std::string get_file_contents(const char *filename);

    vector<std::string> split(const std::string s, const std::string &delim);
    vector<std::string> &split(const std::string &s, char delim, vector<std::string> &elems);
    vector<std::string> split(const std::string &s, char delim);
    bool exists(const std::string& name);

    std::string addOutputDirectory(std::string filename);
    int splitAtIndices(const std::string &s, vector<int> &positions, vector<std::string> &elems);
};

#endif /* defined(__GameDriver__FileReader__) */
