//
//  FileReader.cpp
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

#include "FileReader.h"

#include "misc.h"
#include <fstream>
#include <sstream>
#include <cerrno>
#include <sys/stat.h>
#include "FileParser.h"
#include "Logger.h"
#include <sys/types.h>
#include <dirent.h>

vector<std::string> &FileReader::split(const std::string &s, char delim, vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<std::string> FileReader::split(const std::string s, const std::string &delim)
{
    vector<std::string> elems;
    std::string rest = s;

    int count = 0;
    bool finished = false;

    if (s.length() == 0)
        return elems;

    while (!finished)
    {
        count++;
        size_t index = rest.substr(1, rest.length() - 1).find(delim);

        if (index == std::string::npos)
        {
            index = rest.length() - 1;
            finished = true;
        }

        std::string cutout = rest.substr(1, index + 1);
        elems.push_back(cutout);

        rest = rest.substr(index + 1, s.length() - index - 1);

        if (index == 0)
            break;
    }

    return elems;
}

vector<std::string> FileReader::split(const std::string &s, char delim) {
    vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

int FileReader::splitAtIndices(const std::string &s, vector<int> &positions, vector<std::string> &elems) {

    for (int i = 0; i < positions.size() - 1; i++)
    {
        int start = positions[i];
        int length = positions[i + 1] - start;

        if (s.size() < positions[i + 1])
            return 0;

        std::string segment = s.substr(start, length);
        elems.push_back(segment);
    }

    return 1;
}

bool FileReader::exists(const std::string& name)
{
        struct stat buffer;
        return (stat(name.c_str(), &buffer) == 0);
}

std::string FileReader::get_file_contents(const char *filename)
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);

    if (in)
    {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize((unsigned long)in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return(contents);
    }

    std::string errString = "Could not get file contents for file " + std::string(filename);
    Logger::mainLogger->addString(errString);
    Logger::mainLogger->addString(strerror(errno));

    throw(errno);
}


std::string FileReader::addOutputDirectory(std::string filename)
{
    std::string directory = FileParser::getKey("OUTPUT_DIRECTORY", std::string(""));
    bool outputIndividualCycles = FileParser::getKey("OUTPUT_INDIVIDUAL_CYCLES", false);

    if (!directory.length() && !outputIndividualCycles)
    {
        return filename;
    }
    else
    {
    DIR *dir = opendir(directory.c_str());

        if (dir)
        {
            closedir(dir);
        }
        else if (ENOENT == errno)
        {
            mkdir(directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }
    }

    if (outputIndividualCycles)
    {
        int cycleNum = FileParser::getKey("CYCLE_NUMBER", 0);
        std::string combinedDir;

        if (directory.length())
        {
            combinedDir = directory + "/cycle_" + i_to_str(cycleNum);
        }
        else
        {
            combinedDir = "cycle_" + i_to_str(cycleNum);
        }

        DIR *dir2 = opendir(combinedDir.c_str());

        if (dir2)
        {
            closedir(dir2);
        }
        else if (ENOENT == errno)
        {
            mkdir(combinedDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }

        return combinedDir + "/" + filename;
    }
    else
    {
        return directory + "/" + filename;
    }
}
