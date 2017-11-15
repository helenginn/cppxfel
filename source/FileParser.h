/*
 * FileParser.h
 */
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

#ifndef FILEPARSER_H_
#define FILEPARSER_H_

#include <string>
#include "parameters.h"
#include <iostream>
#include "LoggableObject.h"
class MtzRefiner;
typedef std::map<std::string, int> CodeMap;

typedef std::map<std::string, std::vector<std::string> > CategoryMap;
typedef std::map<std::string, CategoryMap > CategoryTree;


class FileParser: public LoggableObject
{
protected:
    static ParserMap parserMap;
        static ParametersMap parameters;
    static CategoryTree categoryTree;
    static std::map<std::string, std::string> deprecatedList;
    static std::map<std::string, std::string> helpMap;
    static std::map<std::string, CodeMap> codeMaps;

        static void simpleFloat(ParametersMap *map, std::string command,
                        std::string rest);
        static void simpleBool(ParametersMap *map, std::string command,
                        std::string rest);
        static void simpleString(ParametersMap *map, std::string command,
                        std::string rest);
        static void simpleInt(ParametersMap *map, std::string command,
                        std::string rest);
        static void doubleVector(ParametersMap *map, std::string command,
                        std::string rest);
        static void intVector(ParametersMap *map, std::string command,
                        std::string rest);
    static void stringVector(ParametersMap *map, std::string command,
                             std::string rest);

    std::vector<std::string> extras;
    static char splitCharMajor;
    static char splitCharMinor;
    static int threadsFound;
        std::string filename;
    void generateCodeList();
    void generateHelpList();
    void generateCategories();
    void generateFunctionList();
    void generateDeprecatedList();
        ParserFunction splitLine(std::string line, std::string &command, std::string &rest);
        bool checkSpaces(std::string line);
    static std::ostringstream log;
public:

    static int getMaxThreads();

        template<typename Value>
    static void setKey(std::string key, Value newValue)
    {
        ParameterVariant container = newValue;
        parameters[key] = container;
    }

    template<typename Value>
        static Value getKey(std::string key, Value defaultValue)
        {
        if (hasKey(key))
                {
                        ParameterVariant container = parameters[key];
                        ParameterVariant container2 = defaultValue;

                        if (container.which() != container2.which())
                        {
                                log << "Incorrect type for " << key << std::endl;
                        }

            try
            {
                Value value = boost::get<Value>(parameters[key]);
                return value;
            }
            catch (boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::bad_get> > &ex)
            {
                std::cout << "Failed to get value from array" << std::endl;
                std::cout << "Key: " << key << std::endl;

                throw ex;
            }
                }
                else {
                        return defaultValue;
                }
        }
        static bool hasKey(std::string key);
    FileParser(void);
    FileParser(std::string filename, std::vector<std::string> someExtras = std::vector<std::string>());
    static void printCommandInfo(std::string command);
    static void printAllCommands();
        virtual ~FileParser();
    virtual void parseFromPython() { parse(true); };
        virtual void parse(bool fromPython) {};
};

#endif /* FILEPARSER_H_ */
