/*
 * InputFileParser.h
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#ifndef INPUTFILEPARSER_H_
#define INPUTFILEPARSER_H_

#include <boost/shared_ptr.hpp>
#include "FileParser.h"
//#include <boost/python.hpp>
#include "MtzRefiner.h"


class InputFileParser : public FileParser
{
private:
    boost::shared_ptr<MtzRefiner> refiner;
    int processOptions(std::vector<std::string> lines);
    bool dry;

public:
    void refine(int maxCycles);

    virtual void parseFromPython();
        virtual void parse(bool fromPython = false);

        InputFileParser(std::string filename, std::vector<std::string> extras = std::vector<std::string>());
        virtual ~InputFileParser();

    void setDry(bool newDry)
    {
        dry = newDry;
    }

    bool isDry()
    {
        return dry;
    }

    boost::shared_ptr<MtzRefiner> getRefiner()
    {
        return refiner;
    }

    void integrate();
};

#endif /* INPUTFILEPARSER_H_ */
