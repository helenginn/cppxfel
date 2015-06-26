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

class MtzRefiner;

class InputFileParser : public FileParser
{
private:
    boost::shared_ptr<MtzRefiner> refiner;
public:
    void refine(int maxCycles);
    vector<MtzPtr> mtzs();
    
	virtual void parse();

	InputFileParser(std::string filename);
	virtual ~InputFileParser();
    
    boost::shared_ptr<MtzRefiner> getRefiner()
    {
        return refiner;
    }
};

#endif /* INPUTFILEPARSER_H_ */
