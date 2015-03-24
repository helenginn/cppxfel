/*
 * InputFileParser.h
 *
 *  Created on: 13 Jan 2015
 *      Author: helenginn
 */

#ifndef INPUTFILEPARSER_H_
#define INPUTFILEPARSER_H_

#include "FileParser.h"

class InputFileParser : public FileParser
{
public:
	virtual void parse();

	InputFileParser(std::string filename);
	virtual ~InputFileParser();
};

#endif /* INPUTFILEPARSER_H_ */
