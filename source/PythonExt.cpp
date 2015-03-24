//
//  PythonExt.cpp
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "PythonExt.h"
#include "InputFileParser.h"
#include "Logger.h"

void setupCppxfel()
{
    time_t startcputime;
    time(&startcputime);
    
    Logger::mainLogger = LoggerPtr(new Logger());
    boost::thread thr = boost::thread(Logger::awaitPrintingWrapper, Logger::mainLogger);
   
    std::cout << "Welcome to Helen's XFEL tasks" << std::endl;
}

void runScriptFromPython(std::string scriptName)
{
    setupCppxfel();
    
    InputFileParser *parser = new InputFileParser(scriptName);
    parser->parse();
    
    delete parser;
}