//
//  Logger.h
//  cppxfel
//
//  Created by Helen Ginn on 18/01/2015.
//  Copyright (c) 2015 Helen Ginn. All rights reserved.
//

#ifndef __cppxfel__Logger__
#define __cppxfel__Logger__

#include <boost/thread/thread.hpp>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <map>
#include <vector>
#include "parameters.h"
#include <mutex>
#include <condition_variable>

typedef std::pair<StreamPtr, LogLevel> LogAndLevel;

typedef std::map<boost::thread::id, vector<LogAndLevel> > StringMap;


class Logger
{
private:
    StringMap stringsToOutput;
    std::mutex mtx;
    std::mutex writing;
    std::condition_variable printBlock;
    static bool ready;
    LogLevel printedLogLevel;
    bool tryLock(std::mutex &lock, int maxTries = 50);
    static bool isReady();
    static bool shouldExit;

public:
    Logger();
    ~Logger();

    static LoggerPtr mainLogger;

    void addString(std::string string, LogLevel level = LogLevelNormal);
    void changePriorityLevel(LogLevel newLevel);
    void awaitPrinting();
    static void awaitPrintingWrapper(LoggerPtr logger);
    void addStream(std::ostringstream *stream, LogLevel level = LogLevelNormal, bool shouldExit = false);

    static void setShouldExit()
    {
        shouldExit = true;
    }

    static LogLevel getPriorityLevel()
    {
        return mainLogger->printedLogLevel;
    }

    static void log(std::ostringstream &stream)
    {
        mainLogger->addStream(&stream);
        stream.str("");
        stream.clear();
    }
};

#endif /* defined(__cppxfel__Logger__) */
