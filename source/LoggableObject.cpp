//
//  LoggableObject.cpp
//  cppxfel
//
//  Created by Helen Ginn on 12/11/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "LoggableObject.h"

void LoggableObject::sendLog(LogLevel priority, bool shouldExit)
{
    Logger::mainLogger->addStream(&logged, priority, shouldExit);
    logged.str("");
    logged.clear();
}

void LoggableObject::sendLogAndExit()
{
    staticLogAndExit(logged);
}

void LoggableObject::staticLogAndExit(std::ostringstream &otherLog, std::string header)
{
    std::ostringstream errorStart;
    errorStart << "**************** " << header << " *****************" << std::endl;
    Logger::mainLogger->addStream(&errorStart, LogLevelNormal, false);

    otherLog << "**************** " << header << " *****************" << std::endl;
    Logger::mainLogger->addStream(&otherLog, LogLevelNormal, true);

    boost::this_thread::sleep(boost::posix_time::seconds(1));
}
