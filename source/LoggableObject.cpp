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