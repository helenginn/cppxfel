//
//  PythonExt.h
//  cppxfel
//
//  Created by Helen Ginn on 24/03/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__PythonExt__
#define __cppxfel__PythonExt__

#include <string>
#include <vector>

void runScriptFromPython(std::string scriptName);
void runCommandLineArgs(int argc, std::vector<std::string> argv[]);

#endif /* defined(__cppxfel__PythonExt__) */
