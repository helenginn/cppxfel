//
//  FreeMillerLibrary.h
//  cppxfel
//
//  Created by Helen Ginn on 17/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__FreeMillerLibrary__
#define __cppxfel__FreeMillerLibrary__

#include <stdio.h>

#include <vector>
#include <string>
#include "LoggableObject.h"

class FreeMillerLibrary
{
private:
    FreeMillerLibrary(std::string filename);
    FreeMillerLibrary(std::string filename, double maxResolution);
    std::vector<vec> freeIndices;
    
    static FreeMillerLibrary library;
    
    int freeIndexCount()
    {
        return (int)freeIndices.size();
    }
    
public:
    static void setup();
    
    void printSummary();
    
    bool isCertainMillerFree(Miller *miller);
    
    static bool isMillerFree(Miller *miller)
    {
        return library.isCertainMillerFree(miller);
    }
    
    static bool active()
    {
        return (library.freeIndexCount() > 0);
    }
    
    FreeMillerLibrary() {};
};

#endif /* defined(__cppxfel__FreeMillerLibrary__) */
