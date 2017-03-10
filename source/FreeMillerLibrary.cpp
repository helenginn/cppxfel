//
//  FreeMillerLibrary.cpp
//  cppxfel
//
//  Created by Helen Ginn on 17/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "FreeMillerLibrary.h"
#include "FileReader.h"
#include "FileParser.h"
#include "csymlib.h"
#include "parameters.h"
#include "UnitCellLattice.h"
#include <sstream>
#include "CSV.h"
#include "Miller.h"

FreeMillerLibrary FreeMillerLibrary::library;

FreeMillerLibrary::FreeMillerLibrary(std::string filename)
{
    
}

FreeMillerLibrary::FreeMillerLibrary(std::string filename, double maxResolution)
{
    
}

void FreeMillerLibrary::setup()
{
    
}

void FreeMillerLibrary::printSummary()
{
    
}

bool FreeMillerLibrary::isCertainMillerFree(Miller *miller)
{
    return false;
}