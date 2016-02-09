//
//  CSV.h
//  cppxfel
//
//  Created by Helen Ginn on 09/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__CSV__
#define __cppxfel__CSV__

#include <stdio.h>
#include <vector>
#include <string>
#include "parameters.h"
#include <cstdarg>

typedef std::vector<double> Entry;

class CSV
{
private:
    std::vector<std::string> headers;
    std::vector<Entry> entries;
    
public:
    CSV(int count, ...);
    
    void addEntry(int dummy, ...);
    void writeToFile(std::string filename);
};

#endif /* defined(__cppxfel__CSV__) */
