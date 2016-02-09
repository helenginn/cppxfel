//
//  CSV.cpp
//  cppxfel
//
//  Created by Helen Ginn on 09/02/2016.
//  Copyright (c) 2016 Division of Structural Biology Oxford. All rights reserved.
//

#include "CSV.h"
#include <fstream>

CSV::CSV(int count, ...)
{
    va_list arguments;
    va_start(arguments, count);
    
    for (int i = 0; i < count; i++)
    {
        std::string header = std::string(va_arg(arguments, char *));
        headers.push_back(header);
    }
    
    va_end(arguments);
}

void CSV::addEntry(int dummy, ...)
{
    va_list arguments;
    va_start(arguments, dummy);
    Entry newEntry;
    
    for (int i = 0; i < headers.size(); i++)
    {
        double value = va_arg(arguments, double);
        newEntry.push_back(value);
    }
    
    entries.push_back(newEntry);
}

void CSV::writeToFile(std::string filename)
{
    std::ofstream csv;
    csv.open(filename.c_str());
    
    for (int i = 0; i < headers.size(); i++)
    {
        csv << headers[i] << ",";
    }
    
    csv << std::endl;
    
    for (int i = 0; i < entries.size(); i++)
    {
        Entry anEntry = entries[i];
        
        for (int j = 0; j < anEntry.size(); j++)
        {
            csv << anEntry[j] << ",";
        }
        
        csv << std::endl;
    }
    
    csv.close();
}

double CSV::valueForEntry(std::string header, int entry)
{
    for (int i = 0; i < headerCount(); i++)
    {
        if (headers[i] == header)
        {
            return entries[entry][i];
        }
    }
}