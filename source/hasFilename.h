//
//  hasFilename.h
//  cppxfel
//
//  Created by Helen Ginn on 22/03/2017.
//  Copyright (c) 2017 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef cppxfel_hasFilename_h
#define cppxfel_hasFilename_h

#include <string>

class hasFilename
{
private:
    std::string filename;

public:
    const std::string& getFilename() const
    {
        return filename;
    }

    void setFilename(const std::string& filename)
    {
        this->filename = filename;
    }

    std::string getBasename()
    {
        int fullStopIndex = (int)filename.rfind(".");
        if (fullStopIndex == std::string::npos)
            return filename;

        std::string basename = filename.substr(0, fullStopIndex);

        return basename;
    }
};

#endif
