//
//  Shoebox.h
//  cppxfel
//
//  Created by Helen Ginn on 26/01/2015.
//  Copyright (c) 2015 Helen Ginn. All rights reserved.
//

#ifndef __cppxfel__Shoebox__
#define __cppxfel__Shoebox__

#include <stdio.h>
#include "parameters.h"
#include <vector>
#include <boost/weak_ptr.hpp>
#include "LoggableObject.h"

typedef vector<vector<double> > Box;

class Shoebox : LoggableObject
{
private:
    boost::weak_ptr<Miller> miller;
    Box shoebox;

    double mmPerPixel;
    double detectorDistance;

    double pixelLeak;
    double bandwidthMultiplier;

    bool even;
    void putBoxOnPaddedBackground(Box &smallBox, Box &newShoebox);
    void putBoxOnBackground(Box &smallBox, Box &newShoebox);
    void centreShoebox(Box &smallBox);
    void chopBox(Box &smallBox);
    void compressBigShoebox(double width, double height, double angle, Box &smallBox);
public:
    Shoebox(MillerPtr parent);
    ~Shoebox();

    void printShoebox();
    void simpleShoebox(int foregroundLength, int neitherLength, int backgroundLength, bool shoeboxEven);
    void complexShoebox(double wavelength, double bandwidth, double radius);
    void clearShoebox();
    void sideLengths(int *slowSide, int *fastSide);
    void centre(int *centreX, int *centreY);

    vector<double>& operator[](std::size_t index)
    {
        return shoebox[index];
    }

    vector<double>& row(int index)
    {
        return shoebox[index];
    }

    int height()
    {
        return (int)shoebox.size();
    }

    MillerPtr getMiller()
    {
        return miller.lock();
    }

    void setMiller(MillerPtr miller)
    {
        this->miller = miller;
    }

    bool isEven()
    {
        return even;
    }
};

#endif /* defined(__cppxfel__Shoebox__) */
