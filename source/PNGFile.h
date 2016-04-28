//
//  Image.h
//  CellCounter
//
//  Created by Helen Ginn on 07/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __CellCounter__Image__
#define __CellCounter__Image__

#include <vector>
#include <map>
#include <stdio.h>
#include <png.h>
#include <string>

typedef enum
{
    FlagNone, FlagThreshold
} Flag;

class PNGFile
{
private:
    png_bytep data;
    int height;
    int centreX;
    int centreY;
    std::string filename;
    std::string rootName;
    int bytesPerPixel;
    int pixelsPerRow;
    int writeImage(std::string filename, int width, int height, std::string title);
    void read_png_file(const char *file_name);
    void pixelAt(int x, int y, png_byte **bytes);
    png_byte valueAt(int x, int y);
    bool pixelReachesThreshold(int x, int y);
    void setPixelColour(int x, int y, png_byte red, png_byte green, png_byte blue);
    void setPixelForChannel(int x, int y, int channel, png_byte value);
    void readImage();
    void moveCoordRelative(int *x, int *y);
public:
    void dropImage();
    void preProcess();
    void writeImageOutput();
    void process();
    void setPixelColourRelative(int x, int y, png_byte red, png_byte green, png_byte blue);
    void drawCircleAroundPixel(int x, int y, float radius, png_byte red, png_byte green, png_byte blue);
    
    PNGFile(std::string filename, int width = 2400, int height = 2400);
    
    void setCentre(int newX, int newY)
    {
        centreX = newX;
        centreY = newY;
    }
};



#endif /* defined(__CellCounter__Image__) */
