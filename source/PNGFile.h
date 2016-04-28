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

typedef png_byte DataPoint;

typedef std::pair<int, int> Coord;
typedef std::vector<Coord> CoordList;
typedef std::vector<CoordList> CoordLists;

typedef std::map<int, DataPoint> DataRow;
typedef std::map<int, DataRow> DataBlock;
typedef std::map<int, std::map<int, Flag> > FlagBlock;

class PNGFile
{
private:
    DataBlock data;
    FlagBlock flags;
    int channel;
    int height;
    int rejectCount;
    std::string filename;
    std::string rootName;
    bool brightness;
    int bytesPerPixel;
    int pixelsPerRow;
    int columnCount;
    int writeImage(std::string filename, int width, int height, std::string title);
    void read_png_file(const char *file_name);
    void pixelAt(int x, int y, png_byte **bytes);
    png_byte valueAt(int x, int y);
    Flag flagAt(int x, int y);
    void setFlag(int x, int y, Flag flag);
    bool pixelReachesThreshold(int x, int y);
    void setPixelColour(int x, int y, png_byte red, png_byte green, png_byte blue);
    void setPixelsForChannel(std::vector<Coord> coords, int channel, png_byte value);
    void setPixelForChannel(int x, int y, int channel, png_byte value);
    int totalPixelCount(std::vector<CoordList> &coordLists);
    void setPixelsColour(CoordList coords, png_byte red, png_byte green, png_byte blue);
    void setFlags(int x, int y, int add, Flag flag);
    void getImage();
public:
    
    void setChannel(int newChannel)
    {
        channel = newChannel;
    }
    
    void setRejectCount(int count)
    {
        rejectCount = count;
    }
    
    void setBrightness(bool bright)
    {
        brightness = bright;
    }
    
    void dropImage();
    void preProcess();
    void writeImageOutput();
    void process();
    PNGFile(std::string filename);
};



#endif /* defined(__CellCounter__Image__) */
