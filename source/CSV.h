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
#include "LoggableObject.h"
#include "StatisticsManager.h"

typedef enum
{
    PlotVerticalLine = '|',
    PlotHorizontalLine = '_',
    PlotHorizontalTickMark = '-',
    PlotVerticalTickMark = '\'',
    PlotBlank = ' ',
    
} PlotChar;

typedef std::vector<double> Entry;
typedef std::map<int, char> Row;
typedef std::map<int, Row > Plot;

class CSV : public LoggableObject
{
private:
    std::vector<std::string> headers;
    std::vector<Entry> entries;
    double minX, minY, maxX, maxY;
    bool didSetMinMaxXY;
    void minMaxCol(int col, double *min, double *max);
    std::string mapToAscii(Plot plot);
    void writeStringToPlot(std::string text, Plot *plot, int x, int y);
public:
    CSV()
    {
        
    }

    CSV(int count, ...)
    {
        didSetMinMaxXY = false;
        
        va_list arguments;
        va_start(arguments, count);
        
        for (int i = 0; i < count; i++)
        {
            std::string header = std::string(va_arg(arguments, char *));
            addHeader(header);
        }
        
        va_end(arguments);
    }
    
    void setupHistogram(double start, double end, double interval, std::string catHeader, int count, ...)
    {
        va_list arguments;
        va_start(arguments, count);
        
        headers.push_back(catHeader);
        
        for (int i = 0; i < count; i++)
        {
            std::string header = std::string(va_arg(arguments, char *));
            headers.push_back(header);
        }
        
        va_end(arguments);
        
        if (interval < 0)
        {
            std::vector<double> bins;
            StatisticsManager::generateResolutionBins(start, end, -interval, &bins);
            
            for (int i = 0; i < bins.size(); i++)
            {
                addPartialEntry(1, bins[i]);
            }
        }
        else if (interval > 0)
        {
            double value = start;
            while (value <= end)
            {
                addPartialEntry(1, value);
                value += interval;
            }
        }
    }
    
    void addOneToFrequency(double category, std::string whichHeader, double weight = 1, std::string categoryHeader = "");
    int findHeader(std::string whichHeader);

    ~CSV();
    
    void addPartialEntry(int dummy, ...);
    void addEntry(int dummy, ...);
    void writeToFile(std::string filename);
    double valueForEntry(std::string header, int entry);
    double valueForHistogramEntry(std::string whichHeader, double value, std::string categoryHeader = "");
    double valueForHistogramEntry(int whichHeader, double value, int categoryHeader = 0);
    void histogram(std::map<double, int> histogram);
    std::string plotColumns(int col1, int col2);
    void resetColumn(std::string header, double value = 0);
    
    void convolutedPeaks(std::string category, std::string origHeader, std::string destHeader, double stdev);
    void addConvolutedPeak(std::string header, double mean, double stdev, double weight, std::string category = "");
    void setValueForEntry(int entry, std::string header, double value);
    void addEntry(std::vector<double> entry)
    {
        entries.push_back(entry);
    }
    
    int entryCount()
    {
        return (int)entries.size();
    }
    
    int headerCount()
    {
        return (int)headers.size();
    }
    
    void addHeader(std::string header)
    {
        headers.push_back(header);
    }
    
    void setMinMaxXY(double _minX, double _minY, double _maxX, double _maxY)
    {
        minX = _minX;
        minY = _minY;
        maxX = _maxX;
        maxY = _maxY;
        didSetMinMaxXY = true;
    }
};

#endif /* defined(__cppxfel__CSV__) */
