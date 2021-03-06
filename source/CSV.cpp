//
//  CSV.cpp
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "CSV.h"
#include <fstream>
#include <float.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include "FileReader.h"
#include "misc.h"
#include "PNGFile.h"

void CSV::histogram(std::map<double, int> histogram)
{
    for (std::map<double, int>::iterator it = histogram.begin(); it != histogram.end(); it++)
    {
        addEntry(0, it->first, (double)it->second);
    }
}

int CSV::findHeader(std::string whichHeader)
{
    int chosenHeader = -1;

    for (int i = 0; i < headers.size(); i++)
    {
        if (headers[i] == whichHeader)
        {
            chosenHeader = i;
            break;
        }
    }


    if (chosenHeader == -1)
    {
        logged << "Error: " << whichHeader << " does not exist in table." << std::endl;
        sendLogAndExit();
    }

    return chosenHeader;
}

double CSV::valueForHistogramEntry(std::string whichHeader, double value, std::string categoryHeader)
{
    if (entries.size() == 0)
        return 0;

    int categoryNum = 0;

    if (categoryHeader.length())
    {
        categoryNum = findHeader(categoryHeader);
    }

    int chosenHeader = findHeader(whichHeader);

    if (categoryNum < 0 || chosenHeader < 0)
    {
        logged << "Headers asked for in histogram entry (" << whichHeader << ", " << categoryHeader << ") are not acceptable." << std::endl;
        sendLog();
    }

    return valueForHistogramEntry(chosenHeader, value, categoryNum);
}

double CSV::valueForHistogramEntry(int chosenHeader, double value, int categoryNum)
{
    bool ascending = (entries[0][categoryNum] < entries[1][categoryNum]);

    std::vector<double>::iterator low;

    if (ascending)
    {
        low = std::lower_bound(histCategories.begin(), histCategories.end(), value, std::less<double>());
    }
    else
    {
        low = std::lower_bound(histCategories.begin(), histCategories.end(), value, std::greater<double>());
    }

    long num = low - histCategories.begin();

    if (num >= histCategories.size())
    {
        return 0;
    }

    return entries[num][chosenHeader];
}

void CSV::addOneToFrequency(double category, int column, double weight, int categoryNum)
{
    bool ascending = (entries[0][categoryNum] < entries[1][categoryNum]);

    for (int j = 0; j < entries.size() - 1; j++)
    {
        if ((ascending && category >= entries[j][categoryNum] && category < entries[j + 1][categoryNum])
            || (!ascending && category < entries[j][categoryNum] && category >= entries[j + 1][categoryNum]))
        {
            entries[j][column] += weight;
            break;
        }
    }
}


void CSV::addOneToFrequency(double category, std::string whichHeader, double weight, std::string categoryHeader)
{
    if (entries.size() == 0)
        return;

    int column = findHeader(whichHeader);

    int categoryNum = 0;

    if (categoryHeader.length())
    {
        categoryNum = findHeader(categoryHeader);
    }

    addOneToFrequency(category, column, weight, categoryNum);
}

void CSV::addPartialEntry(int dummy, ...)
{
    va_list arguments;
    va_start(arguments, dummy);
    Entry newEntry;

    for (int i = 0; i < dummy; i++)
    {
        double value = va_arg(arguments, double);
        newEntry.push_back(value);
    }

    int difference = (int)(headers.size() - newEntry.size());

    while (difference > 0)
    {
        newEntry.push_back(0);
        difference--;
    }

    entries.push_back(newEntry);
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
    std::string outputFile = FileReader::addOutputDirectory(filename);

    csv.open(outputFile.c_str());

    for (int i = 0; i < headers.size(); i++)
    {
        csv << headers[i] << ", ";
    }

    csv << std::endl;

    for (int i = 0; i < entries.size(); i++)
    {
        Entry anEntry = entries[i];

        for (int j = 0; j < anEntry.size(); j++)
        {
                        csv << std::setprecision(14) << anEntry[j] << ", ";
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

    return 0;
}

std::string CSV::mapToAscii(Plot plot)
{
    std::ostringstream stream;

    for (Plot::iterator it = plot.begin(); it != plot.end(); it++)
    {
        Row row = it->second;

        stream << "N: ";

        for (Row::iterator it2 = row.begin(); it2 != row.end(); it2++)
        {
            char charray[2];
            charray[0] = it2->second;
            charray[1] = '\0';

            stream << charray;
        }

        stream << std::endl;
    }

    return stream.str();
}

void CSV::minMaxCol(int col, double *min, double *max, bool round)
{
    *min = FLT_MAX;
    *max = -FLT_MAX;

    for (int i = 0; i < entries.size(); i++)
    {
        if (entries[i][col] > *max)
            *max = entries[i][col];

        if (entries[i][col] < *min)
            *min = entries[i][col];
    }

    if (round)
    {
        double logMax = log10(*max);

        double whole = floor(logMax);
        double remainder = fmod(logMax, 1);

        for (int i = 1; i <= 10; i++)
        {
            const double logi = log10(i);
            if (remainder <= logi)
            {
                whole += logi;
                break;
            }
        }

        *max = pow(10, whole);
    }
}

void CSV::plotPDB(std::string filename, std::string header1, std::string header2, std::string header3)
{
        int nums[3];
        double maxAngle = FileParser::getKey("MAXIMUM_ANGLE_DISTANCE", 0.0);

        nums[0] = findHeader(header1);
        nums[1] = findHeader(header2);
        nums[2] = findHeader(header3);

        int scales[3];

        scales[0] = 100 / maxAngle;
        scales[1] = 100 / maxAngle;
        scales[2] = 100 / 90;

        std::ofstream pdbLog;
        std::string fullName = FileReader::addOutputDirectory(filename);
        pdbLog.open(fullName.c_str());
        pdbLog << "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1  " << std::endl;

        for (int i = 0; i < entryCount(); i++)
        {
                pdbLog << "HETATM";
                pdbLog << std::fixed;
                pdbLog << std::setw(5) << i + 1 << "  O   HOH Z";
                pdbLog << std::setw(4) << 100;

                pdbLog << "    ";

                for (int j = 0; j < 3; j++)
                {
                        double valForEntry = entries[i][nums[j]] * scales[j];
                        pdbLog << std::setprecision(2) << std::setw(8)  << valForEntry;
                }

                pdbLog << "  1.00 30.00           O" << std::endl;
        }

        pdbLog.close();
}

void CSV::writeStringToPlot(std::string text, Plot *plot, int y, int x)
{
    for (int i = 0; i < text.length(); i++)
    {
        (*plot)[y][x + i] = text[i];
    }
}

typedef enum
{
    GraphStyleScatter,
    GraphStyleLine,
    GraphStyleHeatMap,
} GraphStyle;

void CSV::plotPNG(std::map<std::string, std::string> properties)
{
    int count = 0;

    double height = 900.;

    if (properties.count("height"))
    {
        height = atof(properties["height"].c_str());
    }

    double width = 1200.;

    if (properties.count("width"))
    {
        width = atof(properties["width"].c_str());
    }

    const int xIntervals = 5;
    const int yIntervals = 5;
    const double xAxis = 0.8;
    const double yAxis = 0.7;

    std::string filename = "graph_test.png";

    if (properties.count("filename"))
    {
        filename = properties["filename"] + ".png";
    }

    PNGFilePtr png = PNGFilePtr(new PNGFile(filename, width, height));
    png->setPlain();
    png->setCentre(width * (1 - xAxis) / 2, height * (yAxis + 1) / 2);
    png->drawLine(0, 0, width * xAxis, 0, 0, 0, 0, 0);
    png->drawLine(0, 0, 0, -height * yAxis, 0, 0, 0, 0);

    while (true)
    {
        std::string headerKey = "Header" + i_to_str(count);
        if (!properties.count("x" + headerKey) || !properties.count("y" + headerKey))
        {
            break;
        }

        std::string xHeader = properties["x" + headerKey];
        std::string yHeader = properties["y" + headerKey];

        int xCol = findHeader(xHeader);
        int yCol = findHeader(yHeader);

        if (xCol < 0 || yCol < 0)
        {
            logged << "Cannot draw graph. Headers incompatible: " << xHeader << ", " << yHeader << std::endl;
            sendLog();
            continue;
        }

        std::string minKey = "Min" + i_to_str(count);
        std::string maxKey = "Max" + i_to_str(count);
        std::string roundKey = "round" + i_to_str(count);
        double minX = -FLT_MAX;
        double minY = -FLT_MAX;
        double maxX = FLT_MAX;
        double maxY = FLT_MAX;
        bool round = (properties.count(roundKey) > 0);

        if (properties.count("x" + minKey)) minX = atof(properties["x" + minKey].c_str());
        if (properties.count("y" + minKey)) minY = atof(properties["y" + minKey].c_str());
        if (properties.count("x" + maxKey)) maxX = atof(properties["x" + maxKey].c_str());
        if (properties.count("y" + maxKey)) maxY = atof(properties["y" + maxKey].c_str());

        if (minX == -FLT_MAX || maxX == FLT_MAX) minMaxCol(xCol, &minX, &maxX, false);
        if (minY == -FLT_MAX || maxY == FLT_MAX) minMaxCol(yCol, &minY, &maxY, round);

        if (count == 0)
        {
            // draw X axis
            double aValue = minX;
            double aStep = (maxX - minX) / (double)xIntervals;
            double xPrec = 3;
            if (maxX > 10) xPrec = 0;

            for (int i = 0; i < xIntervals + 1; i++)
            {
                double proportion = (double)i / (double)xIntervals;
                double pos = width * xAxis * proportion;
                png->drawLine(pos, 0, pos, 5, 0, 0, 0, 0);

                png->drawText(f_to_str(aValue, xPrec), pos, 30);

                aValue += aStep;
            }

            std::string xTitle = xHeader;
            if (properties.count("xTitle" + i_to_str(count)))
            {
                xTitle = properties["xTitle" + i_to_str(count)];
            }

            png->drawText(xTitle, width * xAxis * 0.5, 85);

            aValue = minY;
            aStep = (maxY - minY) / (double)yIntervals;
            double yPrec = 2;
            if (maxY > 10) yPrec = 0;

            for (int i = 0; i < yIntervals + 1; i++)
            {
                double proportion = (double)i / (double)yIntervals;
                double pos = -height * yAxis * proportion;
                png->drawLine(0, pos, -5, pos, 0, 0, 0, 0);
                png->drawText(f_to_str(aValue, yPrec), -50, pos);

                aValue += aStep;
            }

            std::string yTitle = yHeader;
            if (properties.count("yTitle" + i_to_str(count)))
            {
                yTitle = properties["yTitle" + i_to_str(count)];
            }

            png->drawText(yTitle, 70, -height * (yAxis + (1 - yAxis) / 4));
        }

        GraphStyle style = GraphStyleScatter;
        std::string styleKey = "style" + i_to_str(count);
        if (properties.count(styleKey))
        {
            if (properties[styleKey] == "scatter")
            {
                style = GraphStyleScatter;
            }
            else if (properties[styleKey] == "line")
            {
                style = GraphStyleLine;
            }
            else if (properties[styleKey] == "line")
            {
                style = GraphStyleHeatMap;
            }
        }

        png_byte red = 0;
        png_byte blue = 0;
        png_byte green = 0;

        std::string colourKey = "colour" + i_to_str(count);
        if (properties.count(colourKey))
        {
            if (properties[colourKey] == "blue")
            {
                blue = 255;
            }
            if (properties[colourKey] == "red")
            {
                red = 255;
            }
            if (properties[colourKey] == "green")
            {
                green = 255;
            }
        }

        float transparency = 0.0;
        std::string transKey = "transparency" + i_to_str(count);
        if (properties.count(transKey))
        {
            transparency = atof(properties[transKey].c_str());
        }

        std::vector<Coord> points;

        for (int i = 0; i < entries.size(); i++)
        {
            if (entries[i][xCol] != entries[i][xCol] || entries[i][yCol] != entries[i][yCol])
            {
                continue;
            }

            points.push_back(std::make_pair(entries[i][xCol], entries[i][yCol]));
        }

        if (style != GraphStyleScatter)
        {
            std::sort(points.begin(), points.end(), std::less<Coord>());
        }

        if (points.size() == 0)
        {
            count++;
            continue;
        }

        double lastX = (points[0].first - minX) / (maxX - minX);
        double lastY = (points[0].second - minY) / (maxY - minY);

        for (int i = 0; i < points.size(); i++)
        {
            double x = points[i].first;
            double y = points[i].second;

            double xProp = (x - minX) / (maxX - minX);
            double yProp = (y - minY) / (maxY - minY);

            if (xProp > 1 || yProp > 1 || xProp < 0 || yProp < 0)
            {
                continue;
            }

            xProp *= xAxis * width;
            yProp *= -yAxis * height;

            if (style == GraphStyleScatter)
            {
                const int length = 1;
                png->drawLine(xProp - length, yProp - length, xProp + length, yProp + length, transparency, red, green, blue);
                png->drawLine(xProp - length, yProp + length, xProp + length, yProp - length, transparency, red, green, blue);
            }
            else if (style == GraphStyleLine)
            {
                png->drawLine(lastX, lastY, xProp, yProp, transparency, red, green, blue);
                lastX = xProp;
                lastY = yProp;
            }
        }

        count++;
    }

    png->writeImageOutput();
}

std::string CSV::plotColumns(int col1, int col2)
{
    char *cols = std::getenv("COLUMNS");
    char *linesstr = getenv("LINES");

 //   std::cout << "Plotting " << entries.size() << " points" << std::endl;

    int columns = 100;

    if (cols && strlen(cols))
    {
        columns = atoi(cols);
    }

    int lines = 50;

    if (cols && strlen(linesstr))
    {
        lines = atoi(linesstr);
    }

    if (columns < 80) columns = 80;
    if (lines < 50) lines = 50;

    const int leftMargin = 8;
    const int bottomMargin = 4;

    int graphCols = columns - leftMargin;
    int graphLines = lines - bottomMargin;

    Plot plot;

    for (int i = 0; i < lines; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            plot[i][j] = PlotBlank;
        }
    }

    double xMin, xMax, yMin, yMax;

    if (!didSetMinMaxXY)
    {
        minMaxCol(col1, &xMin, &xMax);
        minMaxCol(col2, &yMin, &yMax);
    }
    else
    {
        xMin = minX;
        xMax = maxX;
        yMin = minY;
        yMax = maxY;
    }

    for (int i = 0; i < graphLines; i++)
    {
        plot[i][leftMargin] = PlotVerticalLine;
    }

    for (int i = 0; i < graphLines; i += ((double)graphLines / 10.))
    {
        plot[i][leftMargin - 1] = PlotHorizontalTickMark;
        std::ostringstream valueStream;
        double value = ((double)i / (double)graphLines) * (yMin - yMax) + yMax;

        valueStream << std::fixed << std::setprecision(2) << value;

        writeStringToPlot(valueStream.str(), &plot, i, 0);
    }

    for (int i = 0; i < graphCols; i++)
    {
        plot[graphLines][leftMargin + i] = PlotHorizontalLine;
    }

    writeStringToPlot(headers[col1], &plot, graphLines + bottomMargin - 1, graphCols / 2);

    for (int i = 0; i < graphCols; i += ((double)graphCols / 10.))
    {
        plot[graphLines + 1][leftMargin + i] = PlotVerticalTickMark;

        std::ostringstream valueStream;
        double value = ((double)i / (double)graphCols) * (xMax - xMin) + xMin;

        valueStream << std::fixed << std::setprecision(2) << value;

        writeStringToPlot(valueStream.str(), &plot, graphLines + 2, leftMargin + i - 2);
    }

    for (int i = 0; i < entries.size(); i++)
    {
        double xValue = entries[i][col1];
        double yValue = entries[i][col2];

        if (xValue < xMin || xValue > xMax || yValue < yMin || yValue > yMax)
            continue;

        double xProportion = (xValue - xMin) / (xMax - xMin);
        double yProportion = (yValue - yMin) / (yMax - yMin);

        int xCharsIn = xProportion * (double)(graphCols - 1);
        int yCharsIn = yProportion * (double)graphLines;

        int xFullCharsIn = xCharsIn + leftMargin + 2;
        int yFullCharsIn = (graphLines - yCharsIn);

        if (yFullCharsIn == graphLines)
            yFullCharsIn -= 1;

        char currentChar = plot[yFullCharsIn][xFullCharsIn];

        if (currentChar == ' ')
        {
            plot[yFullCharsIn][xFullCharsIn] = '.';
        }
        else if (currentChar == '.')
        {
            plot[yFullCharsIn][xFullCharsIn] = ':';
        }
        else if (currentChar == ':')
        {
            plot[yFullCharsIn][xFullCharsIn] = '^';
        }
        else if (currentChar == '^')
        {
            plot[yFullCharsIn][xFullCharsIn] = '*';
        }
        else if (currentChar == '*')
        {
            plot[yFullCharsIn][xFullCharsIn] = 'x';
        }
        else if (currentChar == 'x')
        {
            plot[yFullCharsIn][xFullCharsIn] = '#';
        }
        else if (currentChar == '#')
        {
            plot[yFullCharsIn][xFullCharsIn] = '%';
        }
        else if (currentChar == '%')
        {
            plot[yFullCharsIn][xFullCharsIn] = '@';
        }
    }

    std::ostringstream ascii;
    ascii << std::endl << headers[col2] << std::endl << std::endl;

    ascii << mapToAscii(plot);

    std::ostringstream logged;
    logged << ascii.str() << std::endl;
    Logger::mainLogger->addStream(&logged);

    return ascii.str();
}

void CSV::setValueForEntry(int entry, std::string header, double value)
{
    int column = findHeader(header);
    entries[entry][column] = value;
}

void CSV::addConvolutedPeak(int column, double mean, double stdev, double weight, int categoryNum, ConvolutionType type)
{
    double totalIntervals = type == ConvolutionTypeSuperGaussian ? 300 : 10;
    double stdevMult = type == ConvolutionTypeSuperGaussian ? 10 : 1;
    double step = (stdev * stdevMult) / totalIntervals;

    for (double x = -stdev * stdevMult / 2; x < stdev * stdevMult / 2; x += step)
    {
        double y = 1;

        switch (type) {
            case ConvolutionTypeSuperGaussian:
                y = super_gaussian(x, 0, stdev / 3, 1.0);
                break;
            case ConvolutionTypeUniform:
                y = 1;
                break;
            default:
                break;
        }

        addOneToFrequency(x + mean, column, y * weight, categoryNum);
    }
}

void CSV::convolutedPeaks(std::string category, std::string origHeader, std::string destHeader, double stdev, ConvolutionType type)
{
    int origNum = findHeader(origHeader);
    int categoryNum = findHeader(category);
    int column = findHeader(destHeader);

    if (type == ConvolutionTypeSuperGaussian)
    {
        for (int i = 0; i < entries.size(); i++)
        {
            double mean = entries[i][categoryNum];
            double weight = entries[i][origNum];

            addConvolutedPeak(column, mean, stdev, weight, categoryNum, type);
        }
    }
    else if (type == ConvolutionTypeUniform)
    {
        int stdevInt = (int)stdev;
        for (int i = stdevInt; i < entries.size() - stdevInt; i++)
        {
            for (int j = -stdevInt; j <= stdevInt; j++)
            {
                double weight = entries[i][origNum];

                entries[i + j][column] += weight / (stdevInt * 2 + 1);
            }
        }
    }
}

void CSV::resetColumn(std::string header, double value)
{
    int headerNum = findHeader(header);

    for (int i = 0; i < entries.size(); i++)
    {
        entries[i][headerNum] = value;
    }
}

CSV::~CSV()
{
    headers.clear();
    entries.clear();

    std::vector<std::string>().swap(headers);
    std::vector<Entry>().swap(entries);
}
