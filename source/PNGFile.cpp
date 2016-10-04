//
//  Image.cpp
//  CellCounter
//
//  Created by Helen Ginn on 07/08/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "PNGFile.h"
#include <string>
#include <stdlib.h>
#include <iostream>
#include "misc.h"
#include <iomanip>
#include "FileReader.h"
#include <cmath>
#include "Shoebox.h"

int PNGFile::writeImage(std::string filename, int width, int height, std::string title)
{
    int code = 0;
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    
    // Open file for writing (binary mode)
    fp = fopen(filename.c_str(), "wb");
    if (fp == NULL) {
        fprintf(stderr, "Could not open file %s for writing\n", filename.c_str());
        code = 1;
        goto finalise;
    }
    
    // Initialize write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }
    
    // Initialize info structure
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }
    
    // Setup Exception handling
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }
    
    png_init_io(png_ptr, fp);
    
    // Write header (8 bit colour depth)
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    
    // Set title
    if (title.length() > 0) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = (char *)"Title";
        title_text.text = (char *)title.c_str();
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }
    
    png_write_info(png_ptr, info_ptr);
    
    // Write image data

    for (int y = 0 ; y < height ; y++)
    {
        png_bytep row = &data[y * pixelsPerRow * bytesPerPixel];
        png_write_row(png_ptr, row);
    }
    
    // End write
    png_write_end(png_ptr, NULL);
    
finalise:
    if (fp != NULL)
    {
        fclose(fp);
    }
        
    if (info_ptr != NULL)
    {
        png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    }
    
    if (png_ptr != NULL)
    {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);

    }
    
    return code;
}
/*
void PNGFile::read_png_file(const char *file_name)
{
    png_byte header[8];
    
    // open file and test it for begina  png
    FILE *fp = fopen(file_name, "rb");
    if (!fp)
    {
        // BAD! handle
        std::cout << "Failed to open the file. Are you in the right directory?" << std::endl;
        return;
    }
    fread(header, 1, 8, fp);
    if (png_sig_cmp(header, 0, 8))
    {
        // BAD! handle
        std::cout << std::endl << "File does not appear to be a png file." << std::endl;
        return;
    }
    
    // initialize stuff

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!png_ptr)
    {
        std::cout << "Could not allocate space for PNG." << std::endl;
        return;
        // HANDLE
    }
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        std::cout << "Could not allocate space for INFO." << std::endl;
        return;
//        you know
    }
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        std::cout << "Something went wrong. Really a PNG file?" << std::endl;
        return;
        // you know
    }
    
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);
    
    png_read_info(png_ptr, info_ptr);
    
    png_uint_32 width = png_get_image_width(png_ptr, info_ptr);
    png_uint_32 height = png_get_image_height(png_ptr, info_ptr);
    png_byte color_type = png_get_color_type(png_ptr, info_ptr);
//    png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    
//    int number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);
    
    // read file
    
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        std::cout << "Something went wrong. Really a PNG file?" << std::endl;
        return;
        // HANDLE
    }
    
    png_uint_32 rowBytes = png_get_rowbytes(png_ptr, info_ptr);
    
    png_bytep *row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (int y = 0; y < height; y++)
    {
          row_pointers[y] = (png_bytep)malloc(png_get_rowbytes(png_ptr, info_ptr));
    }
    
    png_read_image(png_ptr, row_pointers);

    int offset = 0;
    
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < rowBytes; x++)
        {
            data[y][x] = row_pointers[y][x];
        }
        
        offset += height;
    }
    
    bytesPerPixel = (double)rowBytes / (double)width;
    
    if (color_type == 3)
    {
        std::cout << "ABORTING - this program does not support palette-based PNG files." << std::endl;
        exit(1);
    }
    
    if (color_type == 2)
        bytesPerPixel = 3;
    if (color_type == 6)
        bytesPerPixel = 4;
}*/

PNGFile::PNGFile(std::string filename, int width, int height)
{
    // set data using libpng...
    // data =
    
    this->height = height;
    bytesPerPixel = 3; // change
    pixelsPerRow = width; // change
    this->filename = filename;
    
    size_t arraySize = sizeof(png_byte) * bytesPerPixel * height * pixelsPerRow;
    
    data = (png_bytep)malloc(arraySize);
    memset(data, 255, arraySize);
    
    rootName = getBaseFilename(filename);
}

void PNGFile::writeImageOutput()
{
    std::string title = rootName;
    std::string file = FileReader::addOutputDirectory(filename);
    
    writeImage(file, pixelsPerRow, height, title);
}

png_byte PNGFile::valueAt(int x, int y)
{
    png_byte *bytes = new png_byte[bytesPerPixel];
    pixelAt(x, y, &bytes);
    int value = 0;
    
    for (int i = 0; i < bytesPerPixel; i++)
    {
        value += bytes[i];
    }
    
    value /= bytesPerPixel;
    
    return value;
}

void PNGFile::pixelAt(int x, int y, png_byte **bytes)
{
    int offset = y * bytesPerPixel * pixelsPerRow +  x * bytesPerPixel;
    
    *bytes = &data[offset];
    
}

void PNGFile::moveCoordRelative(int *x, int *y)
{
    *x += pixelsPerRow / 2 - centreX;
    *y += height / 2 - centreY;
}

void PNGFile::setPixelColourRelative(int x, int y, png_byte red, png_byte green, png_byte blue)
{
    moveCoordRelative(&x, &y);
    
    setPixelColour(x, y, red, green, blue);
}

void PNGFile::setPixelColour(int x, int y, png_byte red, png_byte green, png_byte blue, float transparency)
{
    int offset = y * bytesPerPixel * pixelsPerRow +  x * bytesPerPixel;
    
    if (offset < 0 || offset >= bytesPerPixel * pixelsPerRow * height)
    {
        return;
    }
    
    if (transparency < 1)
    {
        png_byte currentRed = data[offset];
        png_byte currentGreen = data[offset + 1];
        png_byte currentBlue = data[offset + 2];
        
        red = currentRed + (red - currentRed) * transparency;
        green = currentGreen + (green - currentGreen) * transparency;
        blue = currentBlue + (blue - currentBlue) * transparency;
    }
    
    data[offset] = red;
    data[offset + 1] = green;
    data[offset + 2] = blue;
}

void PNGFile::setPixelForChannel(int x, int y, int channel, png_byte value)
{
    int offset = y * bytesPerPixel * pixelsPerRow +  x * bytesPerPixel + channel;
    
    data[offset] = value;
}

void PNGFile::drawCircleAroundPixel(int x, int y, float radius, float transparency, png_byte red, png_byte green, png_byte blue, float thickness)
{
    moveCoordRelative(&x, &y);
    double radiusSqr = radius * radius;
    double minDiff = radiusSqr - pow(radius - thickness / 2, 2);
    
    for (int i = x - radius; i <= x + radius; i++)
    {
        for (int j = y - radius; j <= y + radius; j++)
        {
            double distSqr = pow(i - x, 2) + pow(j - y, 2);
            double diff = fabs(distSqr - radiusSqr);
            
            if (diff < minDiff)
            {
                setPixelColour(i, j, red, green, blue, transparency);
            }
        }
    }
}

void PNGFile::drawShoeboxAroundPixel(int x, int y, ShoeboxPtr shoebox)
{
    moveCoordRelative(&x, &y);
    int fastLength, slowLength;
    int fastCentre, slowCentre;
    shoebox->sideLengths(&slowLength, &fastLength);
    shoebox->centre(&fastCentre, &slowCentre);
   
    for (int i = 0; i < slowLength; i++) // y
    {
        std::vector<double> row = shoebox->row(i);
        
        int slowOffset = i - slowCentre;
        
        for (int j = 0; j < fastLength; j++) // x
        {
            int fastOffset = j - fastCentre;
            
            double value = row[j];
            
            int xPos = x + fastOffset;
            int yPos = y + slowOffset;
            
            if (value < 0)
            {
                setPixelColour(xPos, yPos, 255, 0, 0, 0.5);
            }
            else if (value > 0)
            {
                setPixelColour(xPos, yPos, 0, 0, 255, 0.5);
            }
        }
    }
}

void PNGFile::dropImage()
{
    free(data);
}
