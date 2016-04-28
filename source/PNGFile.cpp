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

int PNGFile::writeImage(std::string filename, int width, int height, std::string title)
{
    int code = 0;
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_bytep row = NULL;
    
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
    
    // Allocate memory for one row (3 bytes per pixel - RGB)
    row = (png_bytep) malloc(3 * width * sizeof(png_byte));
    
    // Write image data
    int x, y;
    for (y=0 ; y<height ; y++) {
        for (x=0 ; x<width * bytesPerPixel; x++) {
            row[x] = data[y][x];
        }
        png_write_row(png_ptr, row);
    }
    
    // End write
    png_write_end(png_ptr, NULL);
    
finalise:
    if (fp != NULL) fclose(fp);
    if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    if (row != NULL) free(row);
    
    return code;
}

void PNGFile::read_png_file(const char *file_name)
{
    png_byte header[8];
    
    /* open file and test it for begina  png */
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
    
    /* initialize stuff */
    
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
    png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    
//    int number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);
    
    /* read file */
    
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
    
    columnCount = height;
    
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
    
    pixelsPerRow = width;
}

PNGFile::PNGFile(std::string filename)
{
    // set data using libpng...
    // data =
    
    height = 0;
    rejectCount = 40;
    bytesPerPixel = 4; // change
    pixelsPerRow = 500; // change
    columnCount = 500; // change
    brightness = false;
    channel = 2;
    this->filename = filename;
    
    std::cout << "Loading " << filename << std::endl;
    getImage();
    
    rootName = getBaseFilename(filename);
}

void PNGFile::writeImageOutput()
{
    std::string title = rootName;
    std::string file = getPath(filename) + "counted_" + rootName + ".png";
    
    writeImage(file, pixelsPerRow, height, title);
}

png_byte PNGFile::valueAt(int x, int y)
{
    getImage();
    
    png_byte *bytes = new png_byte[bytesPerPixel];
    pixelAt(x, y, &bytes);
    int value = 0;
    
    if (!brightness)
    {
        value = bytes[channel];
    }
    else
    {
        for (int i = 0; i < bytesPerPixel; i++)
        {
            value += bytes[i];
        }
        
        value /= bytesPerPixel;
    }
    
    return value;
}

void PNGFile::pixelAt(int x, int y, png_byte **bytes)
{
    int realX = x * bytesPerPixel;
    
    png_byte *pngBytes = new png_byte[bytesPerPixel];
    
    for (int i = 0; i < bytesPerPixel; i++)
    {
        pngBytes[i] = data[y][realX + i];
    }
    
    memcpy(*bytes, pngBytes, bytesPerPixel);
    
  //  std::cout << int((*bytes)[0]) << " " << int((*bytes)[1]) << " " << int((*bytes)[2]) << std::endl;

    delete [] pngBytes;
}

void PNGFile::setPixelsColour(CoordList coords, png_byte red, png_byte green, png_byte blue)
{
    for (int i = 0; i < coords.size(); i++)
    {
        setPixelColour(coords[i].first, coords[i].second, red, green, blue);
    }
}

void PNGFile::setPixelsForChannel(std::vector<Coord> coords, int channel, png_byte value)
{
    for (int i = 0; i < coords.size(); i++)
    {
        setPixelForChannel(coords[i].first, coords[i].second, channel, value);
    }
}

void PNGFile::setPixelColour(int x, int y, png_byte red, png_byte green, png_byte blue)
{
    int realX = x * bytesPerPixel;
    
    data[y][realX] = red;
    data[y][realX + 1] = green;
    data[y][realX + 2] = blue;
}

void PNGFile::setPixelForChannel(int x, int y, int channel, png_byte value)
{
    int realX = x * bytesPerPixel + channel;
    
    data[y][realX] = value;
}

Flag PNGFile::flagAt(int x, int y)
{
    if (flags.count(x) == 0 || flags[x].count(y) == 0)
        return FlagNone;
    else return flags[x][y];
}

void PNGFile::setFlag(int x, int y, Flag flag)
{
    flags[x][y] = flag;
}

void PNGFile::getImage()
{
    if (data.size() == 0)
        read_png_file(filename.c_str());
}

void PNGFile::dropImage()
{
    for (DataBlock::iterator it = data.begin(); it != data.end(); it++)
    {
        data[it->first].clear();
        DataRow().swap(data[it->first]);
    }
    
    data.clear();
    DataBlock().swap(data);
}
