#!/usr/bin/env libtbx.python

import sys
import numpy
import os
import math
import time

axisYX = None
axisYY = None
axisXX = None
axisXY = None
min_fs = None
min_ss = None
max_fs = None
max_ss = None

# line is a string.
# panel_num is an integer
def find_new_value_as_line(line, panel_num):
    global cpp_panels_lines, axisYY, axisYX, axisXX, axisXY
    global min_fs, min_ss, max_fs, max_ss
    global beamX, beamY
    backslash_split = line.split("/")
    equal_split = backslash_split[1].split("=")
    variable = equal_split[0].strip()

    if variable == "min_fs" or variable == "min_ss" or variable == "max_fs" or variable == "max_ss":
        if variable == "min_fs":
            min_fs = float(equal_split[1].strip(" "))
        if variable == "min_ss":
            min_ss = float(equal_split[1].strip(" "))
        if variable == "max_fs":
            max_fs = float(equal_split[1].strip(" "))
        if variable == "max_ss":
            max_ss = float(equal_split[1].strip(" "))

        return line

    line_beginning = backslash_split[0] + "/" + variable + " = "

    #         rotation = - math.atan2(axisXY, axisXX) / math.pi * 180.0
    if variable == "fs" or variable == "ss":
        # get this information from the rotational component
        rotation = float(cpp_panels_lines[panel_num].split(" ")[9])
        rotation = -rotation
        if variable == "ss":
            rotation += 90
        theta = math.radians(rotation)
        dx = math.cos(theta)
        dy = math.sin(theta)

        if variable == "ss":
            axisYX = dx
            axisYY = dy
        elif variable == "fs":
            axisXX = dx
            axisXY = dy

        sign_x = ""
        if dx > 0: sign_x = "+"
        sign_y = ""
        if dy > 0: sign_y = "+"

        line_beginning += sign_x + '{0:.10f}'.format(dx) + "x " + sign_y + '{0:.10f}'.format(dy) + "y"
        return line_beginning

    imageWidth = max_fs - min_fs
    imageHeight = max_ss - min_ss

    if variable == "corner_x" or variable == "corner_y":
        cpp_values = cpp_panels_lines[panel_num].split(" ")

        minRotatedX = float(cpp_values[1])
        maxRotatedX = float(cpp_values[3])
        minRotatedY = float(cpp_values[2])
        maxRotatedY = float(cpp_values[4])

        midX = (maxRotatedX + minRotatedX) / 2
        midY = (maxRotatedY + minRotatedY) / 2
        vecX = -(axisYX * imageHeight + axisXX * imageWidth) / 2
        vecY = -(axisYY * imageHeight + axisXY * imageWidth) / 2
        distance = math.sqrt(pow(imageWidth / 2, 2) + pow(imageHeight / 2, 2))
        vecDist = math.sqrt(pow(vecY, 2) + pow(vecX, 2))
        scale = distance / vecDist
        vecX *= scale
        vecY *= scale
        midX += vecX
        midY += vecY
        newValue = 0

        if variable == "corner_x":
            newValue = midX - beamX
        elif variable == "corner_y":
            newValue = midY - beamY

        return line_beginning + str(newValue)

    cppline = cpp_panels_lines[panel_num]

    return line

template_name = ""
cpp_panels_name = ""
beamX = 0
beamY = 0

new_geom = ""

try:
    template_name = sys.argv[1]
    cpp_panels_name = sys.argv[2]
    beamX = float(sys.argv[3])
    beamY = float(sys.argv[4])
except:
    print "cppxfel.panel2geom <CrystFEL template> <cppxfel panels file> <beamX> <beamY>"
    exit()

print "; Note: this geometry file been converted from cppxfel."
print "; The cppxfel script is extraordinarily hacky!"
print "; Hacky conversion occurred on:", time.strftime("%d/%m/%Y %H:%M")
print

template = open(template_name, "r").read()
cpp_panels = open(cpp_panels_name, "r").read()

template_lines = template.split("\n")
cpp_panels_lines = cpp_panels.split("\n")

current_panel = "non_panel"
panel_number = -1

for tmpline in template_lines:
    backslash_split = tmpline.split("/")

    if len(backslash_split) > 0 and backslash_split[0] == current_panel:
        new_line = find_new_value_as_line(tmpline, panel_number)
        print new_line
    elif (len(backslash_split) > 0 and len(backslash_split[0]) >= 4 and
          backslash_split[0][0] == "q" and backslash_split[0][2] == "a"):
            current_panel = backslash_split[0]
            panel_number += 1
            new_line = find_new_value_as_line(tmpline, panel_number)
            print new_line
    else:
        print tmpline
