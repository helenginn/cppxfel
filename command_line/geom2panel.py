#!/usr/bin/env libtbx.python

import sys
import numpy
import os
import math

geometry = sys.argv[1]

rigidGroupString = "rigid_group_collection_connected"
#print "Examining geometry file: finding", rigidGroupString

panels = []
geomfile = open(geometry, 'r')
beginning = geomfile.tell()

for line in geomfile:
    if line[0] == ";":
        continue
    components = line.split("=")
    if (components < 2):
        continue
    if components[0].strip() == rigidGroupString:
#        print "Found", rigidGroupString
        panels = components[1].strip().split(",")
        break

if len(panels) == 0:
    panelString = "q0a0,q0a1,q0a2,q0a3,q0a4,q0a5,q0a6,q0a7,q0a8,q0a9,q0a10,q0a11,q0a12,q0a13,q0a14,q0a15,q1a0,q1a1,q1a2,q1a3,q1a4,q1a5,q1a6,q1a7,q1a8,q1a9,q1a10,q1a11,q1a12,q1a13,q1a14,q1a15,q2a0,q2a1,q2a2,q2a3,q2a4,q2a5,q2a6,q2a7,q2a8,q2a9,q2a10,q2a11,q2a12,q2a13,q2a14,q2a15,q3a0,q3a1,q3a2,q3a3,q3a4,q3a5,q3a6,q3a7,q3a8,q3a9,q3a10,q3a11,q3a12,q3a13,q3a14,q3a15"

    panels = panelString.strip().split(",")

#print "This group has", len(panels), "panels."

panelInfo = {}
orderedPanels = []

for panel in panels:
    geomfile.seek(beginning)
    orderedPanels.append({})
    for line in geomfile:
        if (line[:len(panel)] == panel and line[len(panel):len(panel)+1] == "/"):
            key_value = line.split("=")
            key = key_value[0].strip()
            value = key_value[1].strip()

            id_property = key.split("/")
            id = id_property[0]
            property = id_property[1]

            if len(value.split(" ")) > 1:
                subProperties = value.split(" ")
                for i in range(len(subProperties)):
                    subProperties[i] = subProperties[i].replace("x", "").replace("y", "")
                    subProperties[i] = float(subProperties[i])
                value = subProperties
            else:
                value = float(value)

#            print "Found value", value, "for panel", panel, ", property", property
            orderedPanels[-1][property] = value

geomfile.close()

for panel in orderedPanels:
        origTopLeftX = panel['min_fs']
        origTopLeftY = panel['min_ss']
        origBottomRightX = panel['max_fs']
        origBottomRightY = panel['max_ss']
        relativeNewX = panel['corner_x']
        relativeNewY = panel['corner_y']
        axisXX = panel['fs'][0]
        axisXY = panel['fs'][1]
        axisYX = panel['ss'][0]
        axisYY = panel['ss'][1]

        # generate panels.txt

        imageWidth = origBottomRightX - origTopLeftX
        imageHeight = origBottomRightY - origTopLeftY
        relativeNewOtherCornerX = relativeNewX + imageWidth * axisXX + imageHeight * axisYX
        relativeNewOtherCornerY = relativeNewY + imageWidth * axisXY + imageHeight * axisYY
        physicalMiddleX = (relativeNewX + relativeNewOtherCornerX) / 2.0
        physicalMiddleY = (relativeNewY + relativeNewOtherCornerY) / 2.0
        imageMiddleX = origTopLeftX + imageWidth / 2.0
        imageMiddleY = origTopLeftY + imageHeight / 2.0
        offsetX = imageMiddleX - physicalMiddleX
        offsetY = imageMiddleY - physicalMiddleY
        rotation = - math.atan2(axisXY, axisXX) / math.pi * 180.0

#        print "#", panel
        print "PANEL %f %f %f %f %f %f 0 0 %f 1" % (min(relativeNewX, relativeNewOtherCornerX), min(relativeNewY, relativeNewOtherCornerY),
                                                    max(relativeNewX, relativeNewOtherCornerX), max(relativeNewY, relativeNewOtherCornerY),
                                                    offsetX, offsetY, rotation)
