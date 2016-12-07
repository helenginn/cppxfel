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
#		print "Found", rigidGroupString
		panels = components[1].strip().split(",")
		break

if len(panels) == 0:
	panelString = "q1,q2,q3,q4,q5,q6,q7,q8"
#	print "Could not find panels. Will default to", panelString
	panels = panelString.strip().split(",")

#print "This group has", len(panels), "panels."

panelInfo = {}

for panel in panels:
	geomfile.seek(beginning)
	panelInfo[panel] = {}
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
	
#			print "Found value", value, "for panel", panel, ", property", property
			panelInfo[panel][property] = value

geomfile.close()

for panel in panelInfo:
        origTopLeftX = panelInfo[panel]['min_fs']
        origTopLeftY = panelInfo[panel]['min_ss']
        origBottomRightX = panelInfo[panel]['max_fs']
        origBottomRightY = panelInfo[panel]['max_ss']
        relativeNewX = panelInfo[panel]['corner_x']
        relativeNewY = panelInfo[panel]['corner_y']
        axisXX = panelInfo[panel]['fs'][0]
        axisXY = panelInfo[panel]['fs'][1]
        axisYX = panelInfo[panel]['ss'][0]
        axisYY = panelInfo[panel]['ss'][1]
        
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

        print "#", panel
        print "PANEL %f %f %f %f %f %f 0 0 %f 1" % (min(relativeNewX, relativeNewOtherCornerX), min(relativeNewY, relativeNewOtherCornerY),
                                                    max(relativeNewX, relativeNewOtherCornerX), max(relativeNewY, relativeNewOtherCornerY),    
                                                    offsetX, offsetY, rotation)

