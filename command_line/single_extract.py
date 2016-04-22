#!/usr/bin/env libtbx.python

import h5py
import cPickle as pickle
import sys
import scitbx_array_family_flex_ext
import numpy
import os
from multiprocessing import Process
from scitbx.array_family.flex import grid


filename = ""
geometry = ""
entry = "data"
beamX = 882.
beamY = 882.
width = 1765
height = 1765
distance = 100
wavelength = 1.4
skip = 0
sacla = False
max = 0

print "CXI to pickle dumper"

for arg in sys.argv[1:]:
	key_value = arg.split("=")
	if (len(key_value) < 2):
		print "Did not understand", arg
		continue
	
	if (key_value[0] == "filename"):
		filename = key_value[1]
	elif (key_value[0] == "entry"):
		entry = key_value[1]
	elif (key_value[0] == "geometry"):
		geometry = key_value[1]
	elif (key_value[0] == "beamx"):
		beamX = float(key_value[1])
	elif (key_value[0] == "beamy"):
		beamY = float(key_value[1])
	elif (key_value[0] == "width"):
		width = int(key_value[1])
	elif (key_value[0] == "height"):
		height = int(key_value[1])
	elif (key_value[0] == "distance"):
		distance = float(key_value[1])
	elif (key_value[0] == "wavelength"):
		wavelength = float(key_value[1])
	elif (key_value[0] == "skip"):
		skip = int(key_value[1])
	elif (key_value[0] == "sacla"):
		sacla = (key_value[1] == "True")
	elif (key_value[0] == "max"):
		max = int(key_value[1])
	else:
		print "Did not understand", arg

if len(filename) == 0:
	print "No filename selected, use filename=xxx.cxi"
	exit()
	
if len(geometry) == 0:
	print "No geometry file selected, use geometry=xxx.geom"
	exit()
	
print "Selected filename:", filename
print "Dumping entry:", entry
print "Applying geometry:", geometry

cxi = h5py.File(filename, 'r')

entries = []

if not sacla:
	entries.append(entry + "/data")
else:
	for key in cxi.keys():
		entries.append(key)

print "All keys:"
print entries

rigidGroupString = "rigid_group_d0"
if sacla:
	rigidGroupString = "rigid_group_collection_connected"

print "Examining geometry file: finding", rigidGroupString

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
		print "Found", rigidGroupString
		panels = components[1].strip().split(",")
		print "This group has", len(panels), "panels."
		break

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
		
			print "Found value", value, "for panel", panel, ", property", property
			panelInfo[panel][property] = value

geomfile.close()

realMin = skip
realMax = max + skip

if max == 0:
	realMax = len(entries)

if (skip + max) > len(entries):
	realMax = len(entries)

if not sacla:
	realMin = 0
	realMax = 1

anEntry = entries[0]
print anEntry
dataEntry = anEntry

if (sacla):
	dataEntry = anEntry + "/data"

num_images = 1

if sacla:
	num_images = 1

print "\nEntry has", num_images, "images.\n"
sample_file = open('sample.pickle', 'rb')
sample = pickle.load(sample_file)

for i in range(0, 1):
	image = None
	identifier = anEntry
	pixelSize = 0.05
	if not sacla:
		image = cxi[entry + "/data"]
		identifierList = filename.split(".")[:-1]
		identifier = "".join(identifierList)
		
		if identifier[0] == "/":
			identifier = identifier[1:]
	else:
		image = cxi[anEntry + "/data"]
		wavelength = cxi[anEntry + "/photon_wavelength_A"].value
	
	print "Dumping", identifier
	alldata = []
	count = 0
	rowSize = len(image[0])
	if sacla:
		alldata = numpy.concatenate(image[:])
	else:
		alldata = numpy.concatenate(image[:])
	totalPicklePixels = width * height
	print "Generating canvas for pickle image with", totalPicklePixels, "pixels."

	picklePixels = scitbx_array_family_flex_ext.int()
	picklePixels.resize(totalPicklePixels)

	for panel in panelInfo:
		origTopLeftX = int(panelInfo[panel]['min_fs'])
		origTopLeftY = int(panelInfo[panel]['min_ss'])
		origBottomRightX = int(panelInfo[panel]['max_fs'])
		origBottomRightY = int(panelInfo[panel]['max_ss'])
		relativeNewX = panelInfo[panel]['corner_x']
		relativeNewY = panelInfo[panel]['corner_y']
		absoluteNewX = int(beamX + relativeNewX)
		absoluteNewY = int(beamY + relativeNewY)
		axisXX = int(round(panelInfo[panel]['fs'][0]))
		axisXY = int(round(panelInfo[panel]['fs'][1]))
		axisYX = int(round(panelInfo[panel]['ss'][0]))
		axisYY = int(round(panelInfo[panel]['ss'][1]))

		print panel + "\t" + str(origTopLeftX) + "\t" + str(origTopLeftY) + "\t" + str(origBottomRightX) + "\t" + str(origBottomRightY) + "\t" + str(absoluteNewX) + "\t"  + str(absoluteNewY) + "\t" + str(axisXX) + "\t" + str(axisXY) + "\t" + str(axisYX) + "\t" + str(axisYY)

		for k in range(origTopLeftX, origBottomRightX):
			xOffset = k - origTopLeftX
			for j in range(origTopLeftY, origBottomRightY):
				yOffset = j - origTopLeftY
				newX = absoluteNewY + axisXX * yOffset + axisXY * xOffset
				newY = absoluteNewX + axisYX * yOffset + axisYY * xOffset
		
				newValue = int(alldata[j * rowSize + k])
				picklePixels[newX * width + newY] = newValue

	string = picklePixels.as_numpy_array().tostring()
	newFile = open(identifier + '.img', 'wb')
	newFile.write(string)
	newFile.close()

	sample['DISTANCE'] = distance
	sample['SIZE1'] = width
	sample['SIZE2'] = height
	sample['ACTIVE_AREAS'] = scitbx_array_family_flex_ext.int([0, 0, width, height])
	sample['PIXEL_SIZE'] = pixelSize
	sample['BEAM_CENTER_X'] = width / 2 * pixelSize
	sample['BEAM_CENTER_Y'] = height / 2 * pixelSize
	sample['WAVELENGTH'] = wavelength
# 
	picklePixels.reshape(grid(width, height))
# 
	sample['DATA'] = picklePixels

	pickleName = identifier + '.pickle'
	new_pickle = open(pickleName, 'wb')
	pickle.dump(sample, new_pickle, protocol=-1)
	print "Dumped pickle", pickleName, "and img", identifier + ".img"

	newFile.close()
	new_pickle.close()


print "Done."