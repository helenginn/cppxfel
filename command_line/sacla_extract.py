#!/usr/bin/env libtbx.python

import h5py
import cPickle as pickle
import sys
import scitbx_array_family_flex_ext
import numpy
import os
from multiprocessing import Process
from scitbx.array_family.flex import grid
from multiprocessing import Lock

"""
cxiKey is a h5py file key. allKeys is a list of keys found.
"""
def findTags(cxiKey, allKeys):
        subKeys = None
        try:
                subKeys = cxiKey.keys()

                for subKey in subKeys:
                        if subKey[:4] == "tag-":
                                if subKey[5].isnumeric():
                                        allKeys.append(cxiKey[subKey])
                                        print "Found key", subKey
                        else:
                                findTags(cxiKey[subKey], allKeys)
        except Exception, e:
                pass


filename = ""
geometry = ""
beamX = 1100.
beamY = 1100.
width = 2200
height = 2200
skip = 0
sacla = False
max = 0
wavelength = 1.24
distance = 70

print "CXI to pickle dumper"

for arg in sys.argv[1:]:
        key_value = arg.split("=")
        if (len(key_value) < 2):
                print "Did not understand", arg
                continue

        if (key_value[0] == "filename"):
                filename = key_value[1]
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
        elif (key_value[0] == "skip"):
                skip = int(key_value[1])
        elif (key_value[0] == "max"):
                max = int(key_value[1])
        elif (key_value[0] == "wavelength"):
                wavelength = int(key_value[1])
        elif (key_value[0] == "distance"):
                distance = int(key_value[1])
        else:
                print "Did not understand", arg

if len(filename) == 0:
        print "No filename selected, use filename=xxx.cxi"
        exit()

if len(geometry) == 0:
        print "No geometry file selected, use geometry=xxx.geom"
        exit()

print "Selected filename:", filename
print "Applying geometry:", geometry

cxi = h5py.File(filename, 'r')

entries = []

findTags(cxi, entries)

print "Found", len(entries), "images."

binnedEntries = {}

for tag in entries:
        tagNameOnly = tag.name.split("/")[-1]
        if tagNameOnly in binnedEntries:
                binnedEntries[tagNameOnly].append(tag)
        else:
                binnedEntries[tagNameOnly] = []

print "Found", len(binnedEntries), "unique images."

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
                break

if len(panels) == 0:
        panelString = "q1,q2,q3,q4,q5,q6,q7,q8"
        print "Could not find panels. Will default to", panelString
        panels = panelString.strip().split(",")

print "This group has", len(panels), "panels."

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

lock = Lock()

class EntryDumper():
        def __init__(self):
                self.lock = Lock()

        def dumpEntry(self, anEntry):
                global beamX, beamY, width, height, wavelength, distance
                num_images = 1

                sample_file = open('sample.pickle', 'rb')
                sample = pickle.load(sample_file)

                realEntryMin = 0
                realEntryMax = 1

                pixelSize = 0.05

                identifier = anEntry.name[1:]
                print "Dumping", identifier

                totalPicklePixels = width * height
                print "Generating canvas for pickle image with", totalPicklePixels, "pixels."

                picklePixels = scitbx_array_family_flex_ext.int()
                picklePixels.resize(totalPicklePixels)
                count = 0

                for panel in panelInfo:
                        count += 1

                        test = "Test"
                        self.lock.acquire()
                        print "Lock acquired"
                        try:
                                dataEntry = cxi[anEntry.name + "/data"]
                                image = dataEntry
                                rowSize = len(image[0])
                                alldata = []
                                alldata = numpy.concatenate(image[:])
                        finally:
                                self.lock.release()
                                print "Lock released"


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

                                        try:
                                                newValue = int(alldata[j * rowSize + k])
                                        except IndexError:
                                                print "Geometry file specifies coordinates out of range of hdf5 data."
                                                exit()

                                        try:
                                                picklePixels[newX * width + newY] = newValue
                                        except IndexError:
                                                print "Error: image canvas is not large enough."
                                                print "Tried to access point", newX, newY, "- please increase canvas size."
                                                print "Specify width=x height=y on the command line."
                                                exit()

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

                picklePixels.reshape(grid(width, height))

                sample['DATA'] = picklePixels

                pickleName = identifier + '.pickle'
                new_pickle = open(pickleName, 'wb')
                pickle.dump(sample, new_pickle, protocol=-1)
                print "Dumped pickle", pickleName, "and img", identifier + ".img"

                newFile.close()
                new_pickle.close()

#               compressCommand = "cxi.image2pickle " + pickleName
#               os.system(compressCommand)

        def dumpEntryThread(self, threadNum):
                global realMin, realMax, entries
                for i in range(realMin + threadNum, realMax, maxThreads):
                        print "Dumping entry", i

                        self.dumpEntry(entries[i])

maxThreads = int(os.getenv('NSLOTS', 4))

threads = []
dumper = EntryDumper()

for thread in range(0, maxThreads):
        aThread = Process(target=dumper.dumpEntryThread, args=(thread, ))
        threads.append(aThread)
        aThread.start()

print "Done."
