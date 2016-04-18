from xfel.cxi.cspad_ana.cspad_tbx import dpack
import numpy as np
import h5py
import sys
import cPickle as pickle
from dials.array_family import flex
from dxtbx.format.FormatPYunspecifiedStill import *
from dxtbx.imageset import MemImageSet
from dxtbx.datablock import DataBlockFactory

from libtbx.phil import parse
from dials.util.options import OptionParser
phil_scope = parse('''
include scope dials.algorithms.spot_finding.factory.phil_scope
''', process_includes=True)
parser = OptionParser(phil=phil_scope)
params, options, all_paths = parser.parse_args(show_diff_phil=False, return_unhandled=True)

inp = all_paths[0]

f = h5py.File(inp, "r")
tags = [tag for tag in f if tag.startswith("tag-")]
arr = np.zeros((8192, 512), dtype='uint32')

for tag in tags:
    print "Processing " + tag,

    f[tag + "/data"].read_direct(arr)

    data = dpack(data=flex.int(arr),
                 distance=50.0,
                 pixel_size=0.05,
                 wavelength=1.77,
                 beam_center_x=0,
                 beam_center_y=0,
                  # distance, pixel_size, wavelength, beam_center are not used by 
                  # dials.find_spots unless you specify resolution cutoff
                 ccd_image_saturation=65535,
                 saturated_value=65535,
                 timestamp=None
                )


    img = FormatPYunspecifiedStillInMemory(data)
    # This writes "This is not a file; assume the data are in the defined dictionary format", which is true.
    # I want to get rid of the message.
    imgset = MemImageSet([img])
    datablock = DataBlockFactory.from_imageset(imgset)[0]
    observed = flex.reflection_table.from_observations(datablock, params)

    #pickle.dump(data, open(tag + ".pickle", "wb"), protocol=-1)
    pickle.dump(observed, open("_" + tag + "_strong.pickle", "wb"), protocol=-1)
    print " found %d spots" % len(observed)
    