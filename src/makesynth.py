#
#	makesynth.py
#
#	Makes synthetic HST images using tinytim psfs.
#
#	Benjamin Proudfoot
#	01/14/22
#

### THIS CODE RUNS WITH TEST DATA ###

import numpy as np
import sys
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import commentjson as json
from tqdm import tqdm
from getpsf import *
from analysis import *
import tinytim_psfs.make_psf
import os.path
import datetime
from astropy.io import fits
from insertpsf import *


class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

### Load run_props from JSON file ###
runprops = ReadJson("synthprops.txt").outProps()
output_filename = runprops.get("output_filename")

# Getting inputs for Tiny Tim PSFs
focus = runprops.get("focus")
xpos = runprops.get("xpos")
ypos = runprops.get("ypos")
size_psf = runprops.get("psf_size")
sample_factor = runprops.get("sample_factor")
nchip = runprops.get("chip")
ndet = runprops.get("det_int")
filter = runprops.get("filter")

# Making TinyTim PSFs
filename = "modelpsfs/wfc3psf_" + str(ndet) + "_" + str(nchip) + "_" + filter + "_" + str(xpos) + "_" + str(ypos) + "_" + str(round(focus,1)) + "_" + str(size_psf) + "_" + str(sample_factor) + ".fits"
if not os.path.exists(filename):
	tinytim_psfs.make_psf.make_subsampled_model_psf(filename, psf_position = (xpos, ypos), focus = focus, chip = nchip, detector = ndet,
					   filter_name = filter, psf_size = size_psf, subsampling_factor = sample_factor)
else:
	print("File exists")

# Load in one PSF to check the size and make array for PSFs
ttpsf = getpsf_hst(filename)

# Pulling parameters from synthprops
image_size = runprops.get("image_size")
skynoise = runprops.get("skynoise")
psf_x = np.array(runprops.get("psf_x"))
psf_y = np.array(runprops.get("psf_y"))
psf_h = np.array(runprops.get("psf_h"))
background = np.random.poisson(lam = skynoise, size = (image_size,image_size))

# Creating and saving synthetic image
synthimage = insertpsf_n(image = background, psf = ttpsf, xcens = psf_x, ycens = psf_y, heights = psf_h, psfscale = runprops.get("sample_factor"))
if os.path.exists('../data/' + runprops.get("output_name") + ".fits"):
	os.remove('../data/' + runprops.get("output_name") + ".fits")
hdu = fits.PrimaryHDU(synthimage)
hdu.writeto('../data/' + runprops.get("output_name") + ".fits")
