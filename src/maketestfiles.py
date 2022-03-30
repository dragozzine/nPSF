#
#	maketestfiles.py
#
#	Makes test files for nPSF. This loops over many of the settings internal to nPSF
#	and makes synthetic test cases to rigorously test nSPF's functionality.
#
#	Benjamin Proudfoot
#	03/29/22

import numpy as np
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import commentjson as json
from tqdm import tqdm
from getpsf import *
from analysis import *
import tinytim_psfs.make_psf
import os.path
import shutil
import datetime
from astropy.io import fits
from insertpsf import *

# Enable runprops
class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

# These are the grid locations which we want to test nPSF. Subject to change.
# If you don't want to loop over a certain parameter, only leave one value in the
# following arrays.
seps = [1,2,4,8,16]
xys = ["x","y","xy"]
hrats = [1,5,10,100]
#focuses = [-4,-2,0,2,4]
focuses = [1]

h1 = 10000	# the height of the first psf... could make this a multiple of the sky counts???

testbatch = "bright_sep_test"
if not os.path.exists(testbatch):
	os.makedirs(testbatch)

# Load in synthprops.txt which is used to make synthetic images.
runprops = ReadJson("synthprops.txt").outProps()

# Make the TinyTim psfs required. Uses the setting within synthprops (except for focus!)
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
for i, focus in enumerate(focuses):
	filename = "modelpsfs/wfc3psf_" + str(ndet) + "_" + str(nchip) + "_" + filter + "_" + str(xpos) + "_" + str(ypos) + "_" + str(round(focus,1)) + "_" + str(size_psf) + "_" + str(sample_factor) + ".fits"
	if round(focus,1) == 0.0:
		focus = 0.001
	if not os.path.exists(filename):
		tinytim_psfs.make_psf.make_subsampled_model_psf(filename, psf_position = (xpos, ypos), focus = focus, chip = nchip, detector = ndet,
						   filter_name = filter, psf_size = size_psf, subsampling_factor = sample_factor)
	else:
		pass	# File exists... skip this iteration

# Putting PSFs into array
testpsf = getpsf_hst(filename)
numpsfs = len(focuses)
psfs = np.empty((testpsf.shape[0],testpsf.shape[1],numpsfs))
for i, focus in enumerate(focuses):
	filename = "modelpsfs/wfc3psf_" + str(ndet) + "_" + str(nchip) + "_" + filter + "_" + str(xpos) + "_" + str(ypos) + "_" + str(round(focus,1)) + "_" + str(size_psf) + "_" + str(sample_factor) + ".fits"
	psfs[:,:,i] = getpsf_hst(filename)

# Getting charge diffusion kernel from tiny tim psf header. Only one needs to be sampled.
head = fits.open(filename)[0].header
line1 = [float(x) for x in head["COMMENT"][3].split()]
line2 = [float(x) for x in head["COMMENT"][4].split()]
line3 = [float(x) for x in head["COMMENT"][5].split()]
cd_kernel = np.array([line1,line2,line3])
runprops["cd_kernel"] = cd_kernel

# Pulling test image settings from synthprops
image_size = runprops.get("image_size")
skynoise = runprops.get("skynoise")

# Make random noise background image
background = np.random.poisson(lam = skynoise, size = (image_size,image_size))

# Begin looping over appropriate arrays
for i,sep in enumerate(seps):
	for j,xy in enumerate(xys):
		for k,hrat in enumerate(hrats):
			for l,focus in enumerate(focuses):
				# Retrieve the correct PSF from our psf array
				ttpsf = psfs[:,:,l]

				# Find the center of the image
				xcen = float(image_size)/2
				ycen = float(image_size)/2

				# Calculate x and y positions of the 2 psfs based on xy and sep
				if xy == "x":
					psf_x = [xcen, xcen + sep]
					psf_y = [ycen, ycen]
				elif xy == "y":
					psf_x = [xcen, xcen]
					psf_y = [ycen, ycen + sep]
				elif xy == "xy":
					psf_x = [xcen, xcen + sep/np.sqrt(2)]
					psf_y = [ycen, ycen + sep/np.sqrt(2)]

				# Calculate PSF heights
				psf_h = [h1, h1/hrat]

				# Make synthetic image
				#synthimage = insertpsf_n(image = background, psf = ttpsf, xcens = np.array(psf_x), ycens = np.array(psf_y), heights = np.array(psf_h), psfscale = runprops.get("sample_factor"))
				synthimage = insertpsf_n(image = background, psf = ttpsf, xcens = np.array(psf_x), ycens = np.array(psf_y), heights = np.array(psf_h), psfscale = runprops.get("sample_factor"), runprops = runprops)

				# Save said image
				outputname = "../data/" + testbatch + "/sep" + str(sep) + "_" + xy + "_hrat" + str(hrat) + "_focus" + str(round(focus,1)) + ".fits"
				if os.path.exists(outputname):
					os.remove(outputname)
				hdu = fits.PrimaryHDU(synthimage)
				hdu.writeto(outputname)
				#print(outputname)
				#print(psf_x)
				#print(psf_y)
				#print(psf_h)

shutil.copy("synthprops.txt", "../data/" + testbatch + "/othersettings.txt")








