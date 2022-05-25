#
# multi_run.py
#
# Takes a folder of HST images and runs npsf on them
#
# Benjamin Proudfoot
# 01/24/22
#

import numpy as np
import emcee
import sys
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import corner
import commentjson as json
from tqdm import tqdm
from params import params_to_fitarray
from params import fitarray_to_paramsdf
from init_guess import init_guess
from likelihood import * # Jarrod's code for EnsembleSampler
from guess_test import make_test_df
from getpsf import *
from analysis import *
from clustering import *
import tinytim_psfs.make_psf
from csv import writer
import os.path
import datetime
from astropy.io import fits
# from params.py import df_to_array, npsf_init_guess()
from schwimmbad import MPIPool
import shutil
from multi_run import npsf_run

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data


def make_map(psf1params, resultspath, runprops):
    # Loading in image to be solved and making a small postage stamp version
    #x = runprops.get("stamp_x")
    #y = runprops.get("stamp_y")
    x = 0
    y = 0
    size = runprops.get("stamp_size")

    f = runprops.get('input_image')
    #imageraw = getpsf_hst(f)
    imageraw = getimage_hst(f)

    # Clean cosmic rays from image (maybe this can be removed when guesses are good enough?)
    # This may also be irrelevant if we move to simultaneous pair fitting
    import ccdproc
    cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim"), gain_apply = False)
    image = np.array(cr_cleaned)[x:x+size,y:y+size]

    # Make sure that the CR algorithm hasn't marked the target as a cosmic ray!
    reset = False
    while True:
        crredo = False
        if crmask[x:x+size,y:y+size][int(psf1params[0]),int(psf1params[1])]:
            crredo = True
        if crredo:
            cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim") + 1.5, gain_apply = False)
            if reset:
                print("CR rejection algorithm is flagging your target as a CR. Aborting run.") 
                print("Consider increasing the objlim in runprops.")
                #sys.exit()
                return -1, resultspath
            reset = True
        else:
            break
    image = np.array(cr_cleaned)[x:x+size,y:y+size][x:x+size,y:y+size]

    # Ensure there are no negative pixels. Add constant offset, which will be corrected in model images.
    if np.nanmin(image) < 0:
        image = image - np.floor(np.nanmin(image)) + 1.0

    # Calculating image noise characteristics for priors
    runprops["med_noise"] = np.median(image)
    runprops["std_noise"] = np.std(image)
    #hmin = (( runprops.get("med_noise") + runprops.get("std_noise")*runprops.get("noise_cutoff") )/0.1)
    hmin = runprops.get("med_noise")

    # Retrieve the psf we'll be using
    xpos = runprops.get("xpos")
    ypos = runprops.get("ypos")
    size_psf = runprops.get("psf_size")
    sample_factor = runprops.get("sample_factor")
    nchip = runprops.get("chip")
    ndet = runprops.get("det_int")
    filter = runprops.get("filter")

    focus = psf1params[-1]
    filename = "modelpsfs/wfc3psf_" + str(ndet) + "_" + str(nchip) + "_" + filter + "_" + str(xpos) + "_" + str(ypos) + "_" + str(round(focus,1)) + "_" + str(size_psf) + "_" + str(sample_factor) + ".fits"
    psf = getpsf_hst(filename)

    # Getting charge diffusion kernel from tiny tim psf header
    head = fits.open(filename)[0].header
    line1 = [float(x) for x in head["COMMENT"][3].split()]
    line2 = [float(x) for x in head["COMMENT"][4].split()]
    line3 = [float(x) for x in head["COMMENT"][5].split()]
    cd_kernel = np.array([line1,line2,line3])
    runprops["cd_kernel"] = cd_kernel

    # Set up a sampling grid. Initially, have a dx,dy of 1 pixel. Set up for 100 bins of dh
    dx = 0.5
    nx = int(size/dx)
    dy = 0.5
    ny = int(size/dy)
    dh = psf1params[2] - hmin
    nh = 100

    # x and y grids have the 2nd psf placed at the center of the pixel. Not too accurate, but remember, we're looking
    # for relatively rough upper limits
    xgrid = np.linspace(3.0, size - 3.0, num = int(nx))
    ygrid = np.linspace(3.0, size - 3.0, num = int(ny))
    hgrid = np.logspace(np.log10(hmin), np.log10(psf1params[2]), num = int(nh))

    # Set up arrays holding the likelihood map. These are made to look exactly like the flatchain from emcee
    llhoods = np.empty(int(nx*ny*nh))
    grid = np.empty((int(nx*ny*nh), 7))

    # Begin looping
    index = 0
#    for i in tqdm(range(nx)):
#        for j in tqdm(range(ny)):
    for i in tqdm(range(nx)):
        for j in range(ny):
            psf2loc = np.array([xgrid[i],ygrid[j]])
            h_llhoods = log_likelihood_map(psf1params[0:3], psf2loc, hgrid, image, psf, runprops)
            #h_llhoods = np.ones(nh)
            for k in range(nh):
                llhoods[index] = h_llhoods[k]
                grid[index,0] = psf1params[0]	# x1
                grid[index,1] = psf1params[1]   # y1
                grid[index,2] = psf1params[2]   # h1
                grid[index,3] = xgrid[i]        # x2
                grid[index,4] = ygrid[j]        # y2
                grid[index,5] = hgrid[k]        # h2
                grid[index,6] = psf1params[3]   # focus
                #print(round(100*index/(nx*ny*nh),4),h_llhoods[k])
                index += 1
    return grid, llhoods


# Machinery for running from a folder of images
import os
import glob
cwd = os.getcwd()
print(cwd)

# Make sure you are in a data directory
if "data" in cwd:
    folder = cwd
    images = glob.glob("*.fits")
    os.chdir("../../src")
else:
    print("When running multi_map.py, make sure you are located in the data folder")
    sys.exit()

print(os.getcwd())
# Loop over each image
for image in images:
    fullpath = folder + "/" + image
    print(fullpath)
    runprops = ReadJson("runprops.txt").outProps()
    runprops["input_image"] = fullpath
    runprops["image_name"] = image[:-5]
    runprops["npsfs"] = 1

    # Run nPSF for 1 psf
    #bestfit, resultspath = npsf_run(runprops)
    #np.array([x1, y1, h1, focus]) from nPSF run results
    bestfit = np.array([29.173379990355034,29.421632830659068,11393.147542136518,-3.4615050680606383])
    resultspath = "../results/testmap"

    # Make likelihood map
    grid, llhoods = make_map(bestfit, resultspath, runprops)

    # Make plots
    likelihood_map(grid, llhoods, resultspath, runprops)

