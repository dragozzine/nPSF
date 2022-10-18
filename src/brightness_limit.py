#	brightness_limit.py
#
#	Intent of this code was to determine the brightness detection limit of an object in a specific image. Code was never finished due to simpleness of just getting the information from SAOImageDS9 itself.
#
#	William Giforos
#	08/2022


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
from likelihood import log_likelihood # Jarrod's code for EnsembleSampler
from guess_test import make_test_df
from getpsf import *
from analysis import *
from clustering import *
import tinytim_psfs.make_psf
from csv import writer
import os
import datetime
from astropy.io import fits
# from params.py import df_to_array, npsf_init_guess()
from schwimmbad import MPIPool
import shutil
import astropy.io.fits

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

cwd = os.getcwd()
print(cwd)    
 
# This code is designed to be run from within an object folder in the data directory.    
    
# Make sure you are in a data directory
if "data" in cwd:
    print("Loading the runprops and starting guess")
else:
    print("When running npsf_run.py, make sure you are located in the object folder in the data directory")
    sys.exit()
    
### Load run_props from JSON file ###
runprops = ReadJson("runprops.txt").outProps()

nwalkers = runprops.get("nwalkers")
nsteps = runprops.get("nsteps")
nburnin = runprops.get("nburnin")
output_filename = runprops.get("output_filename")

runprops["best_likelihood"] = -np.inf

# Create results folder
x = datetime.datetime.now()
date = str(x.strftime("%Y"))+"-"+str(x.strftime("%m"))+"-"+str(x.strftime("%d"))+"_"+str(x.strftime("%H"))+"."+str(x.strftime("%M"))+"."+str(x.strftime("%S"))
resultspath = "../../results/"+runprops.get('image_path')[8:]+date+"_"+runprops.get('input_image')[:-5]+"_"+str(runprops.get("npsfs"))+"psf"
if not os.path.exists(resultspath):
    os.makedirs(resultspath)

shutil.copy("runprops.txt", resultspath + "/runprops.txt")
shutil.copy(runprops.get("starting_guess"), resultspath + "/startguess.csv")

# Get initial guess (with init_guess.py)
# This uses test data from guess_test.py
params_df = pd.read_csv(runprops.get("starting_guess"),sep=',',index_col=0)

p0_df = init_guess(params_df, nwalkers)
paramnames = params_df.columns.tolist()
p0, change_dict = params_to_fitarray(p0_df)

# Get ndim
ndim = np.shape(p0)[1]
print("ndim:",ndim)

# Reroute to src directory and recreate the resultspath to run from src directory
os.chdir("../../src")

resultspath = "../results/"+runprops.get('image_path')[8:]+date+"_"+runprops.get('input_image')[:-5]+"_"+str(runprops.get("npsfs"))+"psf"
runprops["resultspath"] = resultspath

# Loading in image to be solved and making a small postage stamp version
x = runprops.get("stamp_x")
y = runprops.get("stamp_y")
size = runprops.get("stamp_size")

f = runprops.get('image_path') + runprops.get('input_image')
#imageraw, filter, nchip, bunits = getimage_hst(f)
# Loading in the .fits object
fitsfile = astropy.io.fits.open(f)
# Getting data out of object
psf = fitsfile[1].data
imageraw = psf

# Clean cosmic rays from image (maybe this can be removed when guesses are good enough?)
# This may also be irrelevant if we move to simultaneous pair fitting
import ccdproc
if runprops.get('cr_clean')==True:
    print("cr_clean set to True")
    cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim"), gain_apply = False)
    image = np.array(cr_cleaned)[x:x+size,y:y+size]

    plt.figure()
    plt.imshow(image, cmap = "hot", interpolation = "nearest", origin = "lower")
    plt.colorbar()
    plt.scatter(p0[:,1].flatten(), p0[:,0].flatten(), color = "blue", marker = "x", s = 5, alpha = 0.2)
    plt.savefig(resultspath + "/cleanedimage.png")
    plt.figure()
    plt.imshow(crmask[x:x+size,y:y+size], cmap = "hot", interpolation = "nearest", origin = "lower")
    plt.colorbar()
    plt.scatter(p0[:,1].flatten(), p0[:,0].flatten(), color = "blue", marker = "x", s = 5, alpha = 0.2)
    plt.savefig(resultspath + "/crmask.png")
    plt.close("all")


    # Make sure that the CR algorithm hasn't marked the target as a cosmic ray!
    reset = False
    maskcheck = runprops.get("maskcheck")

    while maskcheck:
        crredo = False
        for i in range(nwalkers):
            #print(i,crmask[int(p0[i,0]),int(p0[i,1])])
            if crmask[x:x+size,y:y+size][int(p0[i,0]),int(p0[i,1])]:
                crredo = True
        if crredo:
            cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim") + 1.5, gain_apply = False)

            plt.figure()
            plt.imshow(cr_cleaned[x:x+size,y:y+size], cmap = "hot", interpolation = "nearest", origin = "lower")
            plt.colorbar()
            plt.scatter(p0[:,1].flatten(), p0[:,0].flatten(), color = "blue", marker = "x", s = 5, alpha = 0.2)
            plt.savefig(resultspath + "/cleanedimage.png")
            plt.figure()
            plt.imshow(crmask[x:x+size,y:y+size], cmap = "hot", interpolation = "nearest", origin = "lower")
            plt.colorbar()
            plt.scatter(p0[:,1].flatten(), p0[:,0].flatten(), color = "blue", marker = "x", s = 5, alpha = 0.2)
            plt.savefig(resultspath + "/crmask.png")
            plt.close("all")

            if reset:
               print("CR rejection algorithm is flagging your target as a CR. Aborting run.") 
               print("Consider increasing the objlim in runprops.")
               sys.exit()
            reset = True
        else:
            break
    image = np.array(cr_cleaned)[x:x+size,y:y+size]

# Print out image that hasn't been cleaned
else:
    print("cr_clean set to False")
    image = imageraw[x:x+size,y:y+size]
    
    plt.figure()
    plt.imshow(image, cmap = "hot", interpolation = "nearest", origin = "lower")
    plt.colorbar()
    plt.savefig(resultspath + "/cleanedimage.png")
    plt.close("all")

# Ensure there are no negative pixels. Add constant offset, which will be corrected in model images.
if np.nanmin(image) < 0:
    image = image - np.floor(np.nanmin(image)) + 1.0

# Save the image array in case of the need to reconstruct the plot_best_fit plots
#np.save(resultspath + '/img_arr.npy',image)

# Calculating image noise characteristics for priors
runprops["med_noise"] = np.median(image)
print("med_noise:",runprops["med_noise"])
#sys.exit()
