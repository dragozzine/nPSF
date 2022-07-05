#
# remake_plots.py
#
# Remakes plots using the chain.h5 file
#
# Benjamin Proudfoot
# 01/24/22
#

import glob
import os
import sys
import commentjson as json
import emcee
import numpy as np

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

# Ensure that its being run from the right place
if 'results' in os.getcwd():
    runprops = ReadJson('runprops.txt').outProps()
else:
    print("remake_plots needs to be run from a results directory.")
    sys.exit()
    
# Load in chain file and image/psfs arrays
resultspath = os.getcwd()
backend = emcee.backends.HDFBackend(resultspath + '/chain.h5')
#image = np.load('img_arr.npy')
#psfs = np.load('psfs_arr.npy')

# Calculating image noise characteristics for priors in plot_best_fit
#runprops["med_noise"] = np.median(image)
#runprops["std_noise"] = np.sqrt(runprops["med_noise"])
#print(runprops["std_noise"])
#print((runprops.get("std_noise")*runprops.get("noise_cutoff"))/0.1)

#fmin = runprops.get("fmin")
#fmax = runprops.get("fmax")
#fspace = runprops.get("fspace")
#focuses = np.arange(fmin, fmax + fspace, fspace)

#os.chdir("../../src")
os.chdir("../../../src")

# Rerun plots
from analysis import *
#plot_best_fit(backend, image, psfs, focuses, runprops)
plots(backend, resultspath, runprops)
