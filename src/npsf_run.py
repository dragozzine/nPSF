# npsf_run.py
#
# Rochelle Steele
# April 8, 2020
#
# This code takes information about the run, likelihood calculations and parameters and
# runs emcee to get the posteriors in a chain.
#
# Input: run_prop (JSON dictionary), npsf_likelihood (Jarrod), parameters df/array (Ian)
# Output: Chains of parameters from run

""" Run nPSF
Runs emcee with likelihood calculations and parameters from other .py files (likelihood.py
and params.py respectively) and with properties from a dictionary (not written yet).
Inputs:
    JSON file with properties of run
    log-likelihood calculations (likelihood.py)
    parameters as df/array (params.py)

Outputs:
    chain from the run with posteriors
    makes corner plot for parameters
    WILL save a csv file with flat chain (does not have that functionality yet)
"""

### THIS CODE RUNS WITH TEST DATA ###

from statsmodels.regression.linear_model import yule_walker
import numpy as np
import emcee
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import commentjson
#import tqdm
from params import params_to_fitarray
from params import fitarray_to_paramsdf
from init_guess import init_guess
from likelihood import log_likelihood # Jarrod's code for EnsembleSampler
from guess_test import make_test_df
from getpsf import *
from analysis import *
import tinytim_psfs
# from params.py import df_to_array, npsf_init_guess()
np.random.seed(42)

json_file = "run_props.json"
# npsf_likelihood -- from Jarrod?

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

def burnin(p0, nburnin):
    """ This code determines if running the burnin is necessary (based on
    run_props) and then either runs a burnin and returns the resultant state
    or simply returns the initial guess, p0

    Inputs: p0 (initial guess), nburnin (number of steps to run the burnin)

    Output: the state where the sampler will start, either p0 or the state
      after the burnin, depending on whether or not the burnin was run
    """
    print("Beginning the burn in")
    if nburnin == 0:
        return p0
    else:
        state = sampler.run_mcmc(p0, nburnin, progress = True)
        sampler.reset()
        return state

def make_hist(flatchain, param):
    """ Makes histogram of for one parameter. I will probably not call this
    anymore now that I can make corner plots, but it's still accessible.

    Input: flatchain (2D version of results from sampler), param (the index
      of the parameter to be used.

    Output: Saves an image of the histogram. Returns nothing.

    From: https://emcee.readthedocs.io/en/stable/tutorials/quickstart/
    """
    plt.hist(flatchain[:, param], 100, color="k", histtype="step")
    plt.xlabel(r"$\theta_1$")
    plt.ylabel(r"$p(\theta_1)$")
    plt.gca().set_yticks([]);
    plt.savefig('../results/emcee_hist_test.png')
    return

def make_corner_plot(flatchain):
    """ Makes a corner plot from the data from sampler. I may continue to
     make adjustments to this function.

    Input: flatchain (2D version of data from sampler)

    Output: Saves an image of the corner plot. Returns nothing.
    """
    fig = corner.corner(flatchain, labels = np.array(["x1","y1","h1","x2","y2","h2"]),
                        plot_datapoints = False, color = "blue", fill_contours = True,
                        truths = [49.5, 50.2, 2000, 43.0, 57.3, 800], show_titles = True
                        )
    fig.savefig('../results/emcee_corner_test.png')
    return

def make_walker_plots(chain):
    numparams = chain.shape[2]
    numwalkers = chain.shape[1]
    numgens = chain.shape[0]
    names = np.array(["x1","y1","h1","x2","y2","h2"])
    for i in range(numparams):
        plt.figure()
        for j in range(numwalkers):
            plt.plot(np.reshape(chain[0:numgens,j,i], numgens))
        plt.ylabel(names[i])
        plt.xlabel("Generation")
        plt.savefig("../results/" + names[i] + "_walkers.png")
        plt.close()

### Load run_props from JSON file ###
runprops = ReadJson("runprops.txt").outProps()

nwalkers = runprops.get("nwalkers")
nsteps = runprops.get("nsteps")
nburnin = runprops.get("nburnin")
output_filename = runprops.get("output_filename")

# Get initial guess (with init_guess.py)
# This uses test data from guess_test.py
params_df = make_test_df()
p0_df = init_guess(params_df, nwalkers)
p0, change_dict = params_to_fitarray(p0_df)

print(p0.shape)

# Get ndim
ndim = np.shape(p0)[1]

# Loading in image to be solved
image = np.loadtxt("../data/testimage_2objs.txt")

# Creating Tiny Tim PSFs
fmin = runprops.get("fmin")
fmax = runprops.get("fmax")
fspace = runprops.get("fspace")
focuses = np.arange(fmin, fmax + fspace, fspace)

xpos = runprops.get("xpos")
ypos = runprops.get("ypos")

size_psf = runprops.get("psf_size")

sample_factor = runprops.get("sample_factor")

nchip = runprops.get("chip")
ndet = runprops.get("det_int")

filter = runprops.get("filter")

numpsfs = focuses.size
psfs = np.empty((x,y,numpsfs))

for i,focus in enumerate(focuses):
    filename = "model_psf_" + str(focus) + ".fits"
    tinytim_psfs.make_subsampled_model_psf(filename, psf_position = (xpos, ypos), focus = focus, chip = nchip, detector = ndet,
					   filter_name = filter, psf_size = size_psf, subsampling_factor = sample_factor)
    psfs[:,:,i] = getpsf_hst(filename)

# check for bad walkers / bad initial guesses
reset = 0
for i in range(nwalkers):
    llhood = log_likelihood(p0[i,:], image, psfs, focuses)
    print(i, p0[i,:], llhood)
    while (llhood == -np.Inf):
        if(reset % 500 == 0) and (reset != 0):
            print("ERROR: Initial guesses for walkers may be bad.")
            print("Initial guesses have been reset " + str(reset) + " times")
            abort = input("Abort script? (yes/no) ")
            while (abort != "yes") and (abort != "no"):
                print("Invalid input")
                abort = input("Abort script? (yes/no) ")
            if abort == "yes":
                sys.exit()
        # Resetting parameters
        p0_dfreset = init_guess(params_df, 1)
        p0_reset, change_dict = params_to_fitarray(p0_dfreset)
        p0[i,:] = p0_reset
        llhood = log_likelihood(p0[i,:]*(-1), image, psfs, focuses)
        print(p0_reset, llhood, i, reset)
        reset += 1

# Create sampler object
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood, args = [image, psfs, focuses])
#moves = emcee.moves.StretchMove())

# Burnin (optional, default 0 steps)
state = burnin(p0, nburnin)

# Sample
print("Beginning sample")
sampler.run_mcmc(state, nsteps, progress = True)

# Get chain
chain = sampler.get_chain()
flatchain = sampler.get_chain(flat = True)

print(chain.shape)
print(flatchain.shape)

# Make corner plot
make_corner_plot(flatchain)

make_walker_plots(chain)

#testconvergence_geweke(sampler)

# save chain
# Convert chain to df (Ian)
# Save df as csv file
