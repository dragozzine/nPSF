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
import tinytim_psfs.make_psf
from csv import writer
# from params.py import df_to_array, npsf_init_guess()
#np.random.seed(42)

from schwimmbad import MultiPool as Pool

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

def make_corner_plot(flatchain, names, filename):
    """ Makes a corner plot from the data from sampler. I may continue to
     make adjustments to this function.

    Input: flatchain (2D version of data from sampler)

    Output: Saves an image of the corner plot. Returns nothing.
    """
    fig = corner.corner(flatchain, labels = names,
                        plot_datapoints = False, color = "blue", fill_contours = True, show_titles = True,
                        bins = 40)
    fig.savefig('../results/' + filename)
    return

def make_walker_plots(chain):
    numparams = chain.shape[2]
    numwalkers = chain.shape[1]
    numgens = chain.shape[0]
    names = np.array(["x1","y1","h1","x2","y2","h2","f"])
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

runprops["best_likelihood"] = -np.inf

# Get initial guess (with init_guess.py)
# This uses test data from guess_test.py
params_df = pd.read_csv(runprops.get("starting_guess"),sep=',',index_col=0)

print(params_df)

p0_df = init_guess(params_df, nwalkers)
p0, change_dict = params_to_fitarray(p0_df)

print(p0.shape)

# Get ndim
ndim = np.shape(p0)[1]

# Loading in image to be solved
f = runprops.get('input_image')
imageraw = getimage_hst(f)[200:300,200:300]

# Clean cosmic rays from image (maybe this can be removed when guesses are good enough?)
import ccdproc
cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5, gain_apply = False)
image = np.array(cr_cleaned)

# Getting inputs for Tiny Tim PSFs
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

filename = "model_psf_1.5.fits"

# Making TinyTim PSFs
for i,focus in enumerate(focuses):
    filename = "model_psf_" + str(round(focus,1)) + ".fits"
    tinytim_psfs.make_psf.make_subsampled_model_psf(filename, psf_position = (xpos, ypos), focus = focus, chip = nchip, detector = ndet,
					   filter_name = filter, psf_size = size_psf, subsampling_factor = sample_factor)

# Load in one PSF to check the size and make array for PSFs
testpsf = getpsf_hst(filename)
psfs = np.empty((testpsf.shape[0],testpsf.shape[1],numpsfs))

# Putting PSFs into array
for i, focus in enumerate(focuses):
    filename = "model_psf_" + str(round(focus,1)) + ".fits"
    psfs[:,:,i] = getpsf_hst(filename)

# check for bad walkers / bad initial guesses
reset = 0
for i in range(nwalkers):
    llhood = log_likelihood(p0[i,:], image, psfs, focuses, runprops)
    print(i, p0[i,:], llhood)
    while (llhood == -np.Inf):
        # Resetting parameters
        p = random.random()
        p0[i,:] = p*p0[random.randrange(nwalkers),:] + (1-p)*p0[random.randrange(nwalkers),:]
        llhood = log_likelihood(p0[i,:], image, psfs, focuses, runprops)
        #print(reset, llhood, i, reset)
        reset += 1
        if reset > 10000:
            sys.exit()

# Create sampler object
with Pool(4) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood, args = [image, psfs, focuses, runprops], pool = pool, moves = emcee.moves.StretchMove(a = 2.0))

    # Burnin (optional, default 0 steps)
    state = burnin(p0, nburnin)

    # Sample
    print("Beginning sample")
    sampler.run_mcmc(state, nsteps, progress = True)

# Get chain
chain = sampler.get_chain()
flatchain = sampler.get_chain(flat = True)
llhoods = sampler.get_log_prob(flat = True)

print(chain.shape)
print(flatchain.shape)

# Convert x,y pixel coordinates to dRA and ddec
#get wcs
import astropy.io
import astropy.wcs
f = astropy.io.fits.open("../data/idsj03wgq_flc.fits")
w = astropy.wcs.WCS(f[2].header)
ra1,dec1 = w.pixel_to_world_values(flatchain[:,0].flatten(), flatchain[:,1].flatten())
ra2,dec2 = w.pixel_to_world_values(flatchain[:,3].flatten(), flatchain[:,4].flatten())

dra = (ra2 - ra1)*3600*np.cos(np.deg2rad(dec1))
ddec = (dec2 - dec1)*3600

# Make corner plot
names = np.array(["x1","y1","h1","x2","y2","h2","f"])
make_corner_plot(flatchain, names, "corner.pdf")

dnames = names.copy()
dfchain = flatchain.copy()

dnames = np.append(dnames, ["dra","ddec"])
dfchain = np.concatenate((dfchain,np.array(dra).reshape((dra.size,1)) ), axis = 1)
dfchain = np.concatenate((dfchain,np.array(ddec).reshape((ddec.size,1)) ), axis = 1)

make_corner_plot(dfchain, dnames, "cornerderived.pdf")

make_walker_plots(chain)

likelihoodspdf = PdfPages("../results/likelihoods.pdf")
ylimmin = np.percentile(llhoods.flatten(), 1)
ylimmax = llhoods.flatten().max() + 1
names = np.array(["x1","y1","h1","x2","y2","h2","f"])
for i in range(dnames.size):
    plt.figure(figsize = (9,9))
    plt.subplot(221)
    plt.hist(dfchain[:,i].flatten(), bins = 40, histtype = "step", color = "black")
    plt.subplot(223)
    plt.scatter(dfchain[:,i].flatten(), llhoods.flatten(),
                c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), nwalkers),
                cmap = "nipy_spectral", edgecolors = "none", rasterized=True, alpha=0.1)
    plt.xlabel(dnames[i])
    plt.ylabel("Log(L)")
    plt.ylim(ylimmin, ylimmax)
    plt.subplot(224)
    llflat = llhoods.flatten()
    plt.hist(llflat[np.isfinite(llflat)], bins = 40, orientation = "horizontal",
             histtype = "step", color = "black")
    plt.ylim(ylimmin, ylimmax)
    likelihoodspdf.savefig()

likelihoodspdf.close()
plt.close("all")





#testconvergence_geweke(sampler)
