# npsf_run_new.py
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

import numpy as np
import emcee
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
# from likelihood import log_probability # Jarrod's code for EnsembleSampler
# from params.py import df_to_array, npsf_init_guess()
np.random.seed(42)


ndim = 5 # will get from Ian's params.py code
# npsf_likelihood -- from Jarrod?
# vvv in run_props vvv
nburnin = 0 # number of steps of burnin (will be deleted)
nsteps = 5000 # number of steps
nwalkers = 32 # number of walkers
filename = "results.csv" # name of output file containing  flatchain


def log_probability(x, mu, cov):
    """ Test code for the log probaiblity before Jarrod's code (likelihood) is ready
    This put gives the log probability of a Gaussian to be used in EnsembleSampler
    Input: x (p0), mu (mean), cov (Sigma)
    Output: log probability of Gaussian
    From: https://emcee.readthedocs.io/en/stable/tutorials/quickstart/
    """
    diff = x - mu
    return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))

def means_cov():
    """ More test code. This is needed for log_probability (test version).
    This code calculated the mean and cov of a Gaussian that will be created
    in log_probability.
    Input: [none]
    Output: means, cov 
    From: https://emcee.readthedocs.io/en/stable/tutorials/quickstart/
    """
    means = np.random.rand(ndim)

    cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
    cov = np.triu(cov)
    cov += cov.T - np.diag(cov.diagonal())
    cov = np.dot(cov, cov)

    return means, cov

def burnin(p0, nburnin):
    """ This code determines if running the burnin is necessary (based on
    run_props) and then either runs a burnin and returns the resultant state
    or simply returns the initial guess, p0
    Inputs: p0 (initial guess), nburnin (number of steps to run the burnin)
    Output: the state where the sampler will start, either p0 or the state
      after the burnin, depending on whether or not the burnin was run
    """
    if nburnin == 0:
        return p0
    else:
        state = sampler.run_mcmc(p0, nburnin)
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
    plt.savefig('../results/emcee_Gauss_hist.png')
    return

def make_corner_plot(flatchain):
    """ Makes a corner plot from the data from sampler. I may continue to
     make adjustments to this function.
    Input: flatchain (2D version of data from sampler)
    Output: Saves an image of the corner plot. Returns nothing.
    """
    fig = corner.corner(flatchain)#, labels=labels, truths=[m_true, b_true, np.log(f_true)])
    fig.savefig('../results/emcee_corner_Gaussian.png')
    return

#######################################
### For use when params.py is ready ###

## Get initial guesses parameters, from Ian
#guesses, change_dict = npsf_init_guess(run_props, nwalkers) # Is this from Ian? What is the output of this function?

## Place in numpy array
#p0 = np.zeros((ndim, walkers))
#for i in range(ndim):
#    p0[i,:] = guesses[i+1,:] # why is it i+1?

########################################

# Get initial guesses (for test case)
p0 = np.random.rand(nwalkers, ndim) # From: https://emcee.readthedocs.io/en/stable/tutorials/quickstart/

# Get means and cov (for test code)
means, cov = means_cov()

###########################################################
### This is where the code could check for bad walkers. ### 
### I have not written this functionality yet.          ###
###########################################################

# Create sampler object
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args = [means, cov] )

# Burnin (optional, default 0 steps)
state = burnin(p0, nburnin)

# Sample
sampler.run_mcmc(state, nsteps)

# Get chain
chain = sampler.get_chain()
flatchain = sampler.get_chain(flat = True)

# Make corner plot
make_corner_plot(flatchain)

# save chain
# Convert chain to df (Ian)
# Save df as csv file
