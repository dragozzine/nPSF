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
Inputs:
    JSON file with properties of run
    log-likelihood calculations (Jarrod)
    parameters as df/array (Ian)

Outputs:
    chain from the run with posteriors

"""
### THIS CODE IS NOT FUNCTIONING YET ###


import numpy as np
import emcee
import sys
import pandas as pd
from likelihood import log_probability # Jarrod's code for EnsembleSampler
# from Ian import df_to_array, npsf_init_guess()


# Define varaibles
# ndim - Ian ?
ndim = 3 # will get from Ian's params.py code
# npsf_likelihood -- from Jarrod?
# vvv in run_props vvv
nburnin = 0 # number of steps of burnin (will be deleted)
nsteps = 5000 # number of steps
walkers = 50
filename = "results.csv"


def test_log_probability():
""" Test code to use for log_probability before Jarrod's code is ready
    NOT YET FUNCTIONAL
    Will output a Gaussian log_probability
"""
    return 0.5


# Get initial guesses parameters, from Ian
guesses, change_dict = npsf_init_guess(run_props, nwalkers) # Is this from Ian? What is the output of this function?

guesses = np.zeros((

# Place in numpy array
p0 = np.zeros((ndim, walkers))
for i in range(ndim):
    p0[i,:] = guesses[i+1,:] # why is it i+1?

# Check for bad walkers if you want

# Create sampler object
sampler = emcee.EnsembleSampler(walkers, ndim, log_probability, args = )

# Burnin (optional, default 0 steps)
state = sampler.run_mcmc(p0, nburnin)
sampler.reset()

# Sample
sampler.run_mcmc(state, nsteps)

# Get chain
chain = sampler.get_chain()
flatchain = sampler.get_chain(flat = True)

# save chain (how?)
# Convert chain to df (Ian)
# Save df as csv file

