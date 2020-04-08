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
# from Jarrod import npsf_likelihood
# from Ian import df_to_array


# Define varaibles
# ndim - Ian ?
# npsf_likelihood -- from Jarrod?
# vvv in run_props vvv
nburnin = 0 # number of steps of burnin (will be deleted)
nsteps = 5000 # number of steps
walkers = 50
filename = "results.csv"

# Get initial guesses - from who???
# p0

# parameters from Ian

# Check for bad walkers

# Create sampler object
sampler = emcee.EnsembleSampler(walkers, ndim, npsf_likelihood, args = )

# Burnin (optional, default 0 steps)
state = sampler.run_mcmc(p0, nburnin)
sampler.reset()

# Sample
sampler.run_mcmc(state, nburnin)

# Get chain
chain = sampler.get_chain()
flatchain = sampler.get_chain(flat = True)

# save chain (how?)
# Convert chain to df (Ian)
# Save df as csv file

