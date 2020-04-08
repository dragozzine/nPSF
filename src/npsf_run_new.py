# npsf_run_new.py
#
# Rochelle Steele
# April 8, 2020
#
# Input:
# Output: Chains from run

### THIS CODE IS NOT FUNCTIONING YET ###

import numpy as np
import emcee
import sys
import pandas as pd
# from Jarrod import npsf_likelihood
# from Ian import df_to_array



# Define varaibles?
walkers = 50
 # ndim
 # npsf_likelihood -- from Jarrod?
nburnin = 0
nsteps = 500

# Get initial guesses

# Check for bad walkers

# Decide number of steps

# How parameters are there?

# Create sampler object

sampler = emcee.EnsembleSampler(walkers, ndim, npsf_likelihood, args = )

# Burnin (optional, default 0 steps)

state = sampler.run_mcmc(p0, nburnin)

# Sample
