# likelihood.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020
#
# This script includes various functions for calculating likelihoods
#
#Updates by Jarrod Hansen 
#April 22, 2019


# log_likelihood(parameter array) - returns the log of the prior probability for the given parameters
# LATER: log_likelihood_df(parameters dataframe, priors dataframe,observations dataframe) 
# - returns the model dataframe
# log_probability(parameter array) - returns the log posterior probability that emcee uses
# LATER: draw_prior - draws from the priors dataframe
# 
# test_log_likelihood - tests log_likelihood


# log_likelihood
# Inputs:
# - parameter vector [which parameters are which and go with which priors to be handled LATER]
# Outputs:
# - the (natural) log likelihood for these parameters
# including -np.inf for unallowed priors


# log_probability
# Inputs: 
# - parameter vector [which parameters are which and go with which priors to be handled LATER]
# Outputs:
# - the (natural) log posterior probability for these parameters
# 
# Example from emcee (https://emcee.readthedocs.io/en/stable/tutorials/line/):
#def log_probability(theta, x, y, yerr):
#    lp = log_prior(theta)
#    if not np.isfinite(lp):
#        return -np.inf
#    return lp + log_likelihood(theta, x, y, yerr)


import math
import numpy as np
from insertpsf import *
import scipy.stats


def log_likelihood(parameters, image, psf):
    xsize = image.shape[0]
    ysize = image.shape[1]
    
    x1, y1, h1, x2, y2, h2 = parameters

    xcens = np.array([x1,x2])
    ycens = np.array([y1,y2])
    heights = np.array([h1,h2])

    #These are the parameters we need to be able to calulate the likelihood
    #I looked at Ians code and could not figure out how each of these are embedded in his dataframe
    #We will need to write a code to extract the parameters i need

    # Priors for now... since priors never got completed
    if x1 < 0 or x1 > xsize or x2 < 0 or x2 > xsize:
        return -np.inf
    if y1 < 0 or y1 > ysize or y2 < 0 or y2 > ysize:
        return -np.inf
    if h1 < 0 or h2 < 0:
        return -np.inf

    psfimage = insertpsf_n(image = np.zeros((xsize,ysize)), psf = psf,
				xcens = xcens, ycens = ycens, heights = heights)
 
    residuals=image-psfimage
    loglike=0
    mu=200
    #We will need to change this to something that represents the noise profile better?
#    for i in range(xsize):
#        for j in range(ysize):
#            add = scipy.stats.poisson.logpmf(int(residuals[i,j]), mu)
#            loglike += add
            #print(i,j,add)
            #print(residuals[i,j])

    loglike = scipy.stats.poisson.logpmf(np.rint(residuals),mu).sum()
    return loglike

def log_probability(image,parameters):
    lp = log_prior(parameters)#I did not see anything in this code yet
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(image,parameters)


# test log_likelihood
# give parameters with a known answer and confirm it is correct

