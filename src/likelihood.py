# likelihood.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020
#
# This script includes various functions for calculating likelihoods
#
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
from getpsf import getpsf_2dgau
from scipy.stats import poisson


def log_likelihood(image,parameters):
    xsize=len(image)
    ysize=len(image[1])
    xpos1=parameters[1]
    ypos1=parameters[2]
    height1=parameters[3]
    xpos2=parameters[4]
    ypos2=parameters[5]
    height2=parameters[6]
    xpsfsize=parameters[7] 
    ypsfsize=parameters[8] 
    #Size could be a fixed value or something decied by emcee based on the image size to make sure
    #that your psf is not bigger then the image and make sure its an odd number
    
    #These are the parameters we need to be able to calulate the likelihood
    #I looked at Ians code and could not figure out how each of these are embedded in his dataframe
    #We will need to write a code to extract the parameters i need
    psf1=height1*getpsf_2dgau((xpsfsize,ypsfsize),np.array([[1,0],[0,1]]), 1)
    psf2=height2*getpsf_2dgau((xpsfsize,ypsfsize),np.array([[1,0],[0,1]]), 1)
    x=0
    y=0
    residuals=image
    while x < (xsize):
        while y < (ysize):
            if (x-xpos1+xpsfsize/2>-1 && y-ypos1+ypsfsize/2>-1 && x-xpos1+xpsfsize/2<xpsfsize && y-ypos1+ypsfsize/2<ypsfsize):
                residuals=residuals[x][y]-psf1[x-xpos1+xsize/2+.5][y-ypos1+ysize/2+.5]
            y++
        x++
    
    x=0
    y=0
    while x < (xsize):
        while y < (ysize):
            if (x-xpos2+xpsfsize/2>-1 && y-ypos2+ypsfsize/2>-1 && x-xpos2+xpsfsize/2<xpsfsize && y-ypos2+ypsfsize/2<ypsfsize):
                 residuals=residuals[x][y]-psf2[x-xpos2+xsize/2+.5][y-ypos2+ysize/2+.5]
            y++
        x++
        
        
    loglike=0
    x=0
    y=0
    mu=1 #We will need to change this to something that represents the noise profile better
    while x < (xsize):
        while y < (ysize):
            loglike=loglike+logpmf(residuals[x][y],mu,loc=0)
            y++
        x++
        
    return loglike
    

    


def log_probability(image,parameters):
    from likelihood import log_likelihood 
	lp = log_prior(paramaters)#I did not see anything in this code yet
	if not np.isfinite(lp):
		return 0.5
		#return -np.inf
	return lp + log_likelihood(image,parameters)


# test log_likelihood
# give parameters with a known answer and confirm it is correct

