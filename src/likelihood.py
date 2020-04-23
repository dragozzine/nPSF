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
    #The way this is set up the sizes must be integers meaning you cannot search for subpixel variation
    #in the quality of the fits. I think this would be possible to do but is beond the scope of the current project.
    
    xpsfsize=parameters[7] 
    ypsfsize=parameters[8] 
    #Size could be a fixed value or something decied by emcee based on the image size to make sure
    #that your psf is not bigger then the image and make sure its an odd number
    
    #These are the parameters we need to be able to calulate the likelihood
    #I looked at Ians code and could not figure out how each of these are embedded in his dataframe
    #We will need to write a code to extract the parameters i need

 
    psfimage=insertpsf_one(image = np.zeros((xsize,ysize)), psf = getpsf_2dgau(), xcen = xpos1, ycen = ypos1,
                  psfscale = 5, psfheight = height1)
    psfimage=psfimage+insertpsf_one(psfimage, psf = getpsf_2dgau(), xcen = xpos2, ycen = ypos2,
                  psfscale = 5, psfheight = height2)   
    residuals=image-psfimage
    loglike=0
    x=0
    y=0
    mu=1 
    #We will need to change this to something that represents the noise profile better
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

