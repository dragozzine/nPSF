# likelihood.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020
#
# This script includes various functions for calculating likelihoods

import math
import numpy as np
from insertpsf import *
import scipy.stats
import sys

def log_likelihood(parameters, image, psfs, focuses):
    xsize = image.shape[0]
    ysize = image.shape[1]

    xcens = np.array([x1,x2])
    ycens = np.array([y1,y2])
    heights = np.array([h1,h2])


    if parameters.size == 4:
        x1, y1, h1, focus = parameters
        xcens = np.array([x1,])
        ycens = np.array([y1,])
        heights = np.array([h1,])
    elif parameters.size == 7:
        x1, y1, h1, x2, y2, h2, focus = parameters
        xcens = np.array([x1,x2])
        ycens = np.array([y1,y2])
        heights = np.array([h1,h2])
    elif parameters.size == 10:
        x1, y1, h1, x2, y2, h2, x3, y3, h3, focus = parameters
        xcens = np.array([x1,x2,x3])
        ycens = np.array([y1,y2,y3])
        heights = np.array([h1,h2,h3])
    else:
        print("Wrong number of input parameters. Acceptable numbers are 4, 7, or 10. You have" + str(parameters.size) + ". Aborting run")
        sys.exit()


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

    # Get the median pixel value in the image to make a "blank" image with the right associated noise
    noise = np.median(image)

    # Choose appropriate PSF based on the focus value
    rfocus = round(focus,1)
    psfindex = np.where(np.isclose(focuses,rfocus))[0][0]
    psf = psfs[:,:,psfindex]

    psfimage = insertpsf_n(image = np.random.poisson(lam = noise, size = (xsize,ysize)), psf = psf,
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

#    loglike = scipy.stats.poisson.logpmf(np.rint(residuals),mu).sum()
# TODO: FIX LIKELIHOOD
# Above was in nPSF before Winter 2021 Physics 529
# But this seemed totally wrong and mu=200 was arbitrary
# DR added this on 1/25/2021:
    loglike = scipy.stats.poisson.logpmf(np.rint(psfimage),np.rint(image)).sum()

    return loglike

def log_probability(image,parameters):
    lp = log_prior(parameters)#I did not see anything in this code yet
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(image,parameters)


# test log_likelihood
# give parameters with a known answer and confirm it is correct

