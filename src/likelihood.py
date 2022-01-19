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
from csv import writer

def log_likelihood(parameters, image, psfs, focuses, runprops):
    xsize = image.shape[0]
    ysize = image.shape[1]

    xcens = np.empty(2)
    ycens = np.empty(2)
    heights = np.empty(2)

    # Ensure the correct number of input parameters and enforce priors
    if parameters.size == 4:
        x1, y1, h1, focus = parameters
        xcens = np.array([x1,])
        ycens = np.array([y1,])
        heights = np.array([h1,])
        if x1 < 0 or x1 > xsize:
            return -np.inf
        if y1 < 0 or y1 > ysize:
            return -np.inf
        if h1 < 0:
            return -np.inf
    elif parameters.size == 7:
        x1, y1, h1, x2, y2, h2, focus = parameters
        xcens = np.array([x1,x2])
        ycens = np.array([y1,y2])
        heights = np.array([h1,h2])
        if x1 < 0 or x1 > xsize or x2 < 0 or x2 > xsize:
            #print("x")
            return -np.inf
        if y1 < 0 or y1 > ysize or y2 < 0 or y2 > ysize:
            #print("y")
            return -np.inf
        if h1 < 0 or h2 < 0 or h1 < h2:
            #print("h")
            return -np.inf
        #print(h2*runprops.get("noise_cutoff"),( runprops.get("med_noise") + runprops.get("std_noise") ))
        if h2*runprops.get("noise_cutoff") < ( runprops.get("med_noise") + runprops.get("std_noise") ):
            print("noise cutoff")
            return -np.inf
    elif parameters.size == 10:
        x1, y1, h1, x2, y2, h2, x3, y3, h3, focus = parameters
        xcens = np.array([x1,x2,x3])
        ycens = np.array([y1,y2,y3])
        heights = np.array([h1,h2,h3])
        if x1 < 0 or x1 > xsize or x2 < 0 or x2 > xsize or x3 < 0 or x3 > xsize:
            return -np.inf
        if y1 < 0 or y1 > ysize or y2 < 0 or y2 > ysize or y3 < 0 or y3 > ysize:
            return -np.inf
        if h1 < 0 or h2 < 0 or h1 < h2 or h3 < 0 or h1 < h3:
            return -np.inf
        if h2*runprops.get("noise_cutoff") < ( runprops.get("med_noise") + runprops.get("std_noise") ):
            return -np.inf
        if h3*runprops.get("noise_cutoff") < ( runprops.get("med_noise") + runprops.get("std_noise") ):
            return -np.inf
    else:
        print("Wrong number of input parameters. Acceptable numbers are 4, 7, or 10. You have" + str(parameters.size) + ". Aborting run")
        sys.exit()

    if focus < runprops.get("fmin") or focus > runprops.get("fmax"):
        print("f prior", focus)
        return -np.inf

    # Get the median pixel value in the image to make a "blank" image with the right associated noise
    skycounts = runprops.get("med_noise")

    # Choose appropriate PSF based on the focus value
    rfocus = round(focus,1)
    psfindex = np.where(np.isclose(focuses,rfocus))[0][0]
    psf = psfs[:,:,psfindex]

    #rng = np.random.default_rng(42)

    psfimage = insertpsf_n(image = np.ones((xsize,ysize))*skycounts, psf = psf,
				xcens = xcens, ycens = ycens, heights = heights, psfscale = runprops.get("sample_factor"))

#    plt.figure()
#    plt.imshow(psfimage, cmap = "hot", interpolation = "nearest")
#    plt.show()
#    plt.close()
#    plt.figure()
#    plt.imshow(image, cmap = "hot", interpolation = "nearest")
#    plt.show()
#    plt.close()
#    residuals=image-psfimage
#    plt.figure()
#    plt.imshow(residuals, cmap = "hot", interpolation = "nearest")
#    plt.show()
#    plt.close()
#    plt.figure()
#    plt.imshow(psf, cmap = "hot", interpolation = "nearest")
#    plt.show()
#    plt.close()
    loglike=0
#    mu=200
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
    likearray = scipy.stats.poisson.logpmf(np.rint(psfimage),np.rint(image))
    loglike = np.nansum(likearray)
#    print(np.nanmin(likearray),loglike)
#    plt.figure()
#    plt.imshow(likearray, cmap = "hot", interpolation = "nearest")
#    plt.show()
#    plt.close()
    #np.savetxt("loglike.csv", likearray, delimiter = ",")
    #np.savetxt("image.csv", image, delimiter = ",")
    #sys.exit()

    #numnans = np.isnan(scipy.stats.poisson.logpmf(np.rint(psfimage),np.rint(image))).sum()
    #print(np.where)

#    if runprops.get("best_likelihood") < loglike:
#        runprops["best_likelihood"] = loglike
#        with open("../results/best_likelihoods.csv", "a+", newline = "") as write_obj:
#            csv_writer = writer(write_obj, delimiter = ',')
#            out_list = []
#            out_list.insert(0, parameters)
#            out_list.insert(0, loglike)
#            csv_writer.writerow(out_list)

#    print(parameters, loglike)
    return loglike

def log_probability(image,parameters):
    lp = log_prior(parameters)#I did not see anything in this code yet
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(image,parameters)


# test log_likelihood
# give parameters with a known answer and confirm it is correct

