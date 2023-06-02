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
from priors import log_prior
import random

def log_likelihood(parameters, image, psfs, focuses, runprops, plotit = False):
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
	# Noise cutoff. Central pixel of PSF has 10% of total flux. Reject values of h2 where the central pixel
        # of the second PSF would be below the one sigma upper bound on the image's noise. If h2 can be below this
        # it begins fitting to the image's noise, which really borks the run.
        if 0.1*h2 < ( runprops.get("std_noise")*runprops.get("noise_cutoff") ):
            #print("noise cutoff")
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
        if 0.1*h2 < ( runprops.get("std_noise")*runprops.get("noise_cutoff") ):
            return -np.inf
        if 0.1*h3 < ( runprops.get("std_noise")*runprops.get("noise_cutoff") ):
            return -np.inf
    else:
        print("Wrong number of input parameters. Acceptable numbers are 4, 7, or 10. You have" + str(parameters.size) + ". Aborting run")
        sys.exit()

    if focus < runprops.get("fmin") or focus > runprops.get("fmax"):
        #print("f prior", focus)
        return -np.inf

    # Get the median pixel value in the image to make a "blank" image with the right associated noise
    skycounts = runprops.get("med_noise")

    # Choose appropriate PSF based on the focus value
    if focuses.all() != 0:
        rfocus = round(focus,1)
        psfindex = np.where(np.isclose(focuses,rfocus))[0][0]
        psf = psfs[:,:,psfindex]
        
        psfimage = insertpsf_n(image = np.ones((xsize,ysize))*skycounts, psf = psf,
				xcens = xcens, ycens = ycens, heights = heights,
				psfscale = runprops.get("sample_factor"), runprops = runprops)
    else:
        # upload the single input psf (empirical psf), no focus values used
        psf = psfs[:,:,0]
        # insert the psf without convolving the charge diffusion kernel by setting runprops = None
        psfimage = insertpsf_n(image = np.ones((xsize,ysize))*skycounts, psf = psf,
				xcens = xcens, ycens = ycens, heights = heights,
				psfscale = runprops.get("sample_factor"), runprops = None)

    #rng = np.random.default_rng(42)

    loglike=0
    likearray = scipy.stats.poisson.logpmf(np.rint(psfimage),np.rint(image))
    loglike = np.nansum(likearray)

    # Make plots if plotit flag is true
    if plotit:
        # Model image
        plt.figure()
        plt.imshow(psfimage, cmap = "hot", interpolation = "nearest", origin = "lower")
        plt.colorbar()
        plt.savefig(runprops.get("resultspath") + "/bestfit.png", dpi = 300)
        # Residual image
        residuals=image-psfimage
        plt.figure()
        plt.imshow(residuals, cmap = "hot", interpolation = "nearest", origin = "lower")
        plt.colorbar()
        plt.scatter(ycens, xcens, color = "blue", marker = "o")
        plt.xlim(ycens[0]-10, ycens[0]+10)
        plt.ylim(xcens[0]-10, xcens[0]+10)
        plt.savefig(runprops.get("resultspath") + "/residuals.png", dpi = 300)
        # Likelihood image
        plt.figure()
        plt.imshow(likearray, cmap = "hot", interpolation = "nearest", origin = "lower")
        plt.colorbar()
        plt.savefig(runprops.get("resultspath") + "/llhood.png", dpi = 300)
        plt.close("all")
        return

#    if runprops.get("best_likelihood") < loglike:
#        runprops["best_likelihood"] = loglike
#        with open("../results/best_likelihoods.csv", "a+", newline = "") as write_obj:
#            csv_writer = writer(write_obj, delimiter = ',')
#            out_list = []
#            out_list.insert(0, parameters)
#            out_list.insert(0, loglike)
#            csv_writer.writerow(out_list)

    #print(parameters, loglike)
    return loglike

def generate_bestfit_residual(parameters, image, psfs, focuses, runprops, plotit = False):
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
	# Noise cutoff. Central pixel of PSF has 10% of total flux. Reject values of h2 where the central pixel
        # of the second PSF would be below the one sigma upper bound on the image's noise. If h2 can be below this
        # it begins fitting to the image's noise, which really borks the run.
        if 0.1*h2 < ( runprops.get("std_noise")*runprops.get("noise_cutoff") ):
            #print("noise cutoff")
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
        if 0.1*h2 < ( runprops.get("std_noise")*runprops.get("noise_cutoff") ):
            return -np.inf
        if 0.1*h3 < ( runprops.get("std_noise")*runprops.get("noise_cutoff") ):
            return -np.inf
    else:
        print("Wrong number of input parameters. Acceptable numbers are 4, 7, or 10. You have" + str(parameters.size) + ". Aborting run")
        sys.exit()
    
    if focus < runprops.get("fmin") or focus > runprops.get("fmax"):
        #print("f prior", focus)
        return -np.inf

    # Get the median pixel value in the image to make a "blank" image with the right associated noise
    skycounts = runprops.get("med_noise")

    # Choose appropriate PSF based on the focus value
    if focuses.all() != 0:
        rfocus = round(focus,1)
        psfindex = np.where(np.isclose(focuses,rfocus))[0][0]
        psf = psfs[:,:,psfindex]

        psfimage = insertpsf_n(image = np.ones((xsize,ysize))*skycounts, psf = psf,
				xcens = xcens, ycens = ycens, heights = heights,
				psfscale = runprops.get("sample_factor"), runprops = runprops)
    else:
        # upload the single input psf (empirical psf)
        psf = psfs[:,:,0]
        # insert the psf without convolving the charge diffusion kernel
        psfimage = insertpsf_n(image = np.ones((xsize,ysize))*skycounts, psf = psf,
				xcens = xcens, ycens = ycens, heights = heights,
				psfscale = runprops.get("sample_factor"), runprops = None)

    residuals = image - psfimage
    
    return psfimage, residuals, xcens, ycens

def log_likelihood_map(psf1params, psf2loc, psf2heights, image, psf, runprops):
    xsize = image.shape[0]
    ysize = image.shape[1]
    image_int = np.rint(image)

    # Get the median pixel value in the image to make a "blank" image with the right associated noise
    skycounts = runprops.get("med_noise")

    # Make initial image with a single psf within it
    psf1image = insertpsf_one(image = np.ones((xsize,ysize))*skycounts, psf = psf, xcen = psf1params[0],
			      ycen = psf1params[1], psfheight = psf1params[2], psfscale = runprops.get("sample_factor"))
    psf2image_norm = insertpsf_one(image = np.zeros((xsize,ysize)), psf = psf, xcen = psf2loc[0],
                              ycen = psf2loc[1], psfheight = 1.0, psfscale = runprops.get("sample_factor"))

    # Prep arrrays for storage of likelihood values
    llhoods = np.empty(psf2heights.size)

    # Now loop over the heights of the 2nd psf
    # add psf images together (with correct heights) and compare to input image
    for i in range(psf2heights.size):
        modelimage_nocd = psf1image + psf2image_norm*psf2heights[i]
        modelimage = cd_convolve(modelimage_nocd, runprops)

        likearray = scipy.stats.poisson.logpmf(np.rint(modelimage), image_int)
        llhoods[i] = np.nansum(likearray)

    return llhoods


def log_probability(parameters, image, psfs, focuses, runprops):
    llhood = log_likelihood(parameters, image, psfs, focuses, runprops)
    lp = log_prior(parameters) #log_prior is set to return 0, as there are no priors for nPSF.
    if not np.isfinite(lp):
        return -np.inf
    return lp + llhood

            
# test log_likelihood
# give parameters with a known answer and confirm it is correct

