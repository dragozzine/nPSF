# getpsf.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020 
#
# This script includes various functions for getting (or making) PSFs.
# 
# getpsf_2dgau(shape, covariance) - returns a centered 2D Guassian with the given shape and covariance
# LATER: getpsf_hst
# 
# It also includes some test functions
# test_getpsf_2dgau

import numpy as np

# getpsf_2dgau
# Inputs: 
#  - size of psf to return (usually a square, default=(21,21)
#    = remember that you'll probably want to supersample your PSF
#  - covariance matrix that describes the shape of the 2d Guassian, default=((1,0),(0,1))
# Output: a numpy 2d array with the given size that has a *normalized* 2d Gaussian

def getpsf_2dgau(size = 21):
	psf = np.ones(size,size)
	return psf

# Notes: 
# Watch out for how to do the centering when the size in one direction is odd vs. even

# To handle the covariance matrix, see the following
# https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Gaussian2D.html
# https://www.unige.ch/sciences/astro/files/5413/8971/4090/2_Segransan_StatClassUnige.pdf Slide 18
# https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
# and feel free to use the first one for calculations

# Before returning, normalize the PSF by dividing the whole 2d array by the total of the whole array




# test_getpsf_2dgau
# Runs a test of getpsf_2dgau
# inputs values big enough to make an image (300,300) of a rotaated 2d Guassian
# checks that the sum of the returned image is 1
# plots the image and saves to a file
