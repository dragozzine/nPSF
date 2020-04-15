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

import math
import numpy as np
import astropy.modeling.functional_models
import astropy.io.fits
import matplotlib.pyplot as plt

# getpsf_2dgau
# Inputs: 
#  - size of psf to return (usually a square, default=(21,21)
#    = remember that you'll probably want to supersample your PSF
#  - covariance matrix that describes the shape of the 2d Guassian, default=((1,0),(0,1))
# Output: a numpy 2d array with the given size that has a *normalized* 2d Gaussian

def getpsf_2dgau(size = (21,21), cov = np.array([[1,0],[0,1]]), supersample = 5):
	# Creating the gaussian model with center at the center of the image
	gaussian = astropy.modeling.functional_models.Gaussian2D(amplitude = 1, 
	  x_mean = ((supersample*size[0])-1)/2, y_mean = ((supersample*size[1])-1)/2, 
	  cov_matrix = supersample**2 * cov)

	# Evaluating the psf at grid values
	psf = np.zeros((supersample*size[0],supersample*size[1]))
	for i in range(supersample*size[0]):
		for j in range(supersample*size[1]):
			psf[i,j]=gaussian.evaluate(i,j,gaussian.amplitude,gaussian.x_mean,
			gaussian.y_mean,gaussian.x_stddev,gaussian.y_stddev,gaussian.theta)
	return psf/psf.sum()


# getpsf_hst gets the stored tinytim psf from the data directory and loads it in
# watch out for subsampling in the HST psf

def getpsf_hst(filename):
	# Loading in the .fits object
	fitsfile = astropy.io.fits.open(filename)
	
	# Getting data out of object
	psf = fitsfile[0].data

	# Normalizing the psf (remove if TT psfs come normalized) and returning the HST psf
	return psf/psf.sum()

def test_getpsf_2dgau(size = (300,300), cov = np.array([[5,1],[1,5]]), supersample = 1):
	# Getting psf from getpsf_2dgau
	psf = getpsf_2dgau(size = size, cov = cov, supersample = supersample)

	# Check if sum == 1 using math.isclose due to floating point errors
	if not math.isclose(psf.sum(), 1.0, rel_tol = 1e-10):
		print("Error in getpsf_2dgau: Sum of psf is not 1.0")
		print("PSF sum = ",psf.sum())
	# Plot the psf to see what it looks like
	plt.figure()
	plt.imshow(psf, cmap='hot', interpolation='nearest')
	plt.colorbar()
	plt.show()
	#plt.savefig()		Make this have a place to save
	plt.close()

def test_getpsf_hst(filename):
	# Getting the HST psf
	psf = getpsf_hst(filename)

	# Plot the psf to see what it looks like
	plt.figure()
	plt.imshow(psf, cmap='hot', interpolation='nearest')
	plt.colorbar()
	plt.show()
	#plt.savefig()		Make this have a place to save
	plt.close()
	



#test_getpsf_2dgau(size = (20,20), cov = np.array([[5,1.5],[1.5,1.5]]))
#test_getpsf_2dgau(size = (20,20), cov = np.array([[5,1.5],[1.5,1.5]]), supersample = 2)
#test_getpsf_2dgau(size = (40,40), cov = 4*np.array([[5,1.5],[1.5,1.5]]))

#test_getpsf_hst("../data/wfc3psf_248_267_50_F350LP_5_00.fits")
