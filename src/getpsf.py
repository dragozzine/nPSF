# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020 
# 
# Outline completely filled in and all functions tested
# Benjamin Proudfoot 04/20/20



import math
import numpy as np
import astropy.modeling.functional_models
import astropy.io.fits
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

def getpsf_2dgau(size = (21,21), cov = np.array([[1,0],[0,1]]), supersample = 5):
	"""
	Creates a psf with a gaussian profile centered at the center of the output
	array/image

	Inputs:
		size: tuple of size 2 with each entry being the length of an axis

		cov: a covariance matrix specifying the shape of the gaussian profile

		supersample: a ratio at which to super sample the psf
				for example: specifying a 21x21 image with a super sample
					     ratio of 5, will return a 105x105 image

	Outputs:
		array/image of a psf with a gaussian profile

	"""

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



def getpsf_hst(filename):
	"""
	Creates an array/image from a specifed TinyTim file.

	Inputs:
		filename: a string specifying the location of the TinyTim psf

	Outputs:
		an array/image of an HST psf
	"""
	# Loading in the .fits object
	fitsfile = astropy.io.fits.open(filename)

	# Getting data out of object
	psf = fitsfile[0].data

	# Normalizing the psf (remove if TT psfs come normalized) and returning the HST psf
	return psf


def getimage_hst(filename):
	"""
	Creates an array/image from a specifed fits file.

	Inputs:
		filename: a string specifying the location of the image

	Outputs:
		an array/image of anHST image
	"""
	# Loading in the .fits object
	fitsfile = astropy.io.fits.open(filename)

	# Getting data out of object
	psf = fitsfile[1].data

	# Normalizing the psf (remove if TT psfs come normalized) and returning the HST psf
	return psf




def test_getpsf_2dgau(size = (10,10), cov = np.array([[5,1],[1,5]]), supersample = 1):
	"""
	Runs a simple test for getpsf_2dgau and plots an image of the psf.

	Inputs:
		size: tuple of size 2 specifying the lengths of each axis in the final image
	
		cov: a covariance matrix specifying the shape of the 2d gaussian

		supersample: a ratio at which to supersample the psf
	
	Outputs:
		none
	"""

	# Getting psf from getpsf_2dgau
	psf = getpsf_2dgau(size = size, cov = cov, supersample = supersample)

	# Check if sum == 1 using math.isclose due to floating point errors
	if not math.isclose(psf.sum(), 1.0, rel_tol = 1e-10):
		print("Error in getpsf_2dgau: Sum of psf is not 1.0")
		print("PSF sum = ",psf.sum())
	# Plot the psf to see what it looks like
	plt.figure()
	plt.imshow(psf, cmap='hot', interpolation='nearest')
	#plt.colorbar()
	plt.show()
	#plt.savefig()		Make this have a place to save
	plt.close()

def test_getpsf_hst(filename):
	"""
	Runs a simple test for getpsf_hst and plots an image of the psf.

	Inputs:
		filename: a string specifying the location of the TinyTim psf.
	
	Outputs:
		none
	"""

	# Getting the HST psf
	psf = getpsf_hst(filename)

	# Plot the psf to see what it looks like
	plt.figure()
	plt.imshow(psf, cmap='hot', interpolation='nearest')
	plt.colorbar()
	plt.show()
	#plt.savefig()		Make this have a place to save
	plt.close()
	

