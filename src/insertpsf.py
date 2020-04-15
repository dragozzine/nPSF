# insertpsf.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020
#
# This script includes various functions for inserting PSFs into images
#
# insertonepsf(image, psf, xcen, ycen, psfscale, psfheight) - inserts one PSF
#   taken from the psf 2-d numpy array into the image at the position xcen, ycen
#   scaled in size by psfscale and multiplied by psfheight
# LATER: quickinsertonepsf - like insertonepsf but xcen, ycen, and psfscale are all integers
# LATER?: insertpsf_insertNpsfs() - loops over insertonepsf
# 
# test_insertonepsf - tests insertonepsf
# 
# 

import numpy as np
from getpsf import *
import scipy.ndimage
import skimage.measure

# insertonepsf
# Inputs:
# - 2-dimensional numpy array "image" where PSF will be inserted, default = np.zeros(100,100)
# - 2-dimensional numpy array "psf" which will be inserted, default = getpsf_2dgau()
# - xcen: position in the x-direction of the center of the PSF, default=center of image
# - ycen: position in the y-direction of the center of the PSF,  default=center of image
# - psfscale: the super/sub-sampling scale to be used for the PSF, default=1
# - psfheight: the height of the PSF, e.g., the PSF is multipled by this number when inserted, default=1
# Output: 2-d numpy array of the image with the added PSF

# Note, some of the defaults above may not be able to be implemented as actual "defaults"
# in the function call, but if they are "None" you could use these values

def insertpsf_one(image = np.zeros((100,100)), psf = getpsf_2dgau(), xcen = 49.5, ycen = 49.5,
                  psfscale = 5, psfheight = 1):
	
	# Rescaling the psf
	psf = psfheight*psf
	
	# Finding the shift to put center of psf at (xcen, ycen)
	psfxcen = (psf.shape[0] - 1)/2
	psfycen = (psf.shape[1] - 1)/2
	supxcen = xcen*psfscale + (psfscale - 1)*0.5	
	supycen = ycen*psfscale + (psfscale - 1)*0.5
	
	xshift = supxcen - psfxcen
	yshift = supycen - psfycen
	
	# Putting the psf into a full supersampled array
	fullpsf = np.zeros((image.shape[0]*psfscale,image.shape[1]*psfscale))
	fullpsf[0:psf.shape[0], 0:psf.shape[1]] = psf
	
	# Shifting the full supersampled psf to the correct center
	fullpsf = scipy.ndimage.shift(fullpsf, (xshift,yshift))
	
	# Now downsizing the supersampled array
	psfimage = skimage.measure.block_reduce(fullpsf, (psfscale,psfscale))

	# Plotting the PSF image to check (remove this once fully functioning)	
	#plt.figure()
	#plt.imshow(psfimage, cmap = "hot", interpolation = "nearest")
	#plt.colorbar()
	#plt.title("PSF image")
	#plt.show()
	#plt.close()
	
	# Returns passed in image + psf image
	return (image + psfimage)


def test_insertpsf_one(image = np.random.random(size = (100,100))*0.01):
	# Creates composite image of image passed in + psf image
	imageonepsf = insertpsf_one(image = image, xcen = 45.5, psf = getpsf_hst("../data/wfc3psf_248_267_50_F350LP_5_00.fits"))

	# Plots the image
	plt.figure()
	plt.imshow(imageonepsf, cmap = "hot", interpolation = "nearest")
	plt.colorbar()
	plt.title("Composite image")
	plt.show()
	plt.close()


def insertpsf_n(image = np.zeros((100,100)), psf = getpsf_hst("../data/wfc3psf_248_267_50_F350LP_5_00.fits"), 
                xcens = np.array([49.5]), ycens = np.array([49.5]), heights = np.array([1])):
	# Checking if inputs are all the same size
	if not xcens.size == ycens.size == heights.size:
		raise ValueError("xcens, ycens, and heights must be the same length")

	# Looping over parameters to insert psfs
	for i in range(xcens.size):
		image = insertpsf_one(image = image, psf = psf, xcen = xcens[i], 
		                      ycen = ycens[i], psfheight = heights[i])

	# Returns image after convolving with the cd kernel
	return cd_convolve(image)

def test_insertpsf_n(image = np.random.random(size = (100,100))*0.01):
	# Creates composite image of image passed in + psf image
	psf = getpsf_hst("../data/wfc3psf_248_267_50_F350LP_5_00.fits")
	xcen_arr = np.array([45.0,30.45678,47.10])
	ycen_arr = np.array([43.7,60.0,43.5])
	height_arr = np.array([1.0,0.234,0.1])
	image_npsfs = insertpsf_n(image = image, psf = psf, xcens = xcen_arr,
				  ycens = ycen_arr, heights = height_arr)

	# Plots the image
	plt.figure()
	plt.imshow(image_npsfs, cmap = "hot", interpolation = "nearest")
	plt.colorbar()
	plt.title("Composite image, w/ CD")
	plt.show()
	plt.close()

def cd_convolve(image):
	# Defining the two charge diffusion kernels for ir and uv. For now I will use ir,
	# but there could be a way of dealing with these to get a better value
	cd_kernel_ir = np.array([ [0.002, 0.037, 0.002],
				  [0.037, 0.844, 0.037],
				  [0.002, 0.037, 0.002] ])

	cd_kernel_uv = np.array([ [0.027, 0.111, 0.027],
				  [0.111, 0.432, 0.111],
				  [0.027, 0.111, 0.027] ])

	# Convolving the image with the cd kernel
	return scipy.ndimage.convolve( image, cd_kernel_ir, mode='constant', cval=0.0 )

def make_image(paramsdf):
	blank = np.zeros((100,100))

	# Getting all of the inputs for insertpsf_n ready
	# If needed getpsf_hst can be replaced with getpsf_2dgau to fit have gaussian psfs instead
	# In the future, make this read the run_props dictionary to know which function to use
	psf = getpsf_hst("../data/wfc3psf_248_267_50_F350LP_5_00.fits")
	xcen_arr = paramsdf["xcen"].values
	ycen_arr = paramsdf["ycen"].values
	height_arr = paramsdf["heights"].values
	insertpsf_n(blank, psf = psf, xcens = xcen_arr, ycens = ycen_arr, heights = height_arr)

