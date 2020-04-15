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
                  psfscale = 3, psfheight = 1):
	
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
	print(fullpsf.sum())
	
	# Now downsizing the supersampled array
	psfimage = skimage.measure.block_reduce(fullpsf, (psfscale,psfscale))
	print(psfimage.sum())

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
	imageonepsf = insertpsf_one(image = image, xcen = 10.5)

	# Plots the image
	#plt.figure()
	#plt.imshow(imageonepsf, cmap = "hot", interpolation = "nearest")
	#plt.colorbar()
	#plt.title("Composite image")
	#plt.show()
	#plt.close()


def insertpsf_n(image = np.zeros((100,100)), psfs = np.array([getpsf_2dgau()]), 
                xcens = np.array([49.5]), ycens = np.array([49.5]), heights = np.array([1])):
	# Checking if inputs are all the same size
	if not xcens.size == ycens.size == heights.size == psfs.shape[0]:
		raise ValueError("psfs, xcens, ycens, and heights must be the same length")

	# Looping over parameters to insert psfs
	for i in range(xcens.size):
		image = insertpsf_one(image = image, psf = psfs[i], xcen = xcens[i], 
		                      ycen = ycens[i], psfheight = heights[i])

	# Plotting the composite image (remove when fully functioning)
	#plt.figure()
	#plt.imshow(image, cmap = "hot", interpolation = "nearest")
	#plt.colorbar()
	#plt.title("Multi PSF Composite image")
	#plt.show()
	#plt.close()



# Small test case for insertpsf_n
insertpsf_n(image = np.random.random(size = (100,100))*0.01, psfs = np.array([getpsf_2dgau(),
getpsf_2dgau(),getpsf_2dgau()]), xcens = np.array([49.5,47.0,40.0]), ycens = np.array([49.5,49.23,54.98]),
heights = np.array([1.0,0.3,0.2]))



#test_insertonepsf
# makes a simple image and PSF and inserts it at a particular point 
# makes a plot of the image and saves it
