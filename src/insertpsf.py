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
                  psfscale = 1, psfheight = 1):
	
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
	
	plt.figure()
	plt.imshow(psfimage, cmap = "hot", interpolation = "nearest")
	plt.colorbar()
	plt.show()




insertpsf_one(xcen = 70.5, ycen = 30.5, psfscale = 3)
	
#test_insertonepsf
# makes a simple image and PSF and inserts it at a particular point 
# makes a plot of the image and saves it
