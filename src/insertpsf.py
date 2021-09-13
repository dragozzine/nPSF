# insertpsf.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020
#
# Outline completely filled in. Added various functions to insert psfs and make images
# Benjamin Proudfoot 04/20/20

import numpy as np
from getpsf import *
import scipy.ndimage
import skimage.measure
import matplotlib.colors as colors


# Note, some of the defaults above may not be able to be implemented as actual "defaults"
# in the function call, but if they are "None" you could use these values

def insertpsf_one(image = np.zeros((100,100)), psf = getpsf_2dgau(), xcen = 49.5, ycen = 49.5,
                  psfscale = 5, psfheight = 1):

	"""
	Takes a psf created by getpsf.py and inserts it into an input image at the location
	specified with the correct supersampling ratios. 

	Inputs:
		image: an ndarray to insert the psf into (default is 100x100 array of zeros)
	
		psf: an ndarray created by a function in get_psf
	
		xcen: float, the x-position to place the psf

		ycen: float, the y-position to place the psf

		psfscale: int, the super sampling ratio to work with
		
		psfheight: float, the total integrated flux of the psf

	Outputs:
		ndarray of the combined image + psf
	"""
	
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

	# Returns passed in image + psf image
	return (image + psfimage)


def test_insertpsf_one(image = np.random.random(size = (100,100))*0.01):
	"""
	Runs a simple test for insertpsf_one and plots the result

	Inputs:
		image: ndarray of image to insert a test psf into (default is
			a 100x100 array of random noise)
	
	Outputs:
		none
	"""

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
                xcens = np.array([49.5]), ycens = np.array([49.5]), heights = np.array([1]), psfscale = 5):
	"""
	A function which provides a wrapper for looping over insertpsf_one.

	Inputs:
		image: ndarray to start with. Default is 100x100 array of zeros.

		psf: ndarray, created using getpsf.py with which to create the image

		xcens: ndarray with size 1x(number of psfs to add) of the x-positions of the psfs

		ycens: ndarray with size 1x(number of psfs to add) of the y-positions of the psfs

		heights: ndarray with size 1x(number of psfs to add) of the heights of the psfs

	Outputs:
		an ndarray with psfs added to the input image. Note that the output image has
		been convolved with the charge diffusion kernel.
	"""

	# Checking if inputs are all the same size
	if not xcens.size == ycens.size == heights.size:
		raise ValueError("xcens, ycens, and heights must be the same length")

	# Looping over parameters to insert psfs
	for i in range(xcens.size):
		image = insertpsf_one(image = image, psf = psf, xcen = xcens[i], 
		                      ycen = ycens[i], psfheight = heights[i], psfscale = psfscale)

	# Returns image after convolving with the cd kernel
	return cd_convolve(image)

def test_insertpsf_n(image = np.random.random(size = (100,100))*10.0):
	"""
	Runs a simple test for insertpsf_n and plots the result.

	Inputs:
		image: ndarray to insert psfs into. Default is random noise.

	Outouts:
		none
	"""

	# Creates composite image of image passed in + psf image
	psf = getpsf_hst("../data/wfc3psf_248_267_50_F350LP_5_00.fits")
	xcen_arr = np.array([45.0,30.45678,47.10])
	ycen_arr = np.array([43.7,60.0,43.5])
	height_arr = np.array([1000.0,234.0,100.0])
	image_npsfs = insertpsf_n(image = image, psf = psf, xcens = xcen_arr,
				  ycens = ycen_arr, heights = height_arr)

	# Plots the image
	plt.figure()
	plt.imshow(image_npsfs, cmap = "hot", interpolation = "nearest")
	plt.colorbar()
	#plt.title("")
	plt.show()
	plt.close()

def cd_convolve(image):
	"""
	A function to model WC3 charge diffusion effects.

	Inputs:
		image: ndarray, image to apply CD effects to.
	
	Outputs:
		image with CD implemented
	"""

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
	"""
	A function which creates an image based on the parameters dataframe. This scales
	well to adding n psfs in the image.

	Inputs:
		paramsdf: dataframe containing all of the parameters.
	
	outputs:
		final image with psfs added and CD implemented
	"""
	
	# Starting with blankl image
	blank = np.zeros((100,100))

	# Getting all of the inputs for insertpsf_n ready
	# If needed getpsf_hst can be replaced with getpsf_2dgau to fit have gaussian psfs instead
	# In the future, make this read the run_props dictionary to know which function to use
	psf = getpsf_hst("../data/wfc3psf_248_267_50_F350LP_5_00.fits")
	xcen_arr = paramsdf["xcen"].values
	ycen_arr = paramsdf["ycen"].values
	height_arr = paramsdf["heights"].values
	return insertpsf_n(blank, psf = psf, xcens = xcen_arr, ycens = ycen_arr, heights = height_arr)

