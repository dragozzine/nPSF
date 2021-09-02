# TinyTim_PSFs

This repository contains the python module "make_psf" which can be used to generate TinyTim PSFs with the updated set of optical parameters determined by Gillis *et al.* (2019).

Preparation
-----------

1. Download and build the Tiny Tim executable. Version 7.5, which was used in Gillis *et al.* (2019) is available from http://tinytim.stsci.edu/static/tinytim-7.5.tar.gz. Note the directory where you build it.
2. Either clone this repository and add its location to the PYTHONPATH of your project, or download the file make_psf.py and place it in your project.
3. Edit your local copy of make_psf.py to change the default values at the top of the file to suit your purposes. In particular, it is recommended that you at minimum update default_tinytim_path to the location of
   your installed copy of Tiny Tim.
   
Use
---

Within your project, where you desire to generate a PSF, import the function make_subsampled_model_psf with code such as:

		from TinyTim_PSFs.make_psf import make_subsampled_model_psf

You can then call this function with code such as:

		make_subsampled_model_psf(filename="my_subsampled_psf.fits",
		                          psf_position=(14,570),
								  focus=-2.3,
								  chip=1)

This will create a subsampled PSF image and save it to the specified filename.

Other arguments are also allowed for this function; see the documentation of it (e.g. through help(make_subsampled_model_psf) within Python) for a full list.