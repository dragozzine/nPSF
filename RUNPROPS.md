Description of the different input parameters used in runprops.txt files
Example values and names are given for each parameter in this document, but these should be changed in the user's runprops.txt file to reflect their unique run.


## Initialization properties:

"input_image": "idsj03wgq_flc.fits",
    # Insert the name of your _flc.fits file to run nPSf on

"image_path": "../data/2013_FY27/",
    # Here you want the path from src to your image file, i.e. "../data/2013_FY27/"
    
"starting_guess": "startguess.csv",
    # Do not change
    
"output_filename": "results.csv",
    # Do not change
    
"names_dict": {
    "name_1": "2013 FY27",
    "name_2": "Moon1"
    },
    # The names_dict containes a name for each object contained in your image, which will match the number of psfs of your run. For example, if I had a 3 psf run, I would add a ',' after "Moon1" and an additional line that says "name_3": "Moon2" directly between "name_2": "Moon1", and '},'
    
#### Information on walkers (generaly speaking, the values given below are pretty good, but it can be important to vary the nsteps and nburnin):

"nwalkers": 100,
    # The number of walkers you wish to use. Generally speaking, 100 is good
    
"nsteps": 2500,
    # The number of steps you wish each walker to take
    
"nburnin": 2500,
    # The amount of burnin you wish to use (starting steps that will be discarded in the data)
    
"clustering_burnin": 100,
    # How many steps to use during clustering

#### Additional Run Info:

"npsfs": 2,
    # The number of psfs for which you are running
    
"folder": false,
    # Do not change
    
"RunGoal" : "Provide a runprops template",
    # Provides a reminder to future you of why you did this run



## TinyTim input parameters (Used by TinyTim to create the PSF):

"det_int": 22,
    # Specifies the camera being used. 22 is the detector number for WFC3 UVIS. Currently, nPSF is only designed for this camera
    
"xpos": 267,
    # In ds9, this is the y position of the center of your primary target (x and y being flipped for nPSF use)
    
"ypos": 248,
    # In ds9, this is the x position of the center of your primary target
    
"fmin": -8.0,
    # Fmin and fmax are the minimum and maximum extent that it should check focus values (The focus of HST changes over time due to expansion/contraction from heat)
    
"fmax": 8.0,
    
"fspace": 0.1,
    # Leave alone, not worth experimenting with
    
"sample_factor": 8,
    # Must be integer between 1-10. This asks to create the PSF in a higher resolution than the camera actually has (which makes it easier to do stuff with the PSF). A sample factor of 8 gives 64 sub-pixels within each pixel. This factor of 8 seems to work well
    
"psf_size": 2.0,
    # Don’t change (tells tiny tim how big to make the psf)



## nPSF settings:

"processors": 1,
    # Don't change, does nothing
    
"stamp_size": 60,
    # (stamp_size >= 60) e.g. cuts out a chunk that’s “60”x”60” upon which nPSF will perform
    
"stamp_x": 237,
    # Stamp_x and stamp_y are where to start the stamp (bottom left corner, subtract half the stamp_size from xpos and ypos). Will want to change these depending on the characteristics of the object being analyzed. Will want padding around the object as well
    
"stamp_y": 218,
    
"noise_cutoff": 3.0,
    # The central pixel in the PSF has to be a number of sigma above the noise in the image, where the number is the "noise_cutoff" value. (I.e. this drops the noise in the image)
    
"maskcheck": true,
    # Keep this as true
    
"cr_clean": true,
    # Clean cosmic rays from image, true or false (keep it lowercase)   
    
"cr_objlim": 5.0,
    # Checks for cosmic ray pixels and fixes them (it might accidentally think your target is a cosmic ray (if you look at your cr_cleaned image and there is no bright object where your target is, this has happened. If it does, you need to change this value to a larger value)
    
"use_clustering": true,
    # Keep this as true
    
"max_prune": 0.8,
    # Maximum that the clustering algorithm can remove. Probably don't need to change this



## Empirical settings (Ignored unless running an empirical psf program):

#### Settings used to cut out artifacts in your manufactured empirical psf:
"emp_size": 30,
    # Tells nPSF the size of the empirical psf (can make this smaller than actual size to cut out artifacts)

"psf_x": 0,
    # psf_x and psf_y adjust the location of the left corner of the empirical psf (to cut out artifacts). Make sure you change emp_size accordingly so that your range still fits within your empirical psf image.

"psf_y": 0,



## Map settings (Ignored unless running a map program):

"map_identifier": "super_fine",
    # Add an additional map identifier name to avoid replacing old files
    
"sep_grid_adj": 100,
    # Adjust the number of likelihood values that appear on the seperation grid (the bigger the number, the more values that show up)
    
"dx-dy_value": 0.5
    # Set the separation between x and y when running the grid
