# User's Guide to nPSF

### Purpose of the User's Guide

This file is designed to outline the user process of nPSF to a new user and to provide a reference for the various outputs nPSF can generate. Additionally, alternate uses of nPSF are outlined in the `## Alternative Functions of nPSF` section of this document. 

### Purpose of nPSF

    nPSF is designed to be a code that returns a likelihood map for the relative positions of multiple (n) PSFs in an image. It's goal is to provide likelihood information for relative astrometry, e.g., KBO moons.

Files needed to run nPSF:
    1) A flc.fits file of the observational image containing the system you wish to fit psfs to
    2) A csv file of the inital guess for the position of the object 
    3) A txt file for the run properties of the specific run information you wish to accomplish.

The necessary elements of these files are explained below.
  
  
### Observational Image,  "IDnumber"_flc.fits

We generally obtain the flc.fits files we use for nPSF from MAST:`https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html`.
These files can then be added to a newly created directory for the target object within nPSF's `data` directory, e.g. `nPSF/data/2013_FY27`, where `2013_FY27` is the target object.

### Initial Guess, startguess.csv

An initial guess file must also be created within the target object's `data` directory. A base copy of the startguess file can be found within the `src` directory, titled `startguess_src.csv`. Copies of this populated to data directories should be renamed as `startguess.csv`. The parameter set required in the startguess file includes `xpos_#`, `ypos_#`, `height_#`, and `focus` (only a single focus parameter should be included, reguardless of the number of psfs being fitted). A 2 psf fit will require 2 sets of parameters, a 3 psf fit will require 3 sets of parameters, etc. Additional sets of parameters must be added manually by the user. 

The startguess parameters should be supplied by the user as an estimate to where the target objects are located in the HST image (in pixels) according to stamp size restrictions supplied in the `runprops.txt`. The file is arranged into 3 rows, the first identifying the name of each parameter, the second the mean value (or best guess position), and the third the standard deviation of error to initialize the code with.

### The Run properties, "runprops.txt"

All of the variables necessary for running nPSF are defined in the runprops dictionary, which should also be included in the target object's `data` directory. A base file labeled `runprops_src.txt` can be found in `src`. Copies of this populated to data directories should be renamed as `runprops.txt`.The runprops.txt file is set up in such a way to display the whole dictionary in a readable format, and should look familiar to anyone versed in JSON. A guide for all runprops inputs is located in the `RUNPROPS.md` file.

### How to run nPSF

Once all required packages are installed, nPSF is generally run from the command line. From inside the data directory that contains the above 3 necessary files, (e.g. `nPSF/data/2013_FY27`) the following command can be used to run nPSF:

`python ../../src/npsf_run.py`

The first time nPSF is run from a data directory, a new directory will be automatically created in `nPSF/results` with the same name as the utilized data directory. This results directory will hold an additional unique directory which saves all data from your run, and it will continue to be populated by directories saving results from any further runs you may execute.

Additionally, a first time run will include a series of fast moving text that indicates TinyTim is making the PSFs needed for the run. I believe this can take between 5-10 minutes. After that, nPSF will initiate the MCMC process. 

When a run is initiated and after any PSF generation output by TinyTim, a series of values will appear in the terminal:

    ndim: This is an indicator of how many PSFs you are fitting to the image. 4 dimensions indicates a single PSF, 7 dimensions for 2 PSFs, and 10 dimensions for 3 PSFs.
    
    bunits: The type of units the image is in. For the way nPSF is set up, the image units need to be in ELECTRONS for things to work properly.
    
    filter: The filter in which the image was taken. Placing it here makes it easier to identify, rather than pulling up the image header in ds9.
    
    CCD chip: The CCD camera chip number.
    
    cr_clean: Indicates whether you set cosmic ray cleaning to true or false in the runprops. This merely helps you catch whether you put it to what you wanted or not. 
    
    med_moise: The median value of all pixels in the image.
    
    std_noise: The square root of the median noise. 
    
    A single number value follows these results, which I don't know the purpose of, but is the product of the std_noise and the noise_cutoff divided by 0.1.

While nPSF is running, a loading bar will appear showing the progress and estimated time remaining in that step. Once nPSF is done running, a variety of plots will be output to the results folder.

### The Results Folder

If for some reason an expected plot is not created or is created disfunctional, the user may run `remake_plots.py` from the results directory experiencing this issue. Note that a `chain.h5` file for the run must be present to do this.

Convergence Plots:

    walkers_pngs: A folder containing individual files of each parameter's walk.

    walkers.pdf: A single pdf of each parameter's walk over the number of steps taken. 

    walkers_full.pdf: A single pdf of each parameter's walk over the combined burnin, clustering, and number of steps taken.


Image PNGs: 

    bestfit.png: Recreates an image of the psfs using the best fit values found during the nPSF run.

    cleanedimage.png: An image of the psfs after it has been cleaned of cosmic rays. The blue dots represent the starting guess positions of this specific nPSF run. 
    
    crmask.png: The image of the mask that cleans the cosmic rays out of your image to give you the cleanedimage.png.
    
    llhood.png: Ask Ben or Ragozzine
    
    residuals.png: Shows the difference between the real data psf image and the fitted psf image. 
    
    trimmedpsf.png: Unique to empirical psf runs, this shows an image that has been trimmed to negate artifacts in the supplied empirical psf according to the additional parameters "emp_size","psf_x", and "psf_y" included in the empirical settings of the runprops. 

   Work in progress, only appear if you run 3plot_compare.py (trying to get better residuals results, an example of these plots can be found on Haumea in `/home/byu.local/wgiforos/research/nPSF/results/2005_EO304/2022-07-30_10.31.40_idy61hs5q_flc_3psf`):
   
        compare_plots.png: Shows the full image of the cleaned image, the bestfit image, and the residuals image in a row

        compare_plots_normalized.png: Same as compare_plots.png, but was my attempt to normalize the residuals image. I'm pretty sure I did not do this correctly, so it should be fixed if these are sought after again.

        compare_plots_zoomed.png: Shows a grid of images according to the format of compare_plots.png, but zoomed in so that each row shows a different object in the three images.   


Analysis Plots:

    bright_sep_lim.pdf: Not currently functional, needs to be fixed.
    
    corner.pdf: The corner plot shows all the one and two dimensional projections of the posterior probability distributions of your parameters. This is useful because it quickly demonstrates all of the covariances between parameters.

    cornerderived.pdf: A corner plot that includes the additional derived parameters. 
    
    likelihoods.pdf: Plots the Log(likelihood) vs. each parameter of the walk and each derived parameter. The plots to the above and right of the main plot are both histograms that you should want to look as gaussian as possible. Each color represents a different walker. The highest point on the plot represents the best fit value of that parameter. A well explored run will be shown in this plot as a smooth-ish mountain/hill with a good mix of colors (though a good run has a natural green distribution among the other colors).


Information Files: 

    ..._obs_df.csv: Contains the latitude and longitude positions of the primary object and its secondaries. Purposed to be used in Multimoon.

    sigsdf.csv: Contains the posterior values of the run as well as valuable calculated information. 
    
        Posterior values:
        
            x: Positional pixel parameter in the x-direction.
            
            y: Positional pixel parameter in the y-direction.
            
            h: Height or brightness of the psf.
            
            f: Focus of the HST camera. This is always changing, so we include it in our parameters to be solved for.
            
            dra: Derived difference between the RA of the TNOs.
            
            ddec: Derived difference between the declination of the TNOs.
            
            dmag: Derived difference between the brightness magnitude of the TNOs.
            
            sep: Derived separation distance between the TNOs in x and y pixel coordinates.
            
            pa: Position angle, the direction in which one object lies relative to another on the celestial sphere, measured in degrees from north in an easterly direction.
            
            dx: Difference between the TNO x-values.
            
            dy: Difference between the TNO y-values.
            
            lat: Latitude of the primary TNO.
            
            long: Longitude of the primary TNO.
            
            dlat: Difference between the TNO latitudes.
            
            dlong: Difference between the TNO longitudes.
        
        Many of these values will be labeled with either a 1, 2, or 3 to indicate the object in question. For parameters calculating the difference between objects, a 2psf run will yield values with a 1 and nothing more, such as dra1. For 3psf runs and higher, the difference between the primary and secondary is labeled dra2, the difference between primary and tertiary is dra3, etc.

Additional Files:

    chain.h5: Contains all values found during the nPSF run. Can't be opened directly, but is accessible through python.

    runprops.txt: A copy of the `runprops.txt` file used for this run.

    startguess.csv: A copy of the `startguess.csv` file used for this run.
    
    
## Alternative Functions of nPSF

Optimizer:
    Instead of doing a full nPSF run, the user can simply attempt to optimize the likelihood to get a guess of the correct parameter values by marking `use_optimizer` as true in the runprops. The same runprops and startguess files as `npsf_run.py` can be used for this run according to the additional `Optimizer settings` outlined in `RUNPROPS.md`. 
    
    optimize.py: script used by npsf_run.py when use_optimizer is set to true in the runprops.
    
    optimizer_sigsdf.csv: similar to 'sigsdf.csv' defined above, but without the standard deviation errors and with an additional JD time output at the top of the column. 


Empirical PSFs:
    This function uses an empirical psf provided by the user rather than generating psfs with TinyTim. The focus parameter is also ignored in this type of run. The empirical psf must be included in the data directory as a .fits file along with the image, runprops, and startguess. This empirical psf function uses the same runprops and startguess type files as `npsf_run.py`, but additionally makes use of the runprops `Empirical settings`, which are outlined in the `RUNPROPS.md` file.
    
    npsf_empirical.py: To use an empirical psf, use this function rather than npsf_run.py. It runs similarly, but skips TinyTim psf generation and ignores the focus parameter. 
    
    trimmedpsf.png: output of the empirical psf used. This psf is trimmed according to the 'Empirical Settings' as outlined in 'RUNPROPS.md'.
    

Likelihood Map:
    Creates a grid of values of the given image from which a likelihood map is generated using interpolation techniques. This process never reached the interpolation accuracy we desired, and so the project stands unfinished. All map functions use the same runprops and startguess type files as `npsf_run.py`, but using the runprops `Map settings`, which are outlined in the `RUNPROPS.md` file.
    
    map_run.py: To create a likelihood map, run this function instead of npsf_run.py. It takes an HST image in a data directory and creates a likelihood map for it. 
    
    map_run_fixedgrid.py: This script was created to better test whether our interpolation function was working, as map_run.py is set up so that a finer grid doesn't necessarily have the same points in the grid as a less fine grid. Using a fixedgrid helped combat this issue. Besides this, it runs the same as map_run.py.
    
    multi_map.py: Takes a folder of HST images and creates a liklihood map for each one. I don't know if I ever tested this script, so I can't say whether it works now or not, though I did make the same updates to it as I did to map_run.py to try to keep it up to date. 
    
    remake_map.py: Used to regenerate the output plots for a likelihood map run. This script should be run from a results directory. A chain.h5 file must be present in the results directory to do this. 

    ..._map: The results folder created for map runs, where the '...' is filled by the name of your object's data folder as documented in the runprops under 'image_path'.
    
    dx_dy_llhoods.png: The likelihood map plotted in delta x vs. delta y with a color scale of the likelihoods.
    
    dlat_dlong_llhoods.png: The likelihood map plotted in delta latitude vs. delta longitude with a color scale of the likelihoods.
    
    grid.npy: A loadable file of the x and y values of the grid to allow for easier access when testing the effectiveness of the likelihood map functions. 
    
    llhoods.npy: A loadable file of the llhood values of the grid to allow for easier access when testing the effectiveness of the likelihood map functions. 
    

m-a_sep_TNO.py:
    Indicates hierarchical triple detectability and for the input of any TNO system. Outputs a plot of mass ratio vs. semi-major axis. This script should be run from the results directory and uses a normal runprops as input. 

    objectdata.txt: .txt file similar to runprops.txt that is exclusively used for m-a_sep_TNO.py and must be included in the results directory of the object being analyzed (deleting the '_src' bit at the end of the filename).

    sigsdf_mm.csv: a sigsdf file taken from the Multimoon output of the specific object being run. Contains the J2 values and masses needed for this code to run (you may need to rename the Multimoon file as sigsdf_mm.csv).

    mr_x_sma_....pdf: pdf output of the plot generated by m-a_sep_TNO.py. The darkly shaded area indicates where the tertiary object cannot exist in a hierarchical triple system. The lightly shaded region before the approximate dectability limit is where we would expect the tertiary object to be detectable by nPSF. 
    

### Additional Notes

1) The first psf should always be the main body object.


### Additional Instructions for Ragozzine students new to using nPSF (more step-by-step style)

1) Follow the steps in nPSF's README to install TinyTim. Do this in your Jupyterhub terminal. 

    Virtual Environments for beginners: Different software requires different versions of the same packages. The nice thing about a virtual environment is that you can have different environments with the same software packages installed as different versions, so instead of having to reinstall a package everytime you want to use a different software, you can instead just execute a command and you are in the correct environment to use the desired software. 
    To install your environment for nPSF, you will want to go to your user home directory (e.g. mine is `wgiforos@Haumea:~$`) and type the command ... -> Currently, you will need to look up how to do this yourself. Or ask someone for help.
    
For additional info on virtual environments, see: `https://docs.python.org/3/library/venv.html`. 

You can find your .bashrc in your user base directory by typing `ls -a` and pressing enter. 

2) Follow the information outlined in this users guide (ignore `###How to run nPSF` for now). 

    See: `https://docs.google.com/document/d/1XxJBPNaASCmn9kUZJDkNdUUgqLpHbnqxM2cPmaSHNqU/edit?usp=drive_link`
    
    When it comes to getting data from MAST, typically you want to enter the RA and DEC of your object on one of the nights observed in the data set you wish to pull images from. If you are working on KBOs, you can find this information by going to `http://www2.lowell.edu/users/grundy/tnbs`, clicking the hyperlink `status of the known binaries`, and then searching for your target TNO. Once you've clicked on it, you scroll down to Astrometric Observations and choose a night with HST/WFC3 observations, which includes the RA and DEC information for the chosen night. 

    Once you have copied and plugged this info into the MAST search bar, a list of data should show up. You can further refine this list by checking on the left the following: HST in the Mission category, WFC3/UVIS in Instrument, and your object's title in Target Name (additional parameters can be checked by you according to what you need, these just tend to do a pretty good job of narrowing it down). 

    When you've identified the observational data you're looking for, you can check the boxes for each night and click the download button (looks like an orange basket with some green arrows going into it). Then check the box for FLC in Group, select the files you want and hit download (green marker on save file image). Then you can upload those flc files to the correct data directory on nPSF. 
    
3) Become acquainted with the `RUNPROPS.md` file, which is a guide to using `runprops.txt`.

4) Once you have a data directory containing a `_flc.fits` file, a `runprops.txt` file, and a `startguess.csv` file, and they have been completely filled out according to your desired parameters, follow the instructions of `###How to run nPSF` found in this document.

5) Learn what your results mean in `### The Results Folder` section of this document.