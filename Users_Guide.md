# User's Guide to nPSF

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

 While nPSF is running, a loading bar will appear showing the progress and estimated time remaining in that step. Once nPSF is done running, a variety of plots will be output to the results folder.


### The Results Folder

Convergence Plots:

    walkers_pngs: a folder containing individual files of each parameter's walk

    walkers.pdf: a single pdf of each parameter's walk over the number of steps taken. 

    walkers_full.pdf: a single pdf of each parameter's walk over the combined burnin, clustering, and number of steps taken.


Image PNGs: 

    bestfit.png: recreates an image of the psfs using the best fit values found during the nPSF run.

    cleanedimage.png: An image of the psfs after it has been cleaned of cosmic rays. The blue dots represent the starting guess positions of this specific nPSF run. 

    compare_plots.png:

    compare_plots_normalized.png:

    compare_plots_zoomed.png:


Information Files: 

    ..._obs_df.csv: Contains the latitude and longitude positions of the primary object and its secondaries. Purposed to be used in Multimoon.

    sigsdf.csv: contains the posterior values of the run as well as valuable calculated information. 
        Posterior values:
            


chain.h5: Contains all values found during the nPSF run.


bright_sep_lim.pdf:

corner.pdf:

cornerderived.pdf:

crmask.png:

likelihoods.pdf:

llhood.png:

residuals.png:

runprops.txt: a copy of the `runprops.txt` file used for this run.



startguess.csv: a copy of the `startguess.csv` file used for this run.


    trimmedpsf.png:

### Additional Notes

1) The first psf should always be the main body object.


### Additional Instructions for Ragozzine students new to using nPSF (more step-by-step style)

1) Follow the steps in nPSF's README to install TinyTim. Do this in your Jupyterhub terminal. 

    Virtual Environments for beginners: Different software requires different versions of the same packages. The nice thing about a virtual environment is that you can have different environments with the same software packages installed as different versions, so instead of having to reinstall a package everytime you want to use a different software, you can instead just execute a command and you are in the correct environment to use the desired software. 
    To install your environment for nPSF, you will want to go to your user home directory (e.g. mine is `wgiforos@Haumea:~$`) and type the command 
For additional info on virtual environments, see: `https://docs.python.org/3/library/venv.html`. 

You can find your .bashrc in your user base directory by typing `ls -a` and pressing enter. 

2) Follow the information outlined in this users guide (except `###How to run nPSF`). 

    When it comes to getting data from MAST, typically you want to enter the RA and DEC of your object on one of the nights observed in the data set you wish to pull images from. If you are working on KBOs, you can find this information by going to `http://www2.lowell.edu/users/grundy/tnbs`, clicking the hyperlink `status of the known binaries`, and then searching for your target TNO. Once you've clicked on it, you scroll down to Astrometric Observations and choose a night with HST/WFC3 observations, which includes the RA and DEC information for the chosen night. 

    Once you have copied and plugged this info into the MAST search bar, a list of data should show up. You can further refine this list by checking on the left the following: HST in the Mission category, WFC3/UVIS in Instrument, and your object's title in Target Name (additional parameters can be checked by you according to what you need, these just tend to do a pretty good job of narrowing it down). 

    When you've identified the observational data you're looking for, you can check the boxes for each night and click the download button (looks like an orange basket with some green arrows going into it). Then check the box for FLC in Group, select the files you want and hit download (green marker on save file image). Then you can upload those flc files to the correct data directory on nPSF. 
    
3) Become acquainted with the `RUNPROPS.md` file, which is a guide to using `runprops.txt`.

4) Once you have a data directory containing a `_flc.fits` file, a `runprops.txt` file, and a `startguess.csv` file, and they have been completed filled out according to your desired parameters, follow the instructions of `###How to run nPSF` found in this document.

5) Learn what your results mean in `### The Results Folder` section of this document.