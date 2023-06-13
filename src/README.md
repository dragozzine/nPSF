# src-nPSF
nPSF is designed to be a code that returns a likelihood map for the relative positions of multiple (n) PSFs 
in an image. It's goal is to provide likelihood information for relative astrometry, e.g., KBO moons.
Physics 529 Class Project Winter 2020
Darin Ragozzine and Benjamin Proudfoot
Students: Ian Clark, Jarrod Hansen, Jake Jensen, Maya Laker, Rochelle Steele, maybe others
Brigham Young University
contact darin_ragozzine@byu.edu for more information

If you write a new script, you will need to manually add it to this README file.


### src Scripts

3plot_compare: Prints a plot of the cleaned image, bestfit, and residuals in a 3 psf run as individual plots side by side.


analysis: Runs all of the analysis on the chains from npsf_run.

    plots:
    
    plot_best_fit:
    
    plot_comparison:
    
    map_plots:
    
    likelihood_map_chain:
    
    likelihood_map_grid:
    
    latlon_map:
    
    auto_window: Automatically sets a window size for the user
        Inputs:
            taus: 
            c: 
        
        Outputs:
            len(taus)-1 or np.argmin(m): int,
        
    autocorr_new: Sets a new window for auto_correlation
        Inputs:
            y:
            c: int, Default = 5, 
        
        Outputs:
            taus[window]: 
        
    auto_correlation:
    
    spec:
    
    geweke:
    
    testconvergence_geweke:
    
    make_corner_plot:
    
    make_walker_plots:
    
    make_likelihood_plots:
    
    make_bright_sep:
    
    make_obsdf:
    
    make_sigsdf:
    
    
analysis_emp:


clustering: Implements the clustering algorithm described in Hou 2010

    clustering: 
        Inputs:
            sampler:
            state:
            paramnames:
            log_likelihood:
            backend:
            image:
            psfs:
            focuses:
            runprops:
            const: , Default = 50, 
            lag: , Default = 10, 
            max_prune_frac: , Default = 0.9, 
        
        Outputs:
            sampler:
            state:
            

getpsf: 

    getpsf_2dgau: Creates a psf with a gaussian profile centered at the center of the output array/image
        Inputs:
            size: tuple of size 2 with each entry being the length of an axis
            cov: a covariance matrix specifying the shape of the gaussian profile
            supersample: a ratio at which to super sample the psf
                    for example: specifying a 21x21 image with a super sample
                            ratio of 5, will return a 105x105 image
                            
        Outputs:
            psf/psf.sum(): array/image of a psf with a gaussian profile
            
    getpsf_hst: Creates an array/image from a specifed TinyTim file.
        Inputs:
            filename: a string specifying the location of the TinyTim psf
            
        Outputs:        
            psf: an array/image of an HST psf
            
    getimage_hst: Creates an array/image from a specifed fits file.
        Inputs:
            filename: a string specifying the location of the image
            
        Outputs:
            psf: an array/image of anHST image
            filter: filter type from the image header
            nchip: CCD chip number from the image header
            bunits: units from the image header
            
    test_getpsf_2dgau: Runs a simple test for getpsf_2dgau and plots an image of the psf.
        Inputs: 
            size: tuple of size 2 specifying the lengths of each axis in the final image
            cov: a covariance matrix specifying the shape of the 2d gaussian
            supersample: a ratio at which to supersample the psf
            
        Outputs:
            none
            
    test_getpsf_hst: Runs a simple test for getpsf_hst and plots an image of the psf.
        Inputs:
            filename: a string specifying the location of the TinyTim psf.
            
        Outputs:
            none
          
          
guess_test:

    make_test_df:
        Inputs:
            none
        
        Outputs:
            test_df:
            
            
init_guess: 
           
    init_guess: This function will produce the initial guess used in nPSF.
        Inputs: 
            start_guess_df: A dataframe that has the columns as the parameter names, the first row is the mean value for each 
        parameter, and the second row is the standard deviation for each parameter.
            nwalkers: number of walkers. Added by Rochelle to make compatible with run_props in npsf_run.py 4/17/20
            
        Outputs: 
            params_df: A parameters dataframe with the same column names as start_guess_df and nwalker rows drawn from the 
        distribution.
        
        
insertpsf: 
    
    insertpsf_one: Takes a psf created by getpsf.py and inserts it into an input image at the location specified with the correct supersampling ratios. 
        Inputs: 
            image: an ndarray to insert the psf into (default is 100x100 array of zeros)
            psf: an ndarray created by a function in get_psf
            xcen: float, the x-position to place the psf
            ycen: float, the y-position to place the psf
            psfscale: int, the super sampling ratio to work with
            psfheight: float, the total integrated flux of the psf
        
        Outputs:
            (image + psfimage) or (image + cd_convolve(psfimage, runprops)): ndarray of the combined image + psf
            
    test_insertpsf_one: Runs a simple test for insertpsf_one and plots the result
        Inputs:
            image: ndarray of image to insert a test psf into (default is a 100x100 array of random noise)
        
        Outputs:
            none
            
    insertpsf_n: A function which provides a wrapper for looping over insertpsf_one.
        Inputs: 
            image: ndarray to start with. Default is 100x100 array of zeros.
            psf: ndarray, created using getpsf.py with which to create the image
            xcens: ndarray with size 1x(number of psfs to add) of the x-positions of the psfs
            ycens: ndarray with size 1x(number of psfs to add) of the y-positions of the psfs
            heights: ndarray with size 1x(number of psfs to add) of the heights of the psfs
        
        Outputs: 
            image: an ndarray with psfs added to the input image. Note that the output image has been convolved with the charge 
        diffusion kernel.

    test_insertpsf_n: Runs a simple test for insertpsf_n and plots the result.
        Inputs: 
            image: ndarray to insert psfs into. Default is random noise.
            
        Outputs:
            none
            
    cd_convolve: A function to model WC3 charge diffusion effects.
        Inputs: 
            image: ndarray, image to apply CD effects to.
            
        Outputs:
            scipy.ndimage.convolve( image, cd_kernel, mode='constant', cval=0.0 ): image with CD implemented
            
    make_image: A function which creates an image based on the parameters dataframe. This scales well to adding n psfs in the image.
        Inputs: 
            paramsdf: dataframe containing all of the parameters.
        
        Outputs: 
            insertpsf_n(blank, psf = psf, xcens = xcen_arr, ycens = ycen_arr, heights = height_arr): final image with psfs added and CD implemented
            
            
latlon_transform: Converts RA and DEC values to latitude and longitude coordinates for a given object and its moons. 

    convert_to_primary_centric: This function takes a parameter Dataframe in RA/DEC, and converts it to Latitude and Longitude, while also converting the dates to Primary-Centric Julian dates
        Inputs:
            paramsDF: A dataframe of the observed positional data of the KBO in question
            objectName: The name of the object being observed (needed for the Horizons function)
        
        Outputs:
            _obs_df.csv: saves an obs_df file for use in Multimoon
            forsigsDF: returns a dataframe to analysis.py for inclusion in a sigsdf file
        
            
likelihood: This script includes various functions for calculating likelihoods
    
    log_likelihood: 
        Inputs:
            parameters:
            image:
            psfs: 
            focuses:
            runprops:
            plotit: , Default = False, 
        
        Outputs:
            bestfit.png: 
            residuals.png:
            llhood.png:
            loglike: 
            
    generate_bestfit_residual
        Inputs:
            parameters:
            image:
            psfs:
            focuses:
            runprops:
            plotit: , Default = False, 
            
        Outputs: 
            psfimage:
            residuals:
            xcens:
            ycens:
            
    log_likelihood_map:
        Inputs:
            psf1params:
            psf2loc:
            psf2heights:
            image:
            psf:
            runprops:
            
        Outputs:
            llhoods:
            
    log_probability:
        Inputs:
            image:
            parameters:
            
        Outputs:
            lp + log_likelihood(image,parameters): 
             

m-a_sep_TNO: modified from Benjamin Proudfoot's m-a_sep_Borasisi to provide to take multimoon J2 values and provide relevant graphs of hierarchical triple detectability for various TNOs. 

        Inputs:
            runprops: txt file providing the names of the TNOs to be plotted
            objectdata: txt file containing individual TNO characteristic values
            sigsdf_mm: multimoon csv file containing the results data values
            
        Outputs: 
            mr_x_sma_"".pdf plots: plots containing all of the relevant detectability data
            

makesynth:


maketestfiles:


map_run:


map_run_fixedgrid:


multi_map:


multi_run:


npsf_empirical:


npsf_run: This file runs the nPSF program. It pulls together all of the functions and files in each the project and runs them using the information given by the user in the runprops and startguess files.


npsf_test:


nPSFtest:


params:


priors:


remake_map:


remake_plots:


testing_and_examples:


testsupport:
            
            
SUPPORT FILES:
runprops_src.txt: This text file contains all of the variables that will be used in nearly every function of nPSF. It is read in as a dictionary using JSON methods.

startguess_src.csv: This csv file contains the startguesses for the positions of each psf being utilized in the current run. Should be in pixels accordingly fitted to the stamp size designated in runprops.txt. 

objectdata_src.txt: This text file contains additional parameters to be used exclusively by `m-a_sep_TNO.py`. It is read in as a dictionary using JSON methods.

image.csv: 

synthprops.txt: