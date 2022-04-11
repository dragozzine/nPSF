#
# multi_run.py
#
# Takes a folder of HST images and runs npsf on them
#
# Benjamin Proudfoot
# 01/24/22
#

import numpy as np
import emcee
import sys
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import corner
import commentjson as json
from tqdm import tqdm
from params import params_to_fitarray
from params import fitarray_to_paramsdf
from init_guess import init_guess
from likelihood import log_likelihood # Jarrod's code for EnsembleSampler
from guess_test import make_test_df
from getpsf import *
from analysis import *
from clustering import *
import tinytim_psfs.make_psf
from csv import writer
import os.path
import datetime
from astropy.io import fits
# from params.py import df_to_array, npsf_init_guess()
from schwimmbad import MPIPool
import shutil

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data


def npsf_run(runprops):
    ### Load run_props from JSON file ###
    #runprops = ReadJson("runprops.txt"npsf_run(runprops)).outProps()

    nwalkers = runprops.get("nwalkers")
    nsteps = runprops.get("nsteps")
    nburnin = runprops.get("nburnin")
    output_filename = runprops.get("output_filename")

    runprops["best_likelihood"] = -np.inf

    # Create results folder
    x = datetime.datetime.now()
    date = str(x.strftime("%Y"))+"-"+str(x.strftime("%m"))+"-"+str(x.strftime("%d"))+"_"+str(x.strftime("%H"))+"."+str(x.strftime("%M"))+"."+str(x.strftime("%S"))
    resultspath = "../results/"+date+"_"+runprops.get('image_name')+"_"+str(runprops.get("npsfs"))+"psf"
    if not os.path.exists(resultspath):
        os.makedirs(resultspath)
    runprops["resultspath"] = resultspath

    shutil.copy("runprops.txt", resultspath + "/runprops.txt")

    # Get initial guess (with init_guess.py)
    # This uses test data from guess_test.py
    params_df = pd.read_csv(runprops.get("starting_guess"),sep=',',index_col=0)

    p0_df = init_guess(params_df, nwalkers)
    paramnames = params_df.columns.tolist()
    p0, change_dict = params_to_fitarray(p0_df)

    # Get ndim
    ndim = np.shape(p0)[1]

    # Loading in image to be solved and making a small postage stamp version
    x = runprops.get("stamp_x")
    y = runprops.get("stamp_y")
    size = runprops.get("stamp_size")

    f = runprops.get('input_image')
    imageraw = getimage_hst(f)[x:x+size,y:y+size]

    # Clean cosmic rays from image (maybe this can be removed when guesses are good enough?)
    # This may also be irrelevant if we move to simultaneous pair fitting
    import ccdproc
    cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim"), gain_apply = False)
    image = np.array(cr_cleaned)

    # Make sure that the CR algorithm hasn't marked the target as a cosmic ray!
    reset = False
    while True:
        crredo = False
        for i in range(nwalkers):
            #print(i,crmask[int(p0[i,0]),int(p0[i,1])])
            if crmask[int(p0[i,0]),int(p0[i,1])]:
                crredo = True
        if crredo:
            cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim") + 1.5, gain_apply = False)
            if reset:
                print("CR rejection algorithm is flagging your target as a CR. Aborting run.") 
                print("Consider increasing the objlim in runprops.")
                #sys.exit()
                return -1, resultspath
            reset = True
        else:
            break
    image = np.array(cr_cleaned)

    plt.figure()
    plt.imshow(image, cmap = "hot", interpolation = "nearest", origin = "lower")
    plt.colorbar()
    plt.scatter(p0[:,0].flatten(), p0[:,1].flatten(), color = "blue", marker = "x", s = 5, alpha = 0.2)
    plt.savefig(resultspath + "/cleanedimage.png")
    plt.figure()
    plt.imshow(crmask, cmap = "hot", interpolation = "nearest", origin = "lower")
    plt.colorbar()
    plt.scatter(p0[:,0].flatten(), p0[:,1].flatten(), color = "blue", marker = "x", s = 5, alpha = 0.2)
    plt.savefig(resultspath + "/crmask.png")
    plt.close("all")

    # Ensure there are no negative pixels. Add constant offset, which will be corrected in model images.
    if np.nanmin(image) < 0:
        image = image - np.floor(np.nanmin(image)) + 1.0

    # Calculating image noise characteristics for priors
    runprops["med_noise"] = np.median(image)
    runprops["std_noise"] = np.std(image)

    # Getting inputs for Tiny Tim PSFs
    fmin = runprops.get("fmin")
    fmax = runprops.get("fmax")
    fspace = runprops.get("fspace")
    focuses = np.arange(fmin, fmax + fspace, fspace)
    xpos = runprops.get("xpos")
    ypos = runprops.get("ypos")
    size_psf = runprops.get("psf_size")
    sample_factor = runprops.get("sample_factor")
    nchip = runprops.get("chip")
    ndet = runprops.get("det_int")
    filter = runprops.get("filter")
    numpsfs = focuses.size

    # Making TinyTim PSFs (this is skipped if all the PSFs have previously been made)
    for i,focus in enumerate(focuses):
        filename = "modelpsfs/wfc3psf_" + str(ndet) + "_" + str(nchip) + "_" + filter + "_" + str(xpos) + "_" + str(ypos) + "_" + str(round(focus,1)) + "_" + str(size_psf) + "_" + str(sample_factor) + ".fits"
        if not os.path.exists(filename):
            tinytim_psfs.make_psf.make_subsampled_model_psf(filename, psf_position = (xpos, ypos), focus = focus, chip = nchip, detector = ndet,
					   filter_name = filter, psf_size = size_psf, subsampling_factor = sample_factor)
        else:
            pass

    # Putting PSFs into array
    testpsf = getpsf_hst(filename)
    psfs = np.empty((testpsf.shape[0],testpsf.shape[1],numpsfs))
    for i, focus in enumerate(focuses):
        filename = "modelpsfs/wfc3psf_" + str(ndet) + "_" + str(nchip) + "_" + filter + "_" + str(xpos) + "_" + str(ypos) + "_" + str(round(focus,1)) + "_" + str(size_psf) + "_" + str(sample_factor) + ".fits"
        psfs[:,:,i] = getpsf_hst(filename)

    # Getting charge diffusion kernel from tiny tim psf header
    head = fits.open(filename)[0].header
    line1 = [float(x) for x in head["COMMENT"][3].split()]
    line2 = [float(x) for x in head["COMMENT"][4].split()]
    line3 = [float(x) for x in head["COMMENT"][5].split()]
    cd_kernel = np.array([line1,line2,line3])
    runprops["cd_kernel"] = cd_kernel

    # check for bad walkers / bad initial guesses
    reset = 0
    for i in range(nwalkers):
        #print(p0[i,:])
        llhood = log_likelihood(p0[i,:], image, psfs, focuses, runprops)
        #print(i, p0[i,:], llhood)
        while (llhood == -np.Inf):
            # Resetting parameters
            p = random.random()
            p0[i,:] = p*p0[random.randrange(nwalkers),:] + (1-p)*p0[random.randrange(nwalkers),:]
            llhood = log_likelihood(p0[i,:], image, psfs, focuses, runprops)
            #print(reset, llhood, i, p0[i,:])
            reset += 1
            if reset > 10000:
                sys.exit()

    # Setting up emcee sampler object
    backend = emcee.backends.HDFBackend(resultspath + "/chain.h5")
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood, args = [image, psfs, focuses, runprops], backend = backend)

    # Burnin
    print("Beginning burnin")
    state = sampler.run_mcmc(p0, nburnin, progress = True)

    # Clustering
    if runprops.get("use_clustering"):
        sampler, state = clustering(sampler, state, paramnames, log_likelihood, backend, image, psfs, focuses, runprops, max_prune_frac = runprops.get("max_prune"))
    sampler.reset()

    # Sample
    print("Beginning sample")
    sampler.run_mcmc(state, nsteps, progress = True)

    # Run analysis and make plots
    plots(sampler, resultspath, runprops)
    bestfit = plot_best_fit(sampler, image, psfs, focuses, runprops)

    return bestfit, resultspath


# Machinery for running from a folder of images
import os
import glob
cwd = os.getcwd()
print(cwd)

# Make sure you are in a data directory
if "multi_run.py" in sys.argv[0]:
    if "data" in cwd:
        folder = cwd
        images = glob.glob("*.fits")
        os.chdir("../../src")
    else:
        print("When running multi_run.py, make sure you are located in the data folder")
        sys.exit()

    print(sys.argv[0])
    # Loop over each image
    for image in images:
        fullpath = folder + "/" + image
        print(fullpath)
        runprops = ReadJson("runprops.txt").outProps()
        runprops["input_image"] = fullpath
        runprops["image_name"] = image[:-5]

        # Run nPSF
        npsf_run(runprops)