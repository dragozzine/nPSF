#
# 3plot_compare.py
#
# Prints a plot of the residuals in a 3 psf run as individual residual plots side by side. 
# Runs from an object's results directory
#
# William Giforos
# 08/18/22
#

import glob
import os
import sys
import commentjson as json
import emcee
import numpy as np
from getpsf import *
from init_guess import init_guess
from params import params_to_fitarray
from params import fitarray_to_paramsdf
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data
    
#Credit to Joe Kington for writing the subclass (MidpointNormalize), and to Rutger Kassies for pointing out the answer. https://stackoverflow.com/questions/57180317/making-sure-0-gets-white-in-a-rdbu-colorbar
class MidpointNormalize(mpl.colors.Normalize):
    ## class from the mpl docs:
    # https://matplotlib.org/users/colormapnorms.html

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def comparison(image, psfimage, residuals, xcens, ycens):
    # Make the triple image        
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(5, 3))
    im1 = axes[0].imshow(image, cmap = "hot", interpolation = "nearest", origin = "lower")
    axes[0].title.set_text("cleaned image")
    axes[0].tick_params(labelsize=5)
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.025)
    plt.colorbar(im1, cax=cax).ax.tick_params(labelsize=5)

    im2 = axes[1].imshow(psfimage, cmap = "hot", interpolation = "nearest", origin = "lower")
    axes[1].title.set_text("bestfit")
    axes[1].get_yaxis().set_visible(False)
    axes[1].get_xaxis().set_visible(False)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.025)
    plt.colorbar(im2, cax=cax).ax.tick_params(labelsize=5)

    im3 = axes[2].imshow(residuals, cmap = "hot", interpolation = "nearest", origin = "lower")
    axes[2].scatter(ycens, xcens, s = 0.5, color = "blue", marker = "o")
    axes[2].title.set_text("residuals")
    axes[2].get_yaxis().set_visible(False)
    axes[2].get_xaxis().set_visible(False)
    divider = make_axes_locatable(axes[2])
    cax = divider.append_axes("right", size="5%", pad=0.025)
    plt.colorbar(im3, cax=cax).ax.tick_params(labelsize=5)

    fig.tight_layout()
    fig.savefig(runprops.get("resultspath") + "/compare_plots.png", dpi = 300)
    plt.close()

    # Make the normalized triple image        
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(5, 3))
    im1 = axes[0].imshow(image, cmap = "hot", interpolation = "nearest", origin = "lower")
    axes[0].title.set_text("cleaned image")
    axes[0].tick_params(labelsize=5)
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.025)
    plt.colorbar(im1, cax=cax).ax.tick_params(labelsize=5)

    im2 = axes[1].imshow(psfimage, cmap = "hot", interpolation = "nearest", origin = "lower")
    axes[1].title.set_text("bestfit")
    axes[1].get_yaxis().set_visible(False)
    axes[1].get_xaxis().set_visible(False)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.025)
    plt.colorbar(im2, cax=cax).ax.tick_params(labelsize=5)

    im3 = axes[2].imshow(residuals, cmap = "hot", interpolation = "nearest", origin = "lower", norm=MidpointNormalize(midpoint=0))
    axes[2].scatter(ycens, xcens, s = 0.5, color = "blue", marker = "o")
    axes[2].title.set_text("residuals")
    axes[2].get_yaxis().set_visible(False)
    axes[2].get_xaxis().set_visible(False)
    divider = make_axes_locatable(axes[2])
    cax = divider.append_axes("right", size="5%", pad=0.025)
    plt.colorbar(im3, cax=cax).ax.tick_params(labelsize=5)

    fig.tight_layout()
    fig.savefig(runprops.get("resultspath") + "/compare_plots_normalized.png", dpi = 300)
    plt.close()


def make_subplot(data, title, ax, ycens, xcens, colmin, colmax, residuals=False):
    if residuals == False:
        im = ax.imshow(data, cmap = "hot", interpolation = "nearest", origin = "lower", vmin=colmin, vmax=colmax)
    else: 
        im = ax.imshow(data, cmap = "hot", interpolation = "nearest", origin = "lower") 
        ax.scatter(ycens, xcens, s = 5, color = "blue", marker = "o")
    if title != "None":
        ax.title.set_text(title)
    ax.tick_params(labelsize=5)
    ax.set_xlim(ycens-10, ycens+10)
    ax.set_ylim(xcens-10, xcens+10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.025)            
    plt.colorbar(im, cax=cax).ax.tick_params(labelsize=5)
    
def zoomed_comparison(image, psfimage, residuals, xcens, ycens, npsfs):
    # Make the triple image        
    fig, axes = plt.subplots(nrows=npsfs, ncols=3)
    for i in range(npsfs):
        if i ==0:
            make_subplot(image, "cleaned image", axes[0,0], ycens[0], xcens[0], image.min(), image.max())
            make_subplot(psfimage, "bestfit", axes[0,1], ycens[0], xcens[0], image.min(), image.max())
            make_subplot(residuals, "residuals", axes[0,2], ycens[0], xcens[0], image.min(), image.max(), True)
        else:
            make_subplot(image, "None", axes[i,0], ycens[i], xcens[i], image.min(), image.max())
            make_subplot(psfimage, "None", axes[i,1], ycens[i], xcens[i], image.min(), image.max())
            make_subplot(residuals, "None", axes[i,2], ycens[i], xcens[i], image.min(), image.max(), True)
   
    fig.tight_layout()
    fig.savefig(runprops.get("resultspath") + "/compare_plots_zoomed.png", dpi = 300)
    plt.close()
    
    
# Ensure that its being run from the right place
if 'results' in os.getcwd():
    runprops = ReadJson('runprops.txt').outProps()
else:
    print("remake_plots needs to be run from a results directory.")
    sys.exit()
    
# Load in chain file and image/psfs arrays
resultspath = os.getcwd()
runprops["resultspath"] = resultspath
backend = emcee.backends.HDFBackend(resultspath + '/chain.h5')

#os.chdir("../../src")
os.chdir("../../../src")

# Rerun plots
from analysis import *
nwalkers = runprops.get("nwalkers")
nsteps = runprops.get("nsteps")
nburnin = runprops.get("nburnin")
output_filename = runprops.get("output_filename")

runprops["best_likelihood"] = -np.inf

# Get initial guess (with init_guess.py)
# This uses test data from guess_test.py
params_df = pd.read_csv(resultspath + "/startguess.csv",sep=',',index_col=0)

p0_df = init_guess(params_df, nwalkers)
paramnames = params_df.columns.tolist()
p0, change_dict = params_to_fitarray(p0_df)

# Loading in image to be solved and making a small postage stamp version
x = runprops.get("stamp_x")
y = runprops.get("stamp_y")
size = runprops.get("stamp_size")

f = runprops.get('image_path') + runprops.get('input_image')
imageraw, filter, nchip, bunits = getimage_hst(f)

if bunits == 'ELECTRONS':
    print("bunits:", bunits)
else:
    print("Image header units must be in electrons for accurate results. Quitting.")
    sys.exit()

print("filter:", filter)
print("CCD chip:", nchip)

# Clean cosmic rays from image (maybe this can be removed when guesses are good enough?)
# This may also be irrelevant if we move to simultaneous pair fitting
import ccdproc
if runprops.get('cr_clean')==True:
    print("cr_clean set to True")
    cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim"), gain_apply = False)
    image = np.array(cr_cleaned)[x:x+size,y:y+size]

    # Make sure that the CR algorithm hasn't marked the target as a cosmic ray!
    reset = False
    maskcheck = runprops.get("maskcheck")

    while maskcheck:
        crredo = False
        for i in range(nwalkers):
            #print(i,crmask[int(p0[i,0]),int(p0[i,1])])
            if crmask[x:x+size,y:y+size][int(p0[i,0]),int(p0[i,1])]:
                crredo = True
        if crredo:
            cr_cleaned, crmask = ccdproc.cosmicray_lacosmic(imageraw, sigclip=5.0, objlim = runprops.get("cr_objlim") + 1.5, gain_apply = False)

            if reset:
               print("CR rejection algorithm is flagging your target as a CR. Aborting run.") 
               print("Consider increasing the objlim in runprops.")
               sys.exit()
            reset = True
        else:
            break
    image = np.array(cr_cleaned)[x:x+size,y:y+size]

# Print out image that hasn't been cleaned
else:
    print("cr_clean set to False")
    image = imageraw[x:x+size,y:y+size]
    
    plt.figure()
    plt.imshow(image, cmap = "hot", interpolation = "nearest", origin = "lower")
    plt.colorbar()
    plt.savefig(resultspath + "/cleanedimage.png")
    plt.close("all")

# generate the synthetic image
# Ensure there are no negative pixels. Add constant offset, which will be corrected in model images.
if np.nanmin(image) < 0:
    image = image - np.floor(np.nanmin(image)) + 1.0

# Calculating image noise characteristics for priors
runprops["med_noise"] = np.median(image)
runprops["std_noise"] = np.sqrt(runprops["med_noise"])
print(runprops["std_noise"])
print((runprops.get("std_noise")*runprops.get("noise_cutoff"))/0.1)

# Getting inputs for Tiny Tim PSFs
fmin = runprops.get("fmin")
fmax = runprops.get("fmax")
fspace = runprops.get("fspace")
focuses = np.arange(fmin, fmax + fspace, fspace)
xpos = runprops.get("xpos")
ypos = runprops.get("ypos")
size_psf = runprops.get("psf_size")
sample_factor = runprops.get("sample_factor")
#nchip = runprops.get("chip")
ndet = runprops.get("det_int")
#filter = runprops.get("filter")
numpsfs = focuses.size

# Making TinyTim PSFs (this is skipped if all the PSFs have previously been made)
for i,focus in enumerate(focuses):
    filename = "modelpsfs/wfc3psf_" + str(ndet) + "_" + str(nchip) + "_" + filter + "_" + str(xpos) + "_" + str(ypos) + "_" + str(round(focus,1)) + "_" + str(size_psf) + "_" + str(sample_factor) + ".fits"
    if round(focus,1) == 0.0:
        focus = 0.001
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

#import the sigsdf file
sigs_df = pd.read_csv(resultspath +'/sigsdf.csv',sep=',',index_col=0)
xsig1 = sigs_df.at['x1','1sigma']
xsig1minus = sigs_df.at['x1', '-1sigma']
ysig1 = sigs_df.at['y1','1sigma']

print(image.shape)

# residuals = image - psfimage
psfimage, residuals, xcens, ycens = plot_comparison(backend, image, psfs, focuses, runprops)
print(psfimage)
print(psfimage.shape)
xsig1 = np.full((60,60),xsig1)
print(xsig1)
print(xsig1.shape)
ysig1 = np.full((1,60),ysig1)
sigmas = np.array([ysig1,xsig1])
print(sigmas)
print(sigmas.shape)

#print(xsig1)
#sigmasy, sigmasx = np.meshgrid(ysig1, xsig1)
#print(sigmasy)
#print(sigmasy.shape)
#psfimage = psfimage + (ysig1, xsigma1)
#print(psfimage)

# scale the residuals
#uncertainties = 
#residuals = residuals / uncertainties

npsfs = runprops.get("npsfs")
#comparison(image, psfimage, residuals, xcens, ycens)

#zoomed_comparison(image, psfimage, residuals, xcens, ycens, npsfs)

