#	analysis.py
#
#	Runs all of the analysis on the chains from npsf_run
#
#	Benjamin Proudfoot
#	05/05/20
#
#	updated ability to iterate through n psfs
#	William Giforos
#	05/31/22

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import corner
import emcee
import astropy.io
import astropy.wcs
from likelihood import log_likelihood
import pandas as pd
import os.path
from latlon_transform import convert_to_primary_centric

def plots(sampler, resultspath, runprops):
	# Load in info about run and open the image's WCS
	npsfs = runprops.get("npsfs")  
	name_dict = runprops.get("names_dict")
	objectnames = []
	for i in name_dict.values():
		objectnames.append(i)
	print(objectnames)
	f = astropy.io.fits.open(runprops.get('input_image'))
	w = astropy.wcs.WCS(f[2].header)
    
	# Getting the stored chain
	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))  
	chain = sampler.get_chain(discard=int(burnin+clusterburn), flat = False)
	flatchain = sampler.get_chain(discard=int(burnin+clusterburn), flat = True)
	llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn), flat = True)

    
	# Make derived parameters according to number of psfs
	if npsfs == 1:
		#pass	# No derived parameters??
		names = np.array(["x1","y1","h1","f"])
		make_corner_plot(flatchain, names, resultspath + "/corner.pdf")
		make_walker_plots(sampler, chain, names, resultspath, runprops)
		make_likelihood_plots(flatchain, llhoods, names, resultspath, runprops)

	elif npsfs >= 2:
		# Calculating derived parameters
		ralist = []
		for i in range(npsfs):
			ralist = np.append(ralist,'ra' + str(i + 1))        
		declist = []
		for i in range(npsfs):
			declist = np.append(declist,'dec' + str(i + 1))
            
		astromlist = np.append(ralist,declist)       
		ind = list(range(0,len(flatchain)))        
		astromdf = pd.DataFrame(columns = astromlist, index = ind)

		num = 0        
		for i in range(npsfs):
			astromdf[ralist[i]],astromdf[declist[i]] = w.pixel_to_world_values(flatchain[:,num].flatten() + runprops.get("stamp_x"), flatchain[:,(num+1)].flatten() + runprops.get("stamp_y"))
			num += 3
        
		num = npsfs - 2
		num1 = 5
		for i in range(npsfs-1):        
			astromdf['dra('+ str(npsfs - num) + '-1)'] = (astromdf[ralist[i+1]] - astromdf[ralist[0]])*3600*np.cos(np.deg2rad(astromdf[declist[0]]))
			astromdf['ddec('+ str(npsfs - num) + '-1)'] = (astromdf[declist[i+1]] - astromdf[declist[0]])*3600
			astromdf['dmag('+ str(npsfs - num) + '-1)'] = -2.5*np.log10(flatchain[:,num1].flatten()/flatchain[:,2].flatten())   
			astromdf['sep('+ str(npsfs - num) + '-1)'] = np.sqrt(astromdf['dra('+ str(npsfs - num) + '-1)']**2 + astromdf['ddec('+ str(npsfs - num) + '-1)']**2)        
			astromdf['pa('+ str(npsfs - num) + '-1)'] = np.arctan2(astromdf['ddec('+ str(npsfs - num) + '-1)'],astromdf['dra('+ str(npsfs - num) + '-1)'])*57.2958        
			astromdf['dx('+ str(npsfs - num) + '-1)'] = flatchain[:,(num1-2)].flatten() - flatchain[:,0].flatten()        
			astromdf['dy('+ str(npsfs - num) + '-1)'] = flatchain[:,(num1-1)].flatten() - flatchain[:,1].flatten()        
			num -= 1
			num1 += 3

		# Loading derived parameters into arrays
		names = []
		for i in range(npsfs):
			names = np.append(names,"x" + str(i + 1))
			names = np.append(names,"y" + str(i + 1))
			names = np.append(names,"h" + str(i + 1))
			if i == (npsfs - 1):
				names = np.append(names,'f')

		dnames = names.copy() 
		num = 1
		for i in range(npsfs-1):            
			dnames = np.append(dnames,"dra"+ str(npsfs - num))
			dnames = np.append(dnames,"ddec"+ str(npsfs - num))
			dnames = np.append(dnames,"dmag"+ str(npsfs - num))
			dnames = np.append(dnames,"sep"+ str(npsfs - num))
			dnames = np.append(dnames,"pa"+ str(npsfs - num))
			dnames = np.append(dnames,"dx"+ str(npsfs - num))
			dnames = np.append(dnames,"dy"+ str(npsfs - num))
			num -= 1
       
		dfchain = flatchain.copy()  
		num = npsfs - 2
		for i in range(npsfs-1):        
			dfchain = np.concatenate((dfchain,np.array(astromdf['dra('+ str(npsfs - num) + '-1)']).reshape((astromdf['dra('+ str(npsfs - num) + '-1)'].size,1)) ), axis = 1)
			dfchain = np.concatenate((dfchain,np.array(astromdf['ddec('+ str(npsfs - num) + '-1)']).reshape((astromdf['ddec('+ str(npsfs - num) + '-1)'].size,1)) ), axis = 1)
			dfchain = np.concatenate((dfchain,np.array(astromdf['dmag('+ str(npsfs - num) + '-1)']).reshape((astromdf['dmag('+ str(npsfs - num) + '-1)'].size,1)) ), axis = 1)
			dfchain = np.concatenate((dfchain,np.array(astromdf['sep('+ str(npsfs - num) + '-1)']).reshape((astromdf['sep('+ str(npsfs - num) + '-1)'].size,1)) ), axis = 1)
			dfchain = np.concatenate((dfchain,np.array(astromdf['pa('+ str(npsfs - num) + '-1)']).reshape((astromdf['pa('+ str(npsfs - num) + '-1)'].size,1)) ), axis = 1)
			dfchain = np.concatenate((dfchain,np.array(astromdf['dx('+ str(npsfs - num) + '-1)']).reshape((astromdf['dx('+ str(npsfs - num) + '-1)'].size,1)) ), axis = 1)
			dfchain = np.concatenate((dfchain,np.array(astromdf['dy('+ str(npsfs - num) + '-1)']).reshape((astromdf['dy('+ str(npsfs - num) + '-1)'].size,1)) ), axis = 1)            
			num -= 1

		# Making plots
		make_corner_plot(flatchain, names, resultspath + "/corner.pdf")
		make_corner_plot(dfchain, dnames, resultspath + "/cornerderived.pdf")
		make_walker_plots(sampler, chain, names, resultspath, runprops)
		make_likelihood_plots(dfchain, llhoods, dnames, resultspath, runprops)
		print("Making the obsdf file")
		obsdf = make_obsdf(dfchain, astromdf, f, objectnames, npsfs, resultspath)
		print("Converting RA/DEC to LAT/LONG")
		forsigsdf = convert_to_primary_centric(obsdf, objectnames, npsfs, resultspath, 1000)
		print("Making the sigsdf file")
		make_sigsdf(dfchain, forsigsdf, llhoods, dnames, npsfs, resultspath)

	else:
		print("Error, enter a whole number of objects, 1 or greater. Aborting analysis.")
        
        
def likelihood_map(flatchain, llhoods, resultspath, runprops):
	# Set up for making likelihood maps
	# Begin by getting best fit position. Assume this is the median location of parameters
        # Load in info about run and open the image's WCS
	npsfs = runprops.get("npsfs")
	f = astropy.io.fits.open(runprops.get('input_image'))
	w = astropy.wcs.WCS(f[2].header)

	# Calculating derived parameters
	ra1,dec1 = w.pixel_to_world_values(flatchain[:,0].flatten() + runprops.get("stamp_x"), flatchain[:,1].flatten() + runprops.get("stamp_y"))
	ra2,dec2 = w.pixel_to_world_values(flatchain[:,3].flatten() + runprops.get("stamp_x"), flatchain[:,4].flatten() + runprops.get("stamp_y"))
	dra = (ra2 - ra1)*3600*np.cos(np.deg2rad(dec1))
	ddec = (dec2 - dec1)*3600
	dmag = -2.5*np.log10(flatchain[:,5].flatten()/flatchain[:,2].flatten())
	sep = np.sqrt(dra**2 + ddec**2)
	pa = np.arctan2(ddec,dra)*57.2958
	dx = flatchain[:,3].flatten() - flatchain[:,0].flatten()
	dy = flatchain[:,4].flatten() - flatchain[:,1].flatten()

	# Loading derived parameters into arrays
	names = np.array(["x1","y1","h1","x2","y2","h2","f"])
	dnames = names.copy()
	dnames = np.append(dnames, ["dra","ddec","dmag","sep","pa","dx","dy"])
	dfchain = flatchain.copy()
	dfchain = np.concatenate((dfchain,np.array(dra).reshape((dra.size,1)) ), axis = 1)
	dfchain = np.concatenate((dfchain,np.array(ddec).reshape((ddec.size,1)) ), axis = 1)
	dfchain = np.concatenate((dfchain,np.array(dmag).reshape((dmag.size,1)) ), axis = 1)
	dfchain = np.concatenate((dfchain,np.array(sep).reshape((sep.size,1)) ), axis = 1)
	dfchain = np.concatenate((dfchain,np.array(pa).reshape((pa.size,1)) ), axis = 1)
	dfchain = np.concatenate((dfchain,np.array(dx).reshape((dx.size,1)) ), axis = 1)
	dfchain = np.concatenate((dfchain,np.array(dy).reshape((dy.size,1)) ), axis = 1)

	# Filtering the samples
	bool_arr = llhoods > llhoods.max() - 10.0

	# Make brightness separation plot
	plt.figure()
	plt.scatter(sep[bool_arr], dmag[bool_arr], c = llhoods[bool_arr], cmap = "viridis", s = 5, alpha = 0.4)
	plt.xlabel("Separation (arcsec)")
	plt.ylabel("delta mag")
	plt.gca().invert_yaxis()
	color_bar = plt.colorbar()
	color_bar.set_alpha(1)
	plt.savefig(resultspath + "/brightness_sep.png", dpi = 300)
	plt.close()




def plot_best_fit(sampler, image, psfs, focuses, runprops):
	# Get the best parameter set from the chain
	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))  
	llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn), flat = True)
	flatchain = sampler.get_chain(discard=int(burnin+clusterburn), flat = True)
	ind = np.argmax(llhoods)
	params = flatchain[ind,:].flatten()

	# Run likelihhod function with plotit set to true
	log_likelihood(params, image, psfs, focuses, runprops, plotit = True)
	return params


def auto_window(taus, c):
	m = np.arange(len(taus)) < c * taus
	if np.any(m):
		return np.argmin(m)
	return len(taus) - 1


def autocorr_new(y, c = 5.0):
	f = np.zeros(y.shape[1])
	for yy in y:
		f += emcee.autocorr.function_1d(yy)
	f /= len(y)
	taus = 2.0 * np.cumsum(f) - 1.0
	window = auto_window(taus, c)
	return taus[window]

def autocorrelation(sampler, runprops):
	# Getting chain for the first parameter to calculate different values
	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))  
	chain = sampler.get_chain(discard=int(burnin+clusterburn), flat = False)
	
	#nwalkers = sampler.nwalkers
	nwalkers = runprops.get("nwalkers")
	ndims = sampler.ndim
	nsteps = chain.shape[0]
	
	# Converting parameters df to array of names
	names = np.array(["x1","y1","h1","x2","y2","h2"])	# placeholder for now

	# Calculating values to calculate tau for
	# This chould be changed eventually
	N = np.exp(np.linspace(np.log(100), np.log(nsteps), 10)).astype(int)

	# Setting up array for tau estimates
	tau = np.empty( (len(N), ndims) )

	# Calculating tau for each value in N for each parameter
	for i in range(ndims):
		thischain = chain[:,:,i].T
		for j, n in enumerate(N):
			tau[j,i] = autocorr_new(thischain[:, :n])

	# Setting up to plot curves for tau in a grid
	x = 3
	y = ndims
	nrows = 0
	ncols = 3
	while x <= y:
		y = y - x
		nrows += 1

	# Plotting
	fig, ax = plt.subplots(nrows = nrows, ncols = ncols, sharey=True, 
			       gridspec_kw={'wspace': 0},
			       figsize = (6.4*(ncols-1),4.8*(nrows)), 
			       squeeze = False)
	fig.suptitle("Autocorrelation estimates")
	fig.add_subplot(111, frameon=False)
	plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	plt.grid(False)
	plt.xlabel("number of samples, $N$")
	plt.ylabel(r"$\tau$ estimates")
	for i in range(nrows):
		for j in range(ncols):
			dim = i*ncols + j
			taus = ax[i,j].loglog(N, tau[:,dim], "o-", label="new")
			line = ax[i,j].plot(N, N / 50.0, "--k", label=r"$\tau = N/50$")
			ax[i,j].text(N[0], 100, names[dim])
			if i == 0 and j == 0:
				 plt.legend([line],[r"$\tau = N/50$"],fontsize = 14,
					    loc="lower right")

	fig.savefig("../results/autocorr.png")
	# Outputting the emcee calculated autocorrelation time as an additional check
	print(
		"Mean autocorrelation time: {0:.3f} steps".format(
			np.mean(sampler.get_autocorr_time(discard=int(burnin+clusterburn)))
		)
	)
	print("Estimated autocorrelation time for each parameter:")
	print(sampler.get_autocorr_time(discard=int(burnin+clusterburn)))

def spec(x, order=2):
    beta, sigma = yule_walker(x, order)
    return sigma**2 / (1. - np.sum(beta))**2


def geweke(x, first=.1, last=.5, intervals=20, maxlag=20):
    """Return z-scores for convergence diagnostics.
    Compare the mean of the first % of series with the mean of the last % of
    series. x is divided into a number of segments for which this difference is
    computed. If the series is converged, this score should oscillate between
    -1 and 1.
    Parameters
    ----------
    x : array-like
      The trace of some stochastic parameter.
    first : float
      The fraction of series at the beginning of the trace.
    last : float
      The fraction of series at the end to be compared with the section
      at the beginning.
    intervals : int
      The number of segments.
    maxlag : int
      Maximum autocorrelation lag for estimation of spectral variance
    Returns
    -------
    scores : list [[]]
      Return a list of [i, score], where i is the starting index for each
      interval and score the Geweke score on the interval.
    Notes
    -----
    The Geweke score on some series x is computed by:make_likelihood_plots(dfchain, llhoods, dnames, resultspath)
      .. math:: \frac{E[x_s] - E[x_e]}{\sqrt{V[x_s] + V[x_e]}}
    where :math:`E` stands for the mean, :math:`V` the variance,
    :math:`x_s` a section at the start of the series and
    :math:`x_e` a section at the end of the series.
    References
    ----------
    Geweke (1992)
    """

    if np.ndim(x) > 1:
        return [geweke(y, first, last, intervals) for y in np.transpose(x)]

    # Filter out invalid intervals
    if first + last >= 1:
        raise ValueError(
            "Invalid intervals for Geweke convergence analysis",
            (first, last))

    # Initialize list of z-scores
    zscores = [None] * intervals

    # Starting points for calculations
    starts = np.linspace(0, int(len(x)*(1.-last)), intervals).astype(int)

    # Loop over start indices
    for i,s in enumerate(starts):

        # Size of remaining array
        x_trunc = x[s:]
        n = len(x_trunc)

        # Calculate slices
        first_slice = x_trunc[:int(first * n)]
        last_slice = x_trunc[int(last * n):]

        z = (first_slice.mean() - last_slice.mean())
        z /= np.sqrt(spec(first_slice)/len(first_slice) +
                     spec(last_slice)/len(last_slice))
        zscores[i] = len(x) - n, z

    return zscores

def testconvergence_geweke(sampler, runprops):
	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))  
	#nwalkers = sampler.nwalkers
	nwalkers = runprops.get("nwalkers")
	flatchain = sampler.get_chain(discard=int(burnin+clusterburn), flat = True)
	names = np.array(["x1","y1","h1","x2","y2","h2"])
	for i in range(flatchain.shape[1]):
		zscores = geweke(np.reshape(flatchain[:,i], flatchain.shape[0]), first = 0.1,
				last = 0.5, intervals = 20)
		means = np.mean(sampler.get_chain(discard=int(burnin+clusterburn))[:,:,i].T, axis = 0)
		zscoresmean = geweke(means.flatten(), first = 0.1, last = 0.5, intervals = 20)
		plt.figure()
		plt.plot(np.array(zscores).T[0]/nwalkers, np.array(zscores).T[1], "o-")
		plt.plot(np.array(zscoresmean).T[0], np.array(zscoresmean).T[1], "o-")
		plt.hlines([-1,1], -100, np.array(zscores).T[0].max()/nwalkers + 100, linestyles = "dashed")
		plt.xlabel("First iteration")
		plt.ylabel("Geweke scores")
		plt.xlim((0,np.array(zscores).T[0].max()/nwalkers + 100))
		plt.savefig("../results/" + names[i] + "_geweke.png")
		plt.close()

def make_corner_plot(flatchain, names, filename):
	""" Makes a corner plot from the data from sampler. I may continue to
	 make adjustments to this function.

	Input: flatchain (2D version of data from sampler)

	Output: Saves an image of the corner plot. Returns nothing.
	"""
	fig = corner.corner(flatchain, labels = names, plot_datapoints = False, color = "blue", 
			    fill_contours = True, show_titles = True, bins = 40, title_fmt = ".6f")
	fig.savefig(filename)

def make_walker_plots(sampler, chain, names, path, runprops):
	if not os.path.exists(path + "/walkers_pngs"):
		os.makedirs(path + "/walkers_pngs")
	numparams = chain.shape[2]
	numwalkers = chain.shape[1]
	numgens = chain.shape[0]
	for i in range(numparams):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[0:numgens,j,i], numgens))
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		plt.savefig(path + "/walkers_pngs/" + names[i] + "_walkers.png")
		plt.close()

	# Now make the walker plot PDFs
	from matplotlib.backends.backend_pdf import PdfPages

	walkerpdf = PdfPages(path + "/walkers.pdf")
	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))  
	shortchain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False)
	print('Shortchain:', shortchain.shape)    
	#shortchain = sampler.get_chain(flat = False)
	numparams = shortchain.shape[2]
	numwalkers = shortchain.shape[1]
	numgens = shortchain.shape[0]
	del shortchain
	for i in range(numparams):
		plt.figure(dpi = 50)
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[0:numgens,j,i], numgens), alpha=0.2)
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		#plt.savefig(runprops.get('results_folder')+"/walker_"+names[i]+".png")
		walkerpdf.attach_note(names[i])
		walkerpdf.savefig()
		#plt.close()

	walkerpdf.close()
	plt.close("all")
    
	fullwalkerpdf = PdfPages(path + "/walkers_full.pdf")
	backend = emcee.backends.HDFBackend('chain.h5')    

	full_chain = sampler.get_chain(discard=0, flat = False)  
	fullgens = full_chain.shape[0]
	#print(thin_plots, fullgens, full_chain.shape)
	runnedthin = 1    
	thin_plots = 1  
	##if runprops.get('thin_run'):
		##runnedthin = runprops.get('nthinning') 
	for i in range(numparams):
		plt.figure(dpi = 50)
		for j in range(numwalkers):
			plt.plot(np.reshape(full_chain[0:fullgens,j,i], fullgens), alpha=0.2)
		plt.axvline(x=int(burnin/thin_plots/runnedthin))
		plt.axvline(x=int(clusterburn/thin_plots/runnedthin+burnin/thin_plots/runnedthin))
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		#plt.savefig(runprops.get('results_folder')+"/walker_"+names[i]+".png")
		fullwalkerpdf.attach_note(names[i])
		fullwalkerpdf.savefig()
		#plt.close()

	fullwalkerpdf.close()
	plt.close("all")        

def make_likelihood_plots(dfchain, llhoods, dnames, resultspath, runprops):
	nwalkers = runprops.get("nwalkers")
	likelihoodspdf = PdfPages(resultspath + "/likelihoods.pdf")
	ylimmin = np.percentile(llhoods.flatten(), 1)
	ylimmax = llhoods.flatten().max() + 1
	for i in range(dnames.size):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(dfchain[:,i].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(dfchain[:,i].flatten(), llhoods.flatten(),
			    c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), nwalkers),
			    cmap = "nipy_spectral", edgecolors = "none", rasterized=True, alpha=0.1)
		plt.xlabel(dnames[i])
		plt.ylabel("Log(L)")
		#plt.xlim(-0.25, 0.25)
		plt.ylim(ylimmin, ylimmax)
		plt.subplot(224)
		llflat = llhoods.flatten()
		plt.hist(llflat[np.isfinite(llflat)], bins = 40, orientation = "horizontal",
			 histtype = "step", color = "black")
		plt.ylim(ylimmin, ylimmax)
		likelihoodspdf.savefig()
	likelihoodspdf.close()
	plt.close("all")

def make_bright_sep(sep, dmag, resultspath, weights = None):
	if weights == None:
		weights = np.ones(sep.size)

	sep = np.array(sep)
	dmag = np.array(dmag)

	colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']

	num_bins = 20
	log_bins = np.logspace(-2,0,num_bins,endpoint=True)
	bin_indices_log = np.digitize(sep,log_bins)

	sig1vals = []
	sig2vals = []
	sig3vals = []
	dmag_bin = np.transpose([dmag,bin_indices_log])
	weights_bin = np.transpose([weights,bin_indices_log])

	from statsmodels.stats.weightstats import DescrStatsW

	for i in range(1,num_bins+1):
		bin1 = np.extract(dmag_bin[:,1]==i,dmag_bin[:,0])
		weight1 = np.extract(weights_bin[:,1]==i,weights_bin[:,0])
		if len(bin1)!=0:
			sigs = DescrStatsW(data = bin1,weights=weight1).quantile(probs = [0.3173,0.0455,0.0027],return_pandas=False)
			sig1vals.append(sigs[0])
			sig2vals.append(sigs[1])
			sig3vals.append(sigs[2])
		else:
			sig1vals.append(np.nan)
			sig2vals.append(np.nan)
			sig3vals.append(np.nan)

	fig = plt.figure(figsize=(8,6))
	ax1 = fig.add_subplot(111)
	ax1.plot(log_bins,sig3vals,color = colorcycle[0], marker = '.', markersize = 10, label = r'3$\sigma$')
	ax1.plot(log_bins,sig2vals,color = colorcycle[1], marker = '.', markersize = 10, label = r'2$\sigma$')
	ax1.plot(log_bins,sig1vals,color = colorcycle[2], marker = '.', markersize = 10, label = r'1$\sigma$')
	ax1.set_xlabel("Separation (arcsec)")
	ax1.set_ylabel(r"$\Delta$ mag")
	ax1.set_xscale("log")

	#plt.xlabel("Separation (arcsec)")
	#plt.ylabel("delta mag")
	plt.gca().invert_yaxis()
	#plt.xscale('log')
	ax1.legend(loc='upper right')
	plt.savefig(resultspath + "/bright_sep_lim.pdf")
	plt.close()

    
def make_obsdf(dfchain, astromdf, f, objectnames, npsfs, resultspath):
	# Get time from image header. We want the midtime of exposure, so ((EXPEND - EXPSTART) / 2)
	# Where EXPEND and EXPSTART are in Modified Julian Date (MJD)
	expend = f[0].header['expend']
	expstart = f[0].header['expstart']
	timeMJD = expstart + ((expend - expstart) / 2)
	timeJD = timeMJD + 2400000.5
	print("Observation time JD:", timeJD)
	print("size of data frame:", astromdf.shape)
        
    #Make obsdf dataframe
	#obsdf = astromdf[0:10].copy()
	obsdf = astromdf.copy()
    
	obsdf.insert(0, 'time',pd.Series(["{:10.4f}".format(timeJD)], index=[0]))
	obsdf['time'] = obsdf['time'].astype('float64')    
	for i in range(npsfs):
		obsdf.drop(columns=['ra' + str(i+1),'dec' + str(i+1)], inplace=True)
            
	num = npsfs - 2                  
	for i in range(npsfs-1):
		obsdf.drop(columns=['dmag('+ str(npsfs - num) + '-1)','sep('+ str(npsfs - num) + '-1)','pa('+ str(npsfs - num) + '-1)','dx('+ str(npsfs - num) + '-1)','dy('+ str(npsfs - num) + '-1)'], inplace=True)
		obsdf.rename(columns={'dra('+ str(npsfs - num) + '-1)': 'Delta-RA_'+objectnames[i+1], 'ddec('+ str(npsfs - num) + '-1)': 'Delta-DEC_'+objectnames[i+1]}, inplace=True)
		num -= 1

    #Calculate errors for dra and ddec for each moon
	num = 0                      
	for k in range(npsfs-1):  
		num_dra = dfchain[0:len(obsdf),(3*npsfs)+1+num]
		num_ddec = dfchain[0:len(obsdf),(3*npsfs)+2+num]
            
		neg1sig_dra = np.percentile(num_dra,15.866, axis = None)
		pos1sig_dra = np.percentile(num_dra,84.134, axis = None)
		neg1sig_ddec = np.percentile(num_ddec,15.866, axis = None)
		pos1sig_ddec = np.percentile(num_ddec,84.134, axis = None)
		ralist = []     
		declist = []     
        
		for i in range(len(obsdf)):
         
			error_neg_dra = neg1sig_dra - obsdf['Delta-RA_'+objectnames[k+1]].iloc[i]
			error_pos_dra = pos1sig_dra - obsdf['Delta-RA_'+objectnames[k+1]].iloc[i]  
			error_temp_dra = (np.abs(error_neg_dra) + np.abs(error_pos_dra))/2
			error_neg_ddec = neg1sig_ddec - obsdf['Delta-DEC_'+objectnames[k+1]].iloc[i]
			error_pos_ddec = pos1sig_ddec - obsdf['Delta-DEC_'+objectnames[k+1]].iloc[i]  
			error_temp_ddec = (np.abs(error_neg_ddec) + np.abs(error_pos_ddec))/2
      
			ralist = np.append(ralist,error_temp_dra)        
			declist = np.append(declist,error_temp_ddec)        

		obsdf['Delta-RA_'+objectnames[k+1]+'-err'] = ralist
		obsdf['Delta-DEC_'+objectnames[k+1]+'-err'] = declist           
		num += 7
        
	return obsdf      
        
    
def make_sigsdf(dfchain, forsigsdf, llhoods, dnames, npsfs, resultspath):
	num = 1
	for i in range(npsfs-1):            
		dnames = np.append(dnames,"lat"+ str(npsfs - num))
		dnames = np.append(dnames,"long"+ str(npsfs - num))
		dnames = np.append(dnames,"dlat"+ str(npsfs - num))
		dnames = np.append(dnames,"dlong"+ str(npsfs - num))
		num -= 1

	#dfchain = dfchain[0:10,:]
	num = 1
	for i in range(npsfs-1):        
		dfchain = np.concatenate((dfchain,np.array(forsigsdf['lat'+ str(npsfs - num)]).reshape((forsigsdf['lat'+ str(npsfs - num)].size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(forsigsdf['long'+ str(npsfs - num)]).reshape((forsigsdf['long'+ str(npsfs - num)].size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(forsigsdf['dlat'+ str(npsfs - num)]).reshape((forsigsdf['dlat'+ str(npsfs - num)].size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(forsigsdf['dlong'+ str(npsfs - num)]).reshape((forsigsdf['dlong'+ str(npsfs - num)].size,1)) ), axis = 1)
		num -= 1        
        
        
	ind = np.argmax(llhoods)
	#print(dnames[(3*npsfs)+1])
	#print(dnames[(3*npsfs)+2])
    
	sigsdf = pd.DataFrame(columns = ['-3sigma','-2sigma','-1sigma','median','1sigma','2sigma','3sigma', 'mean', 'best fit'], index = dnames)
	j = 0
   
	for i in range(len(dfchain[0])):
		num = dfchain[:,i]

		median = np.percentile(num,50, axis = None)
		neg3sig= np.percentile(num,0.37, axis = None)
		neg2sig = np.percentile(num,2.275, axis = None)
		neg1sig = np.percentile(num,15.866, axis = None)
		pos1sig = np.percentile(num,84.134, axis = None)
		pos2sig = np.percentile(num,97.724, axis = None)
		pos3sig = np.percentile(num,99.63, axis = None)
		mean = np.mean(num)
		bestfit = dfchain[ind,:].flatten()[i]
		sigsdf['-3sigma'].iloc[i] = neg3sig-median
		sigsdf['-2sigma'].iloc[i] = neg2sig-median
		sigsdf['-1sigma'].iloc[i] = neg1sig-median
		sigsdf['median'].iloc[i] = median
		sigsdf['1sigma'].iloc[i] = pos1sig-median
		sigsdf['2sigma'].iloc[i] = pos2sig-median
		sigsdf['3sigma'].iloc[i] = pos3sig-median
		sigsdf['mean'].iloc[i] = mean
		sigsdf['best fit'].iloc[i] = bestfit
        
	print(sigsdf)
	filename = 'sigsdf.csv'    
	sigsdf.to_csv(resultspath + "/" + filename)
        