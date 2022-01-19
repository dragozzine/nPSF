#	analysis.py
#
#	Runs all of the analysis on the chains from npsf_run
#
#	Benjamin Proudfoot
#	05/05/20
#

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import corner
import emcee
import astropy.io
import astropy.wcs

def plots(sampler, resultspath, runprops):
	# Load in info about run and open the image's WCS
	npsfs = runprops.get("npsfs")
	f = astropy.io.fits.open(runprops.get('input_image'))
	w = astropy.wcs.WCS(f[2].header)

	# Getting the stored chain
	chain = sampler.get_chain()
	flatchain = sampler.get_chain(flat = True)
	llhoods = sampler.get_log_prob(flat = True)

	# Make derived parameters according to number of psfs
	if npsfs == 1:
		#pass	# No derived parameters??
		names = np.array(["x1","y1","h1","f"])
		make_corner_plot(flatchain, names, resultspath + "/corner.pdf")
		make_walker_plots(chain, resultspath)
		make_likelihood_plots(flatchain, llhoods, names, resultspath, runprops)

	elif npsfs == 2:
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

		# Making plots
		make_corner_plot(flatchain, names, resultspath + "/corner.pdf")
		make_corner_plot(dfchain, dnames, resultspath + "/cornerderived.pdf")
		make_walker_plots(chain, resultspath)
		make_likelihood_plots(dfchain, llhoods, dnames, resultspath, runprops)

	elif npsfs == 3:
		# Calculating derived parameters
		ra1,dec1 = w.pixel_to_world_values(flatchain[:,0].flatten() + runprops.get("stamp_x"), flatchain[:,1].flatten() + runprops.get("stamp_y"))
		ra2,dec2 = w.pixel_to_world_values(flatchain[:,3].flatten() + runprops.get("stamp_x"), flatchain[:,4].flatten() + runprops.get("stamp_y"))
		ra3,dec3 = w.pixel_to_world_values(flatchain[:,6].flatten() + runprops.get("stamp_x"), flatchain[:,7].flatten() + runprops.get("stamp_y"))
		dra2 = (ra2 - ra1)*3600*np.cos(np.deg2rad(dec1))
		ddec2 = (dec2 - dec1)*3600
		dra3 = (ra3 - ra1)*3600*np.cos(np.deg2rad(dec1))
		ddec3 = (dec3 - dec1)*3600
		dmag2 = -2.5*np.log10(flatchain[:,5].flatten()/flatchain[:,2].flatten())
		dmag3 = -2.5*np.log10(flatchain[:,8].flatten()/flatchain[:,2].flatten())
		sep2 = np.sqrt(dra2**2 + ddec2**2)
		sep3 = np.sqrt(dra3**2 + ddec3**2)
		pa2 = np.arctan2(ddec2,dra2)*57.2958
		pa3 = np.arctan2(ddec3,dra3)*57.2958
		dx2 = flatchain[:,3].flatten() - flatchain[:,0].flatten()
		dy2 = flatchain[:,4].flatten() - flatchain[:,1].flatten()
		dx3 = flatchain[:,6].flatten() - flatchain[:,0].flatten()
		dy3 = flatchain[:,7].flatten() - flatchain[:,1].flatten()

		# Loading derived parameters into arrays
		names = np.array(["x1","y1","h1","x2","y2","h2","x3","y3","h3","f"])
		dnames = names.copy()
		dnames = np.append(dnames, ["dra2","ddec2","dmag2","sep2","pa2","dx2","dy2","dra3","ddec3","dmag3","sep3","pa3","dx3","dy3"])
		dfchain = flatchain.copy()
		dfchain = np.concatenate((dfchain,np.array(dra2).reshape((dra2.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(ddec2).reshape((ddec2.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(dmag2).reshape((dmag2.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(sep2).reshape((sep2.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(pa2).reshape((pa2.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(dx2).reshape((dx2.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(dy2).reshape((dy2.size,1)) ), axis = 1)

		dfchain = np.concatenate((dfchain,np.array(dra3).reshape((dra3.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(ddec3).reshape((ddec3.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(dmag3).reshape((dmag3.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(sep3).reshape((sep3.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(pa3).reshape((pa3.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(dx3).reshape((dx3.size,1)) ), axis = 1)
		dfchain = np.concatenate((dfchain,np.array(dy3).reshape((dy3.size,1)) ), axis = 1)

		# Making plots
		make_corner_plot(flatchain, names, resultspath + "/corner.pdf")
		make_corner_plot(dfchain, dnames, resultspath + "/cornerderived.pdf")
		make_walker_plots(chain, resultspath)
		make_likelihood_plots(dfchain, llhoods, dnames, resultspath, runprops)
	else:
		print("Only 1-3 PSFs are currently supported. Aborting analysis.")

	# Insert more analysis here?

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

def autocorrelation(sampler):
	# Getting chain for the first parameter to calculate different values
	chain = sampler.get_chain()
	
	nwalkers = sampler.nwalkers
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
			np.mean(sampler.get_autocorr_time())
		)
	)
	print("Estimated autocorrelation time for each parameter:")
	print(sampler.get_autocorr_time())

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

def testconvergence_geweke(sampler):
	nwalkers = sampler.nwalkers
	flatchain = sampler.get_chain(flat = True)
	names = np.array(["x1","y1","h1","x2","y2","h2"])
	for i in range(flatchain.shape[1]):
		zscores = geweke(np.reshape(flatchain[:,i], flatchain.shape[0]), first = 0.1,
				last = 0.5, intervals = 20)
		means = np.mean(sampler.get_chain()[:,:,i].T, axis = 0)
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

def make_walker_plots(chain, path):
	numparams = chain.shape[2]
	numwalkers = chain.shape[1]
	numgens = chain.shape[0]
	names = np.array(["x1","y1","h1","x2","y2","h2","f"])
	for i in range(numparams):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[0:numgens,j,i], numgens))
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		plt.savefig(path + "/" + names[i] + "_walkers.png")
		plt.close()

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
		plt.ylim(ylimmin, ylimmax)
		plt.subplot(224)
		llflat = llhoods.flatten()
		plt.hist(llflat[np.isfinite(llflat)], bins = 40, orientation = "horizontal",
			 histtype = "step", color = "black")
		plt.ylim(ylimmin, ylimmax)
		likelihoodspdf.savefig()
	likelihoodspdf.close()
	plt.close("all")
