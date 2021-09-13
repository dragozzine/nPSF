#	analysis.py
#
#	Runs all of the analysis on the chains from npsf_run
#
#	Benjamin Proudfoot
#	05/05/20
#

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import emcee
from statsmodels.regression.linear_model import yule_walker


def plots(sampler, parameters):
			# Here parameters is whatever file/object will have the run params
	flatchain = sampler.get_chain(flat = True)
	chain = sampler.get_chain(flat = False)
	# First start by converting the paramaters into an array of strings
	# code here
	names = list(paramaters)

	fig = corner.corner(chain, bins = 40, labels = names, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".4f")
	fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	for ax in fig.get_axes():
		ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	#fig.savefig(place to save the corner plot)

	
	# Now make the walker plots
	numparams = chain.shape[2]
	numwalkers = chain.shape[1]
	numgens = chain.shape[0]
	for i in range(numparams):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[0:numgens,j,i], numgens))
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		plt.savefig
		plt.close()

	# Likelihood plots??
	llhoods = sampler.get_log_prob(flat = True)
	for i in range(numparams):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(chain[i,:,:].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(flatchain[i,:].flatten(), llhoods.flatten())
		plt.xlabel(names[i])
		plt.ylabel("Log(L)")
		plt.subplot(224)
		plt.hist(llhoods.flatten(), bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		#plt.savefig(#place to save this)
		plt.close("all")
	#end

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
    The Geweke score on some series x is computed by:
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

