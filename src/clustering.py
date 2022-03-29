#
#	clustering.py
#
#	Implements the clustering algorithm described in Hou 2010
#
#	Benjamin Proudfoot
#	01/20/22
#

# load in all required packages and functions
import numpy as np
import pandas as pd
import random
import commentjson as json
import emcee

def clustering(sampler, state, paramnames, log_likelihood, backend, image, psfs, focuses, runprops, const = 50, lag = 10, max_prune_frac = 0.9):
	nwalkers = runprops.get("nwalkers")
	reburnin = runprops.get("clustering_burnin")
	if reburnin == 0:
		return sampler, state

	# Getting important values from the chain
	avllhood = np.mean(sampler.get_log_prob()[-lag:,:], axis = 0)
	lastparams = sampler.get_chain()[-1,:,:]
	ngens = sampler.get_chain().shape[0]

	if ngens < lag:
		print("Chain too short for clustering algorithm, clustering not performed")
		return sampler, state

	# Sorting the walkers by likelihood values
	llhoodparam = pd.DataFrame(columns = ['llhood'] + paramnames)
	for i in range(nwalkers):
		llhoodparam.loc[i] = np.concatenate([np.array([avllhood[i]]),lastparams[i,:]])
	llhoodparam = llhoodparam.sort_values(by=['llhood'], ascending = False)
	llhoodparam = llhoodparam.values

	# Performing rejection tests
	reject = np.zeros(nwalkers)
	for i in range(1,nwalkers-1):
		term1 = -llhoodparam[i+1,0] + llhoodparam[i,0]
		term2 = const*(-llhoodparam[i,0] + llhoodparam[0,0])/(i)
		print(term1, term2)
		if term1 > term2:
			reject[(i+1):] = 1
			break
	freject = reject.sum()/nwalkers
	print(freject)
	ndim = np.shape(lastparams)[1]
	# Pruning walkers based on the clusters found,
	# replacing them with random linear combinations of walkers within the cluster
	# Skipping if cluster is not big enough
	if freject < max_prune_frac:
		params = llhoodparam[:,1:]
		for i in range(len(reject)):
			if reject[i] == 1:
				p = random.random()
				c1 = random.randrange(i)
				c2 = random.randrange(i)
				while c1 == c2:
					c2 = random.randrange(i)
				params[i,:] = (p*params[c1,:] + (1-p)*params[c2,:])
		#sampler.reset()
		sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood, backend=backend,
						args = [image, psfs, focuses, runprops])
		state = sampler.run_mcmc(params, reburnin, progress = True)
		return sampler, state
	else:
		print("Cluster not big enough, clustering not performed")
		return sampler, state
