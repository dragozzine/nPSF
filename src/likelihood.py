# likelihood.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020
#
# This script includes various functions for calculating likelihoods
#
# log_likelihood(parameter array) - returns the log of the prior probability for the given parameters
# LATER: log_likelihood_df(parameters dataframe, priors dataframe,observations dataframe) 
# - returns the model dataframe
# log_probability(parameter array) - returns the log posterior probability that emcee uses
# LATER: draw_prior - draws from the priors dataframe
# 
# test_log_likelihood - tests log_likelihood


# log_likelihood
# Inputs:
# - parameter vector [which parameters are which and go with which priors to be handled LATER]
# Outputs:
# - the (natural) log likelihood for these parameters
# including -np.inf for unallowed priors


# log_probability
# Inputs: 
# - parameter vector [which parameters are which and go with which priors to be handled LATER]
# Outputs:
# - the (natural) log posterior probability for these parameters
# 
# Example from emcee (https://emcee.readthedocs.io/en/stable/tutorials/line/):
#def log_probability(theta, x, y, yerr):
#    lp = log_prior(theta)
#    if not np.isfinite(lp):
#        return -np.inf
#    return lp + log_likelihood(theta, x, y, yerr)

def log_probability(parameters):
	lp = log_prior(paramaters)
	if not np.isfinite(lp):
		return 0.5
		#return -np.inf
	#return lp + log_likelihood(parameters, x, y, yerr)
	return 0.5

# test log_likelihood
# give parameters with a known answer and confirm it is correct

