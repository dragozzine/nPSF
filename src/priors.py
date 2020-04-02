# priors.py
# part of nPSF
# outline from Darin Ragozzine
# April 2, 2020
#
# This script includes various functions for calculating priors. 
#
# log_prior(parameter array) - returns the log of the prior probability for the given parameters
# LATER: log_prior_df(parameters dataframe, priors dataframe) - returns the log of the prior probability for the given parameters
# LATER: draw_prior - draws from the priors dataframe
# 
# test_log_prior - tests log_prior


# log_prior
# Inputs:
# - parameter vector [which parameters are which and go with which priors to be handled LATER]
# Outputs:
# - the (natural) log likelihood of the prior for these parameters
# including -np.inf for unallowed priors


# log_prior_df
# Inputs: 
# - parameters dataframe: a pandas dataframe in the "parameters" format that has the value of the parameters
# - priors dataframe: a pandas dataframe in the "priors" format format that explains what the prior probability distribution is
# More information on the formats to come later.
# Outputs: 
# - the (natural) log likelihood of the prior for these parameters
# including -np.inf for unallowed priors
# 
# Notes: 
# use the chosen probability distributions to calculate priors
# see code from Dallin Spencer



# test log_prior
# give parameters and a prior with a known answer and confirm it is correct

