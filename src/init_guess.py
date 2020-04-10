# init_guess.py
# part of Physics 529 nPSF
# created by Ian Clark
# April 10, 2020

# Information from Dr. Ragozzine's email about this function:

# Initial guesses can come from a range of places.  Dr. Ragozzine
# will distinguish the initial guess (what goes into emcee) and
# the starting guess (what the user comes up with).

# A typical starting guess might be a mean and standard deviation
# for every parameter. So, I might specify that the initial 
# guesses for xcen should be nwalkers draws from a Gaussian with
# center 50 and standard deviation 2. For now, you can assume that
# this starting guess information comes from a dataframe where the
# columns are the parameter names, the first row is the mean value 
# and the second row is the standard deviation. The output will be 
# a parameters dataframe with the same column names, nwalker rows 
# drawn from this distribution. As before, write a test function 
# to see that this works.  

import pandas as pd

def init_guess(start_guess_df):
    """This function will produce the initial guess used in nPSF.
    
    Input: 
    
    start_guess_df - A dataframe that has the columns as the
        parameter names, the first row is the mean value for each 
        parameter, and the second row is the standard deviation 
        for each parameter.
    
    Returns:
    
    params_df - A parameters dataframe with the same column names
        as start_guess_df and nwalker rows drawn from the 
        distribution.
    """
    
    
    return params_df