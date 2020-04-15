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
import numpy as np

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
    
    # The mean will be the center of our distribution.
    # The standard deviation will be the width.
    # With these ideas in mind, I can potentially figure out the 
    #    Gaussian distribution, but there might be a pre-existing
    #    function to do that.
    # So there is: np.random.normal(mu, sigma, nwalkers)
    
    # With these tools, we can just make the items for each column.
    # The only thing I'm wondering about is how to choose the number
    #    of walkers.  I'll choose 100 for now, until we can pick a
    #    better value.  As far as I know, they should all be the same
    #    number for the whole dataframe.
    
    # A guess for the number of walkers and taking the data in 
    #    array form.
    nwalkers = 100
    arrSet = start_guess_df.as_matrix()
    
    # Some code to help us get the names for the columns.
    name_dict = {0:"Junk"}
    n = 0
    dist_arr = []
    for col in start_guess_df.columns:
        name_dict[n] = col
        infos = start_guess_df[col].as_matrix()
        mean1, stdev1 = infos[0],infos[1]
        if n == 0:
            dist_arr = np.random.normal(mean1,stdev1,10)
        else:
            dist_arr = np.vstack((dist_arr,np.random.normal(mean1,stdev1,10)))
        
        n += 1
    
    dict_vals = []
    for x in name_dict:
        dict_vals.append(name_dict[x])
    
    params_df = pd.DataFrame(dist_arr, columns = [dict_vals])
    
    
    return params_df