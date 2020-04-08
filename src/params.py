# params.py
# part of Physics 529 nPSF
# created by Darin Ragozzine
# April 6, 2020
# 
# functions used to work with the parameters
# 
# paramsdf_to_fitarray - converts a parameters dataframe
# and information about which parameters are to be held fixed
# or floating and makes an array that emcee can use
# along with a dictionary that tracks which parameter name
# was assigned to which value in the array
#
# fitarray_to_paramsdf - does the reverse of the above

import pandas as pd

def params_to_fitarray(param_df, fixation_vals):
    """Converts a parameters dataframe and information about which
    parameters are to be held fixed or floating and makes an array
    that emcee can use along with a dictionary that tracks which
    parameter name was assigned to which value in the array.
    
    Parameters:
    
    param_df - The parameters dataframe
    
    fixation_vals - Information about which parameters are to be held
        fixed or floating.  But these may not be used.
        
    Returns:
    
    emcee_arr - An array that emcee can use.
    
    change_dict - A dictionary that tracks which parameter name was
        assigned to which value in the array. (Basically the names of the
        dataframe's columns and which index they match)
    """
    
    # To do:
    #    -Figure out how to automatically create the dictionary.
    
    # Notes:
    #    -So it sounds like emcee cannot use a dataframe, so this function
    #        will take the necessary values out of the dataframe and put them
    #        into an array.
    #    -It doesn't sound like I'm doing any data manipulation in here...
    #    -emcee uses a 2-D array.  The dictionary will pretty much hold the columns.
    
    # Questions:
    #    -Just what will the fixation values array look like, if I use those?
    
    # Takes the dataframe and transforms it into a matrix that emcee can use.
    emcee_arr = param_df.as_matrix()
    
    # Here we create the dictionary. We need to load it with a junk value to start.
    # We also include an integer for the purposes of creating the dictionary keys.
    change_dict = {0:"Junk"}
    n = 0
    
    for col in param_df.columns:
        change_dict[n] = col
        n += 1
    
    return emcee_arr, change_dict

def fitarray_to_paramsdf(emcee_arr, change_dict):
    """Converts the array that emcee uses and the dictionary that
    tracks which parameter name was assigned to which value in the
    array to a parameters dataframe and information about which
    parameters are to be held fixed or floating.
    
    Parameters:
    
    emcee_arr - An array that emcee can use.  It's coming back here to
        be remade into a dataframe.
    
    change_dict - A dictionary that tracks which parameter name was
        assigned to which value in the array.
    
    Returns:
    
    param_df - The parameters dataframe
    
    fixation_vals - Information about which parameters are to be held
        fixed or floating.
    """
    param_df = pd.DataFrame(emcee_arr,
                            colums = ["What the columns are"])
    
    return param_df, fixation_vals