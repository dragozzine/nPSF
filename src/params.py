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
        fixed or floating.
        
    Returns:
    
    emcee_arr - An array that emcee can use.
    
    change_dict - A dictionary that tracks which parameter name was
        assigned to which value in the array.
    """
    
    # To do:
    #    -Find out just what these values are that will be input.
    #    -Find out what kind of array emcee uses.
    
    # Notes:
    #    -I'm guessing that someone else will be handling the fixation_vals,
    #        so I don't need to worry about making them, just using them.
    #    -So it sounds like emcee cannot use a dataframe, so this function
    #        will take the necessary values out of the dataframe and put them
    #        into an array.
    #    -It doesn't sound like I'm doing any data manipulation in here...
    
    # Questions:
    #    -What kind of array does emcee use that this function will need to give it?
    #    -For clarity: The change dictionary will only be loaded with things that
    #        the fixation values array will say are floating?  Is it a string array?
    #        Or maybe it's a dictionary in and of itself?
    #    -What all will be in the parameters dataframe?
    
    
    
    # Dictionary assignment will look something like this:
    # change_dict["param"] = array_val
    
    return emcee_arr, change_dict