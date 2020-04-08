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
    #    -So it sounds like emcee cannot use a dataframe, so this function
    #        will take the necessary values out of the dataframe and put them
    #        into an array.
    #    -It doesn't sound like I'm doing any data manipulation in here...
    
    # Questions:
    #    -What kind of array does emcee use that this function will need to give it?
    #        (e.g., what values will we want?)
    #    -For clarity: The change dictionary will only be loaded with things that
    #        the fixation values array will say are floating? So parameter names and
    #        indexes?
    #    -What all will be in the parameters dataframe?
    #    -Just what will the fixation values array look like?
    
    emcee_arr = param_df.as_matrix()
    
    # Dictionary assignment will look something like this:
    # change_dict["param"] = array_val (probably an index)
    
    return emcee_arr, change_dict

def fitarray_to_paramsdf(emcee_arr, change_dict):
    """Converts the array that emcee uses and the dictionary that
    tracks which parameter name was assigned to which value in the
    array to a parameters dataframe and information about which
    parameters are to be held fixed or floating.
    
    Parameters:
    
    emcee_arr - An array that emcee can use.
    
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