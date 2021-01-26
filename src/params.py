# params.py
# part of Physics 529 nPSF
# created by Darin Ragozzine
# functions finished by Ian Clark
# April 9, 2020
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

#def params_to_fitarray(param_df, fixation_vals):
def params_to_fitarray(param_df):
    """This function takes a dataframe and converts it into an array
    that can be processed and later used in emcee.  It also creates
    a dictionary that holds the column names so that the dataframe
    can be recreated later in the fitarray_to_params function.
    
    Inputs:
    
    param_df - The parameters dataframe
        
    Returns:
    
    emcee_arr - An array that emcee can use.
    
    change_dict - A dictionary that tracks which parameter name was
        assigned to which value in the array. (Basically the names of the
        dataframe's columns and which index they match)
    """
        
    # Questions:
    #    -Just what will the fixation values array look like, if I use those?
    
    # Takes the dataframe and transforms it into a matrix that emcee can use.
    emcee_arr = param_df.to_numpy()
    
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
    array to a parameters dataframe.
    
    Inputs:
    
    emcee_arr - An array that emcee can use.  It's coming back here to
        be remade into a dataframe.
    
    change_dict - A dictionary that tracks which parameter name was
        assigned to which value in the array.
    
    Returns:
    
    param_df - The parameters dataframe
    """
    
    # We take the dictionary values and put them in a list form that
    #    the dataframe can use as its columns.
    dict_vals = []
    for x in change_dict:
        dict_vals.append(change_dict[x])
    
    # We make the dataframe
    param_df = pd.DataFrame(emcee_arr,
                            columns = [dict_vals])
    
    #return param_df, fixation_vals
    return param_df
