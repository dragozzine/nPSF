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

