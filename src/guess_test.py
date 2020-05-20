import sys
sys.path.insert(1,".")

# Import the libraries and functions we'll need.
import numpy as np
import pandas as pd

def make_test_df():
    # Create a dataframe of "averages" and "standard deviations" that init_guess
    # will use.
    foob_arr = [[50,50,2000,43,57,800],[1,1,500,1,1,200]]
    #freb_arr = np.array([[23,40,20,40,38,5],[4,2,1,0.4,8,0.5]])
    colies = ["x1","y1","h1","x2","y2","h2"]
    test_df = pd.DataFrame(foob_arr, columns = colies)
    return test_df
