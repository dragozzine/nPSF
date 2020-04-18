import sys
sys.path.insert(1,".")

# Import the libraries and functions we'll need.
import numpy as np
import pandas as pd

from init_guess import init_guess

def make_test_df():
    # Create a dataframe of "averages" and "standard deviations" that init_guess
    # will use.
    foob_arr = [[23,40,20,40,38,5],[4,2,1,0.4,8,0.5]]
    #freb_arr = np.array([[23,40,20,40,38,5],[4,2,1,0.4,8,0.5]])
    colies = ["x1","y1","h1","x2","y2","h2"]
    test_df = pd.DataFrame(foob_arr, columns = colies)
    return test_df

test_df = make_test_df()

print("This is the test dataframe")
print(test_df)

init_guess_df = init_guess(test_df)

print("This is the dataframe of Monte Carlo values")
print(init_guess_df)
