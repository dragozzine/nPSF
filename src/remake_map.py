#
# remake_map.py
#
# Remakes the map_grid plots using the chain.h5, grid.npy, and llhoods.npy files
#
# William Giforos
# 06/24/22
#

import numpy as np
import emcee
import sys
import os
import commentjson as json

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

    
# Ensure that its being run from the right place, i.e. a results directory
cwd = os.getcwd()
if "results" in cwd:
    runprops = ReadJson("runprops.txt").outProps()
    folder = cwd
    image = runprops["input_image"]
else:
    print("remake_map needs to be run from a results directory.")
    sys.exit()

# map_plots setup
fullpath = folder + "/" + image
print(fullpath)
runprops["input_image"] = fullpath
runprops["image_name"] = image[:-5]
runprops["npsfs"] = 1

# Fake the backend parameter if no chain to load
if os.path.isfile('chain.h5') == True:
    backend = emcee.backends.HDFBackend('chain.h5')
else:
    backend = np.zeros((1,1))

# Import the names of the objects
name_dict = runprops.get("names_dict")
objectnames = []
for i in name_dict.values():
    objectnames.append(i)
print(objectnames)

# Set the resultspath
resultspath = "../" + runprops["image_name"] + "_" + runprops["map_identifier"]

# Load the saved grid and likelihoods
grid = np.load('grid.npy')
grid_llhoods = np.load('llhoods.npy')

# Rerun plots
from analysis import *
map_plots(backend, grid, grid_llhoods, resultspath, runprops)
    
print("Done")