#
# remake_plots.py
#
# Remakes plots using the chain.h5 file
#
# Benjamin Proudfoot
# 01/24/22
#

import glob
import os
import sys
import commentjson as json
import emcee
import numpy as np

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

# Ensure that its being run from the right place
if 'results' in os.getcwd():
    runprops = ReadJson('runprops.txt').outProps()
else:
    print("remake_plots needs to be run from a results directory.")
    sys.exit()
    
# Load in chain file and image/psfs arrays
resultspath = os.getcwd()
backend = emcee.backends.HDFBackend(resultspath + '/chain.h5')

#os.chdir("../../src")
os.chdir("../../../src")

# Rerun plots
from analysis import *
plots(backend, resultspath, runprops)
