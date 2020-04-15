# testing_and_examples.py
# Written by Darin Ragozzine
# tests nPSF codes and thereby shows examples of usage
# 
# April 15, 2020


# assume that we are working in the "src" subdirectory of the main nPSF directory
# we need to first import all the files and their functions
# start by adding the current directory to the python path so it knows where to look
import sys
sys.path.insert(1,".")


# import various packages we'll need
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import numpy
import emcee
import corner
import os
import json
import pandas

# there are ways to load all the commands from all the folders automatically
# but I think for now it makes sense to just put these in by hand
from insertpsf import *
from getpsf import *
from init_guess import *
from likelihood import *
from npsf_run import *
from params import *
from priors import * 
from testsupport import *

print("success loading all packages and files")


# Okay, let's try getting the initial conditions set up
startguessdf=pandas.read_csv("../data/idsj03wgq_startguess.csv")

# init guess is not working
#initguessdf=init_guess(startguessdf)


