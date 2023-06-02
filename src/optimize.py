#
#	optimize.py
#
#	Adapted from Multimoon's mm_optimize.py file, which uses an optimizer to get a better initial guess and was written by Benjamin Proudfoot (09/28/20).
#
#	uses an optimizer to guess the bestfit parameters
#
#	William Giforos
#	05/25/23
#

import numpy as np
import pandas as pd
import random
import commentjson as json
import emcee
import scipy
from tqdm import tqdm
import functools
from datetime import datetime
import likelihood

def neg_log_prob(params, image, psfs, focuses, runprops):
	out = -1*likelihood.log_likelihood(params, image, psfs, focuses, runprops)
	return out

def op_map(i, p0, image, psfs, focuses, runprops, maxiter):
    # Minimization of scalar function of one or more variables.
    return scipy.optimize.minimize(neg_log_prob, p0[i,:], args = (image, psfs, focuses, runprops), method='Nelder-Mead', options={'maxiter':maxiter, 'disp':False})

# Copy of the parallelization code in Multimoon. As Haumea can't do this task as far as I know, I have defined func_optimize below
def pool_optimize(nwalkers, p0, image, psfs, focuses, runprops, pool):
    # basically calling op_map as a partial function named optimization
	optimization = functools.partial(op_map, p0=p0, float_names=float_names, fixed_df = fixed_df, total_df_names=total_df_names, fit_scale=fit_scale, runprops=runprops, obsdf=obsdf, geo_obj_pos=geo_obj_pos, best_llhoods=best_llhoods)
	x = tqdm(range(nwalkers))
	begin = datetime.now()  
	#pool.map is a way of parallelizing this process for each walker involved.
	data = pool.map(optimization, x)
	print(datetime.now()-begin)    
	for i in range(len(data)):
		p0[i,:]= data[i].x
    
	return p0 
        
def func_optimize(nwalkers, p0, image, psfs, focuses, runprops, maxiter):    
	op_data = p0[:nwalkers,:]
	llhoods = np.zeros(nwalkers)
	for i in tqdm(range(nwalkers)):
		#print("p0:",p0[i,:])
		result = op_map(i, p0, image, psfs, focuses, runprops, maxiter)
		#print("optimized:",result.x)
		op_data[i,:] = result.x
		llhoods[i] = -1*result.fun
        
	return op_data, llhoods
