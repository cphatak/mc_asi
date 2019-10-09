#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 12:40:13 2019

Script to initiate the afasi parallel run for a given set of lattice
parameters.

@author: cphatak
"""
import numpy as np
import time as time
import matplotlib as mpl
#determine if linux
mpl.use('Agg') # for linux operation.
import multiprocessing as mp
from afasi_mcsims_par import print_sysinfo,run_MC

#
# Running the parallel version.
#
# LAttice and MC parameters.
num_runs = 10 #Number of runs in parallel to perform.
a = 350.0 #Lattice parameter
s = 120.0 #island separation
nx = 3 #num of islands along x
ny = 3 #num of islands along y
mc_iters = 1000 #number of MC iterations
eq_iters = 0 #number of equilibriation iterations
start_temp = 2000.0 #Start temperature
end_temp = 1.0 #end temperature
n_temp = 200 #number of temperature steps
red_fac = 'Linear' #reduction factor - 'Linear' or 'Geometric'
save_file = 2000 #Save file every N steps in the MC iters.
load_file = True #load magnetic config data during MC runs
file_name = "MCrun_mag_run0_199.txt"
verbose = True
display = True

#get system information
print_sysinfo()

#start multiprocessing computation
nproc = np.minimum(num_runs, mp.cpu_count())
print('Using {0:2d} processes.'.format(int(nproc)))
st_time = time.time()
pool = mp.Pool(processes = nproc)
results = [pool.apply_async(run_MC, args=(x, #job ID
                                          a, #Lattice parameter
                                          s, #island separation
                                          nx, #num of islands along x
                                          ny, #num of islands along y
                                          mc_iters, #number of MC iterations
                                          eq_iters, #number of equilibriation iterations
                                          start_temp, #Start temperature
                                          end_temp, #end temperature
                                          n_temp, #number of temperature steps
                                          red_fac, #reduction factor
                                          load_file, #load magnetic config data for MC runs
                                          save_file, #Save intermediate files every N steps in MC iters.
                                          file_name, #file name for reading the data from
                                          verbose,
                                          display)) for x in range(nproc)]

results = [p.get() for p in results]
print('Total time:',time.time()-st_time)
print('Total runs completed:',nproc)


