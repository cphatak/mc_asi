#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct.25.2020.

Script to initiate the rhomb parallel run for a given set of lattice
parameters.

@author: cphatak
"""
import numpy as np
import time as time
import matplotlib as mpl
#determine if linux
mpl.use('Agg') # for linux operation.
import multiprocessing as mp
from rhomb_mcsims_par import print_sysinfo,run_MC
from rhomb_analyzedata import rhomb_analyzedata
import os as os


#Running for mltiple values of a and s
a_vals = np.asarray([350])
latt_types = np.asarray(['dual_kagome'])
n_a = a_vals.size
n_l = latt_types.size

for ia in range(n_a):
    for jl in range(n_l):

        set_num = n_l*ia + jl + 1
        #
        # Running the parallel version.
        #
        # LAttice and MC parameters.
        num_runs = 10 #Number of runs in parallel to perform.
        a = a_vals[ia]#350.0 #Lattice parameter
        latt_type = latt_types[jl] #type of lattice.
        nx = 5 #num of islands along x
        ny = 5 #num of islands along y
        mc_iters = 1000 #number of MC iterations
        eq_iters = 10 #number of equilibriation iterations
        start_temp = 1.0e-5 #Start temperature
        end_temp = 1.0e-11 #end temperature
        n_temp = 300 #number of temperature steps
        red_fac = 'Geometric' #reduction factor - 'Linear' or 'Geometric'
        save_file = 2000 #Save file every N steps in the MC iters.
        load_file = False #load magnetic config data during MC runs
        file_name = "MCrun_mag_run0_199.txt"
        verbose = False
        display = True
        lattice_draw_step = 5
        
        #set working directory
        work_dir = '/Users/cphatak/ANL_work/spinice/rhomb_latt/dual_kag_MC/'+str(nx)+'x'+str(ny)+'_set'+str(set_num)
        
        #get system information
        print_sysinfo()
        
        #start multiprocessing computation
        nproc = np.minimum(num_runs, mp.cpu_count())
        print('Using {0:2d} processes.'.format(int(nproc)))
        st_time = time.time()
        pool = mp.Pool(processes = nproc)
        results = [pool.apply_async(run_MC, args=(x, #job ID
                                                  a, #Lattice parameter
                                                  nx, #num of islands along x
                                                  ny, #num of islands along y
                                                  latt_type, #Type of Rhombille lattice.
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
                                                  work_dir, #set the working directory for output
                                                  display)) for x in range(nproc)]
        
        results = [p.get() for p in results]
        print('Total time:',time.time()-st_time)
        print('Total runs completed:',nproc)
        
        #Once the run is complete, then we run analyze_data to average the results.
        res = rhomb_analyzedata(a = int(a), tot_runs = nproc, run_num = 0, lattice_draw_step = lattice_draw_step,
                                fldr = work_dir+'/')
        
        #then tar all the run output files.
        os.system("tar -zcf "+work_dir+"/results_runs_data.tar.gz "+work_dir+"/run*")
        
        #then delete all the run files
        os.system("rm -rf "+work_dir+"/run*")


