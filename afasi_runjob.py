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
from afasi_analyzedata import afasi_analyzedata
import os as os


#Running for mltiple values of a and s
a_vals = np.array([350, 450])
s_vals = np.array([120, 220])
n_a = a_vals.size
n_s = s_vals.size

for ia in range(n_a):
    for js in range(n_s):

        set_num = n_s*ia + js + 1
        #
        # Running the parallel version.
        #
        # LAttice and MC parameters.
        num_runs = 10 #Number of runs in parallel to perform.
        a = a_vals[ia]#350.0 #Lattice parameter
        s = s_vals[js]#120.0 #island separation
        nx = 3 #num of islands along x
        ny = 3 #num of islands along y
        mc_iters = 50 #number of MC iterations
        eq_iters = 0 #number of equilibriation iterations
        start_temp = 10.0 #Start temperature
        end_temp = 1.0 #end temperature
        n_temp = 20 #number of temperature steps
        red_fac = 'Linear' #reduction factor - 'Linear' or 'Geometric'
        save_file = 2000 #Save file every N steps in the MC iters.
        load_file = False #load magnetic config data during MC runs
        file_name = "MCrun_mag_run0_199.txt"
        verbose = True
        display = True
        lattice_draw_step = 5
        
        #set working directory
        work_dir = '/Users/cphatak/work/af_asi_sims/test/'+str(nx)+'x'+str(ny)+'_set'+str(set_num)
        
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
                                                  work_dir, #set the working directory for output
                                                  display)) for x in range(nproc)]
        
        results = [p.get() for p in results]
        print('Total time:',time.time()-st_time)
        print('Total runs completed:',nproc)
        
        #Once the run is complete, then we run analyze_data to average the results.
        res = afasi_analyzedata(a = int(a), s = int(s), tot_runs = nproc, run_num = 0, lattice_draw_step = lattice_draw_step,
                                fldr = work_dir+'/')
        
        #then tar all the run output files.
        os.system("tar -zcf "+work_dir+"/results_runs_data.tar.gz "+work_dir+"/run*")
        
        #then delete all the run files
        os.system("rm -rf "+work_dir+"/run*")


