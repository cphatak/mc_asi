#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 12:45:23 2019

Script for aggregating the data resulting from the MC sims run in parallel
We will average the data from all runs and save average, and std.dev values
per temperature step.

@author: cphatak
"""
import numpy as np
from matplotlib import pyplot as plt
import datetime
import os
#from draw_lattice import draw_lattice

def distASI_analyzedata(a = 350, #lattice parameter
                      tot_runs = 10, # total number of runs
                      run_num = 1, #Run number to plot magnetization
                      lattice_draw_step = 10, #draw lattice vector
                      fldr = '12x12_set1/results_runs_data/',
                      draw_col_lattice = False #drawing the islands with color map.
                      ):

    #---------------------------------------------------------------------
    #directory base
    subdir_base = 'run'
    fname_base = 'Dipolar_MC1_'
    magname_base = 'MCrun_mag_'
    cenname_base = 'MCrun_lattice_coords_'
    num_runs = tot_runs
    
    #results folder
    #cwd = os.getcwd()
    opdir = fldr + 'avg_results/'
    
    #create directory if it does not exist.
    if not os.path.exists(opdir):
        os.makedirs(opdir)
    
    
    #number of header lines to skip
    header_skip = 19
    
    #process all the data
    fnames = []
    for num_run in range(num_runs):
        fnames.append(fldr + subdir_base + str(num_run) + '/' + 
                      fname_base + subdir_base + str(num_run) + '.txt')
        
    
    data_list = [np.genfromtxt(filename,delimiter=',',skip_header=header_skip) for filename in fnames]

    #we need to trim the data to the min index
    min_ind = np.min([len(listitem) for listitem in data_list])
    data = np.asarray([i[0:min_ind,:] for i in data_list])
    
    #compute the averages
    avg_data = np.mean(data,axis=0)
    stddev_data = np.std(data,axis=0)
    
    nvals, nvars = avg_data.shape
    
    #get the header from one of the files.
    fname = fldr + subdir_base + str(run_num) + '/' + fname_base + subdir_base + str(run_num) + '.txt'
    hdr=[]
    with open(fname,'r') as reader:
        for i in range(header_skip):
            hdr.append(reader.readline())
        
    
    
    #save the averaged and stddev data
    data_file = opdir + 'averaged_MC_data.txt'
    f = open(data_file,"w+")
    d = datetime.datetime.now()
    f.write('# Averaged Data from MC runs.\n')
    f.write('# Created: C. Phatak, ANL \n')
    f.write('#{:%Y-%m-%d %H:%M:%S} \n'.format(d))
    f.write('# Header from MC runs data follows: \n')
    f.write(''.join(hdr))
    for i in range(nvals):
        f.write('{0:.4e}, {1:.4e}, {2:.3f}, {3:.4e}, {4:.4e}, {5:5f}, {6:5f}, {7:5f}\n'.format(avg_data[i,0],avg_data[i,1],avg_data[i,2],avg_data[i,3],avg_data[i,4],avg_data[i,5],avg_data[i,6],avg_data[i,7]))
    
    f.close()
    
    data_file = opdir + 'stddev_MC_data.txt'
    f = open(data_file,"w+")
    d = datetime.datetime.now()
    f.write('# Std Dev Data from MC runs.\n')
    f.write('# Created: C. Phatak, ANL \n')
    f.write('#{:%Y-%m-%d %H:%M:%S} \n'.format(d))
    f.write('# Header from MC runs data follows: \n')
    f.write(''.join(hdr))
    for i in range(nvals):
        f.write('{0:.4e}, {1:.4e}, {2:.3f}, {3:.4e}, {4:.4e}, {5:5f}, {6:5f}, {7:5f}\n'.format(stddev_data[i,0],stddev_data[i,1],stddev_data[i,2],stddev_data[i,3],stddev_data[i,4],stddev_data[i,5],stddev_data[i,6],stddev_data[i,7]))
    
    f.close()
    
    #saving the magnetization vector maps.
    #load the centers data
    cen_file = fldr + subdir_base + str(run_num) + '/' + cenname_base + subdir_base + str(run_num)
    centers = np.genfromtxt(cen_file + '.txt',delimiter=',',skip_header=2)
    
    for i in range(0,nvals,lattice_draw_step):
        mag_file = fldr + subdir_base + str(run_num) + '/' + magname_base + subdir_base + str(run_num) + '_' + str(i)
        mag_data = np.genfromtxt(mag_file + '.txt', delimiter=',', skip_header=1)
        
        #save the plot
        fig, ax1 = plt.subplots(figsize=(8,8))
        plt.rcParams['image.cmap'] = 'Paired'
        ax1.set_title('a = {0:3d}, T = {1:.4e}'.format(a,avg_data[i,0]))
        col_arr = np.arctan2(mag_data[:,1],mag_data[:,0])
        q1 = ax1.quiver(centers[:,0],centers[:,1],mag_data[:,0],mag_data[:,1],col_arr,pivot='mid',scale=25,headwidth=5)
        #plt.draw()
        plt.savefig(opdir+'Lattice_state_'+str(i)+'.jpg',bbox_inches='tight', dpi=150, quality=95)
        plt.close()

        #calling draw lattice.
        #blank = 0
        #if (draw_col_lattice):
        #    res = draw_lattice(blank,cen_fname = cen_file, mag_fname = mag_file, dim = 1000, del_px = 15, Lx = 300, Ly = 100, thk = 10, save_tfs = False, save_lattice = True,
        #                       buff_offset = 200)
    
    print('Averaging and Plotting complete.\n')
    
    #end
    return 0

