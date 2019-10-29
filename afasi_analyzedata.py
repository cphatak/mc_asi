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
from draw_lattice import draw_lattice

def afasi_analyzedata(a = 350, #lattice parameter
                      s = 120, #separation
                      run_num = 0, #Run number to plot magnetization
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
    num_runs = 10
    
    #results folder
    #cwd = os.getcwd()
    opdir = fldr + 'avg_results/'
    
    #create directory if it does not exist.
    if not os.path.exists(opdir):
        os.makedirs(opdir)
    
    
    #number of header lines to skip
    header_skip = 16
    
    #process all the data
    fnames = []
    for num_run in range(num_runs):
        fnames.append(fldr + subdir_base + str(num_run) + '/' + 
                      fname_base + subdir_base + str(num_run) + '.txt')
        
    
    data = np.asarray([np.genfromtxt(filename,delimiter=',',skip_header=header_skip) for filename in fnames])
    
    
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
        ax1.set_title('a = {0:3d}, s = {1:3d}, T = {2:.4e}'.format(a,s,avg_data[i,0]))
        col_arr = np.arctan2(mag_data[:,1],mag_data[:,0])
        q1 = ax1.quiver(centers[:,0],centers[:,1],mag_data[:,0],mag_data[:,1],col_arr,pivot='mid',scale=25,headwidth=5)
        #plt.draw()
        plt.savefig(opdir+'Lattice_state_'+str(i)+'.pdf',bbox_inches='tight')
        plt.close()

        #calling draw lattice.
        blank = 0
        if (draw_col_lattice):
            res = draw_lattice(blank,cen_fname = cen_file, mag_fname = mag_file, dim = 1000, del_px = 15, Lx = 300, Ly = 100, thk = 10, save_tfs = False, save_lattice = True,
                               buff_offset = 200)
    
    print('Averaging and Plotting complete.\n')
    
    #end
    return 0


"""
Function to compute the AF/FM pairing for each island pair and then save the data
as a function of temperature for a given run and save plots to show the spatial
location of the AF/FM pairs.

@author, CD Phatak, ANL, 28.Oct.2019.
"""

def afasi_isingpairs(a = 350, #lattice parameter
                     s = 120, #separation distance
                     run_num = 0, #run number to compute the data for
                     fldr = '12x12_set1/results_runs_data/'
                     ):
    
    #---------------------------------------------------------------------
    #directory base
    subdir_base = 'run'
    fname_base = 'Dipolar_MC1_'
    magname_base = 'MCrun_mag_'
    cenname_base = 'MCrun_lattice_coords_'
    skip_header = 16
    
    #results folder
    #cwd = os.getcwd()
    opdir = fldr + 'isingpair_results/'
    
    #create directory if it does not exist.
    if not os.path.exists(opdir):
        os.makedirs(opdir)
    
    #read the Dipolar MC output data file and get the number of temperature
    #steps and temperature values.
    data_file = fldr + subdir_base + str(run_num) + '/' + fname_base + subdir_base + str(run_num) + '.txt'
    data = np.genfromtxt(data_file, delimiter=',', skip_header = skip_header)
    
    #get the header for saving later
    hdr=[]
    with open(data_file,'r') as reader:
        for i in range(skip_header):
            hdr.append(reader.readline())
    
    #get the number of temperature steps
    ntemp, nvars = data.shape
    
    #load the centers data
    cen_file = fldr + subdir_base + str(run_num) + '/' + cenname_base + subdir_base + str(run_num)
    centers = np.genfromtxt(cen_file + '.txt',delimiter=',',skip_header=2)
    
    #getnumber of islands
    n_isl, vv = centers.shape
    
    #loop over each island pair for a given temperature and determine whether
    #islands are AF or FM.
    isingpair_file = opdir + 'Isingpair_data.txt'
    f = open(isingpair_file,"w+")
    d = datetime.datetime.now()
    f.write('# Ising Pair Data from MC run.\n')
    f.write('# Created: C. Phatak, ANL \n')
    f.write('#{:%Y-%m-%d %H:%M:%S} \n'.format(d))
    f.write('# Header from MC runs data follows: \n')
    f.write(''.join(hdr))
    
    for it in range(ntemp):
        mag_file = fldr + subdir_base + str(run_num) + '/' + magname_base + subdir_base + str(run_num) + '_' + str(it)
        mag_data = np.genfromtxt(mag_file + '.txt', delimiter=',', skip_header=1)
        
        #counter for AF and FM sites.
        af_count = 0
        fm_count = 0
        af_col = 'r'
        fm_col = 'b'
        
        fig, ax1 = plt.subplots(figsize=(8,8))
        ax1.set_title('a = {0:3d}, s = {1:3d}, T = {2:.4e}'.format(a,s,data[it,0]))
        for isl in range(0,n_isl,2):
            if (mag_data[isl,1] == 0):
                #horizontal island
                if (mag_data[isl,0]*mag_data[isl+1,0] < 0):
                    af_count += 1
                    c = af_col
                else:
                    fm_count += 1
                    c = fm_col
            else:
            
                if (np.sign(mag_data[isl,1]*mag_data[isl,0]) == np.sign(mag_data[isl+1,0]*mag_data[isl+1,1])):
                    #60 deg CCW set
                    if (np.sign(mag_data[isl,1]) == np.sign(mag_data[isl+1,1])):
                        fm_count += 1
                        c = fm_col
                    else:
                        af_count += 1
                        c = af_col
                else:
                    #60 deg CW set
                    if (np.sign(mag_data[isl,1]) == np.sign(mag_data[isl+1,1])):
                        fm_count += 1
                        c = fm_col
                    else:
                        af_count += 1
                        c = af_col
            
            #plot the data
            cent_x = (centers[isl,0]+centers[isl+1,0])/2
            cent_y = (centers[isl,1]+centers[isl+1,1])/2
            ax1.scatter(cent_x,cent_y,c=c)
        
        #out of island loop
        plt.savefig(opdir + 'IsingPairMap_'+str(it)+'.pdf',bbox_inches='tight')
        plt.close()
        #save the data into file
        f.write('{0:.4e}, {1:3d}, {2:3d}\n'.format(data[it,0],af_count,fm_count))
    
    #out of temp loop
    f.close()
    
    #end
    return 1


    
        
    