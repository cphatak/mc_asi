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

#directory base
subdir_base = 'run'
fname_base = 'Dipolar_MC1_'
magname_base = 'MCrun_mag_'
cenname_base = 'MCrun_lattice_coords_'
lattice_draw_step = 10 #draw lattice vector maps every N Temp. steps.
num_runs = 10

#results folder
op_folder = 'results/'
cwd = os.getcwd()
opdir = cwd + '/results/'

#create directory if it does not exist.
if not os.path.exists(opdir):
    os.makedirs(opdir)


#number of header lines to skip
header_skip = 16

#process all the data
fnames = []
for num_run in range(num_runs):
    fnames.append(subdir_base + str(num_run) + '/' + 
                  fname_base + subdir_base + str(num_run) + '.txt')
    

data = np.asarray([np.genfromtxt(filename,delimiter=',',skip_header=header_skip) for filename in fnames])


#compute the averages
avg_data = np.mean(data,axis=0)
stddev_data = np.std(data,axis=0)

nvals, nvars = avg_data.shape

#get the header from one of the files.
fname = subdir_base + str(0) + '/' + fname_base + subdir_base + str(0) + '.txt'
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
    f.write('{0:.4e}, {1:.4e}, {2:.3f}, {3:.4e}, {4:.4e}, {5:5f}, {6:5f}, {7:5f}\n'.format(avg_data[:,0],avg_data[:,1],avg_data[:,2],avg_data[:,3],avg_data[:,4],avg_data[:,5],avg_data[:,6],avg_data[:,7]))

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
    f.write('{0:.4e}, {1:.4e}, {2:.3f}, {3:.4e}, {4:.4e}, {5:5f}, {6:5f}, {7:5f}\n'.format(stddev_data[:,0],stddev_data[:,1],stddev_data[:,2],stddev_data[:,3],stddev_data[:,4],stddev_data[:,5],stddev_data[:,6],stddev_data[:,7]))

f.close()

#saving the magnetization vector maps.
#load the centers data
cen_file = subdir_base + str(0) + '/' + cenname_base + subdir_base + str(0) + '.txt'
centers = np.genfromtxt(cen_file,delimiter=',',skip_header=2)

