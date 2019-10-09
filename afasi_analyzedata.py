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

#directory base
dir_base = '/home/cphatak/af_asi_sims/'
latts = '12x12_runs1/'
subdir_base = 'run'
fname_base = 'Dipolar_MC1_'
num_runs = 10

#number of header lines to skip
header_skip = 16

#process all the data
fnames = []
for num_run in range(num_runs):
    fnames.append(subdir_base + str(num_run) + '/' + 
                  fname_base + subdir_base + str(num_run) + '.txt')
    

data = np.asarray([np.genfromtxt(filename,delimiter=',',skip_header=header_skip) for filename in fnames])


#compute the averages
avg_data = []
stddev_data = []

for i in range(0,num_runs*latt_sets,num_runs):
    avg_data.append(np.mean(data[i:i+num_runs-1,:,:],axis=0))
    stddev_data.append(np.std(data[i:i+num_runs-1,:,:],axis=0))

avg_data = np.mean(data,axis=0)
stddev_data = np.std(data,axis=0)

#get the header from one of the files.
hdr = np.genfromtxt(subdir_base + str(0) + '/' +
                    fname_base + subdir_base + str(0) + '.txt', max_rows = header_skip)

#save the averaged and stddev data
data_file = 'averaged_MC_data.txt'
f = open(data_file,"w+")
d = datetime.datetime.now()
f.write(hdr)


    