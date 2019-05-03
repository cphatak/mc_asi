#
#
# This set of programs and routines are to perform Monte Carlo
# simulations on the AF-ASI lattice using dipolar coupling as the
# the dominant energy interactions.
#
# Written, CD Phatak, ANL, 21.Feb.2019.
#
# Modified to use the Dipolar_MC Class file.

#import necessary modules
import numpy as np
import os as os
import sys as sys
import datetime
import time as time
from scipy.spatial import cKDTree as sp_cKDTree
from matplotlib import pyplot as plt
from dipolar_MC import Dipolar_MC

#------------------------------------------------------------------
#
# Function using the kdtree algorithm to find the nearest neighbors.
# USing the method described here - https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates
#
def do_kdtree(combined_x_y_arrays,points,max_nn,max_dd):
    mytree = sp_cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points, k=max_nn, distance_upper_bound=max_dd)
    return indexes

#------------------------------------------------------------------
#
# Init lattice function for AFASI
#
# This function will initialize the lattice for given set of lattice parameters
# and num of islands along x and y. The output will be an array centers consisting
# of positions of each island (x,y), and an array angles with angle of each island
# for magnetization, and an array nn_inds consisting of nearest neighbor indices
# for each island to be considered for dipolar interactions.
#
def init_afasi_latt(a = 350, # lattice parameter
                    s = 120, # island separation
                    nx = 3, # repeat along x
                    ny = 3, # repeat along y
                    max_nn_dist = 500, # max. distance of nearest neighbors
                    max_nn_num = 9 # max. number of nearest neighbors
                    ):

    #compute total number of islands
    n_isl = nx * ny * 6
    centers = np.zeros([2,n_isl])
    angles = np.zeros([n_isl])
    count = 0

    for i in range(nx):
        for j in range(ny):
            #horizontal islands
            angles[count] = 0
            angles[count+1] = 0
            centers[0,count] = i*2*a - a/2 + j*a
            centers[1,count] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 + s/2
            centers[0,count+1] = i*2*a - a/2 + j*a
            centers[1,count+1] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 - s/2
            #first set of rotated islands
            angles[count+2] = 120
            angles[count+3] = 120
            centers[0,count+2] = i*2*a + a/2 + j*a + np.sqrt(3)/4*s
            centers[1,count+2] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 + s/4
            centers[0,count+3] = i*2*a + a/2 + j*a - np.sqrt(3)/4*s
            centers[1,count+3] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 - s/4
            #second set of rotated islands
            angles[count+4] = 60
            angles[count+5] = 60
            centers[0,count+4] = i*2*a + j*a - np.sqrt(3)*s/4
            centers[1,count+4] = j*np.sqrt(3)*a + a*np.sqrt(3)/4 + s/4
            centers[0,count+5] = i*2*a + j*a + np.sqrt(3)*s/4
            centers[1,count+5] = j*np.sqrt(3)*a + a*np.sqrt(3)/4 - s/4

            #increment count
            count += 6

    #now to use the cKDTree method
    comb_xy = centers.transpose()
    p_list = list(comb_xy)
    nn_inds = do_kdtree(comb_xy,p_list,max_nn_num+1,max_nn_dist)

    return [centers,angles,nn_inds]

#
#------------------------------------------------------------------
#
# the main function for defining parameters for MC simulations.

jobID = 'run1',
a = 350.0 #lattice parameter
s = 150.0 #island separation
nx = 5 #repeat along x
ny = 5 #repeat along y,
mc_iters = 1000 #number of MC iterations
eq_iters = 0 #number of equilibriation iterations
start_temp = 100 #Start temperature
end_temp = 1 #end temperature
red_fac = 0.95 #reduction factor
save_file = 500 #save config data during MC runs
verbose = True
display = False

#Set the next nearest neghbors
max_nn_num = 9
max_nn_dist = 1.5 * a

#Compute number of islands.
#6 islands per motif.
n_isl = nx * ny * 6
print('Number of islands:',n_isl)

centers = np.zeros([2,n_isl])
angles = np.zeros([n_isl])
nn_inds = np.zeros([n_isl,max_nn_num])

#next we initialize the lattice.
[centers,angles,nn_inds] = init_afasi_latt(a = a, s = s, nx = nx, ny = ny, max_nn_dist = max_nn_dist,
                          max_nn_num = max_nn_num)#, centers = centers, angles = angles, nn_inds = nn_inds)

#next we initialize the Dipolar_MC Class for MC sims
dipolar_MC1 = Dipolar_MC(centers = centers, angles = angles, nn_inds = nn_inds,
                         max_nn_num = max_nn_num, max_nn_dist = max_nn_dist)

res = dipolar_MC1.Calc_Energy()
print("Current Energy:",dipolar_MC1.energy)

#compute total number of temperature steps
#checking how long one MC run at high temperature takes.
dipolar_MC1.mc_iters = 100
dipolar_MC1.eq_iters = eq_iters
dipolar_MC1.temp = 0.1

start = time.time()
res = dipolar_MC1.MC_move(verbose=True)
end = time.time()
print("Runtime :",end-start)
print("End energy:",dipolar_MC1.energy)
