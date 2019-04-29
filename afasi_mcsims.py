#
#
# This set of programs and routines are to perform Monte Carlo
# simulations on the AF-ASI lattice using dipolar coupling as the
# the dominant energy interactions.
#
# Written, CD Phatak, ANL, 21.Feb.2019.
#

#import necessary modules
import numpy as np
import os as os
import sys as sys
import datetime
import time as time
from scipy.spatial import cKDTree as sp_cKDTree
from matplotlib import pyplot as plt

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
#------------------------------------------------------------------
#
# Calc_Energy function for AFASI
#
def calc_energy(centers = 0, # array of centers of islands
                mag = 0, # array of magnetization of islands
                nn_inds = 0 # array for nearest neighbors for dipolar interaction.
                ):

    # number of islands
    n_isl,num_neighbors = nn_inds.shape

    # Energy variable
    tot_energy = 0

    # loop over each island and compute the neighbors
    for i in range(n_isl):
        for cnt in range(num_neighbors-1):
            j = nn_inds[i,cnt+1].astype('int')
            if ((i != j) and (j != n_isl)):
                si_sj = mag[0,i]*mag[0,j] + mag[1,i]*mag[1,j]
                r_ij = np.sqrt((centers[0,i]-centers[0,j])**2 +
                               (centers[1,i]-centers[1,j])**2)
                if (r_ij == 0):
                        print(i,j)
                si_rij = (centers[0,i]-centers[0,j])*mag[0,i] + (centers[1,i]-centers[1,j])*mag[1,i]
                sj_rji = (centers[0,j]-centers[0,i])*mag[0,j] + (centers[1,j]-centers[1,i])*mag[1,j]
                tot_energy += (si_sj/r_ij**3 - (3.0*si_rij*sj_rji)/r_ij**5)

    #return total energy
    return tot_energy/2.0

#------------------------------------------------------------------
#
# MC_move function to actual run the MC simulation.
#
# This function will take the input parameters for the MC simulation
# then for number of MC iters, randomly select a spin site and flip its
# direction, calculate energy and check whether to accept it or not.
# It will also calculate the thermodynamic parameters such as sp. heat,
# and susceptibility during the MC iters for a given temperature.
#
def MC_move(centers = 0, #centers of islands
            mag = 0, #magnetization of islands
            nn_inds = 0, #array for holding nearest neighbors
            mc_iters = 10000, #number of MC iters per temp. step
            eq_iters = 0, #number of equilibriation iterations.
            temp = 700 #temperature variable for MC run
            ):

    #initialize arrays for holding various quantities
    avg_en = 0.0
    avg_en2 = 0.0
    avg_mag = 0.0
    avg_mag2 = 0.0
    n_lowaccept = 0
    n_highaccept = 0
    n_noaccept = 0
    
    #get current energy
    curr_energy = Calc_Energy(centers = centers, mag = mag, nn_inds = nn_inds)
    
    #beta_value
    mu0 = 4.0*!pi*1e-7
    mult_fac =

    #num of islands
    temp, n_isl = centers.shape

    for nn in range(mc_iters):
        for ii in range(n_isl):
            # pick a random site in the lattice.
            site = np.random.randint(0,n_isl)

            #change the magnetization
            mag[*,site] *= (-1)

            #calculate the energy
            new_energy = Calc_Energy(centers = centers, mag = mag, nn_inds = nn_inds)

            #difference in energy
            dE = new_energy - curr_energy

            #check if we should accept this energy or not
            if (dE lt 0):
                n_lowaccept += 1
                curr_energy = new_energy


#
#------------------------------------------------------------------
#
# the main function for defining parameters for MC simulations.

#def afasi_mcrun(jobID = 'run1',
a = 350.0 #lattice parameter
s = 150.0 #island separation
nx = 5 #repeat along x
ny = 5 #repeat along y,
mc_iters = 1000 #number of MC iterations
eq_iters = 0 #number of equilibriation iterations
start_temp = 1000 #Start temperature
end_temp = 1 #end temperature
red_fac = 0.95 #reduction factor
save_file = 500 #save config data during MC runs
verbose = True
display = False
    #                ):

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


##next we compute the distance matrix.
#distmap = np.zeros([n_isl,n_isl])
#for i in range(n_isl):
#    for j in range(n_isl):
#        distmap[i,j] = np.sqrt((centers[0,i]-centers[0,j])**2 +
#                               (centers[1,i]-centers[1,j])**2)
#
#next we define the magnetization array
mag = np.zeros([2,n_isl])
mag[0,:] = np.cos(np.deg2rad(angles))
mag[1,:] = np.sin(np.deg2rad(angles))

#next we calculate energy
curr_energy = calc_energy(centers = centers, mag = mag, nn_inds = nn_inds)


#compute total number of temperature steps
