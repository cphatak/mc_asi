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

    # next we calculate the distance map for the lattice
    distmap = np.zeros([n_isl,n_isl])
    for i in range(n_isl):
        for j in range(n_isl):
            distmap[i,j] = np.sqrt((centers[0,i]-centers[0,j])**2 +
                                   (centers[1,i]-centers[1,j])**2)

    # next we go through each island and sort the indices of the nearest neighbors.
    nn_inds = np.zeros([n_isl,max_nn_num])
    for i in range(n_isl):
        inds = np.argsort(distmap[i,:])
        for cnt in range(1,max_nn_num+1):
            j = inds[cnt]
            if (distmap[i,j] <= max_nn_dist):
                nn_inds[i,cnt-1] = j

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
        for cnt in range(num_neighbors):
            j = nn_inds[i,cnt].astype('int')
            if (i != j):
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
# the main function for running the MC simulations.

#def afasi_mcrun(jobID = 'run1',
a = 350.0 #lattice parameter
s = 150.0 #island separation
nx = 11 #repeat along x
ny = 11 #repeat along y,
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

centers = np.zeros([2,n_isl])
angles = np.zeros([n_isl])
nn_inds = np.zeros([n_isl,max_nn_num])
#next we initialize the lattice.
[centers,angles,nn_inds] = init_afasi_latt(a = a, s = s, nx = nx, ny = ny, max_nn_dist = max_nn_dist,
                          max_nn_num = max_nn_num)#, centers = centers, angles = angles, nn_inds = nn_inds)


#next we compute the distance matrix.
distmap = np.zeros([n_isl,n_isl])
for i in range(n_isl):
    for j in range(n_isl):
        distmap[i,j] = np.sqrt((centers[0,i]-centers[0,j])**2 +
                               (centers[1,i]-centers[1,j])**2)

#next we define the magnetization array
mag = np.zeros([2,n_isl])
mag[0,:] = np.cos(np.deg2rad(angles))
mag[1,:] = np.sin(np.deg2rad(angles))

#next we calculate energy
curr_energy = calc_energy(centers = centers, mag = mag, nn_inds = nn_inds)


#compute total number of temperature steps
