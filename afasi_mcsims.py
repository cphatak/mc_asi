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
def init_afasi_latt(a = 350, # lattice parameter
                    s = 120, # island separation
                    nx = 3, # repeat along x
                    ny = 3, # repeat along y
                    centers = 0, # array to hold the centers
                    angles = 0, # array to hold angles
                    ):

    #compute total number of islands
    #n_isl = nx * ny * 6
    #centers = np.zeros([2,n_isl])
    #angles = np.zeros([n_isl])
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

    return 1
#------------------------------------------------------------------
#
# Calc_Energy function for AFASI
#
def calc_energy(distmap = 0, # distance map of islands
                centers = 0, # array of centers of islands
                mag = 0, # array of magnetization of islands
                max_nn_dist = 500, # max. distance of nearest neighbors
                max_nn_num = 9, # max. number of nearest neighbors
                ):

    # number of islands
    nm, n_isl = centers.shape

    # Energy variable
    tot_energy = 0

    # loop over each island and compute the neighbors
    for i in range(n_isl):
        inds = np.argsort(distmap[i,:])
        for cnt in range(1,max_nn_num):
            j = inds[cnt]
            if (distmap[i,j] <= max_nn_dist):
                si_sj = mag[0,i]*mag[0,j] + mag[1,i]*mag[1,j]
                r_ij = distmap[i,j]
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
s = 120.0 #island separation
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

centers = np.zeros([2,n_isl])
angles = np.zeros([n_isl])
#next we initialize the lattice.
lattice = init_afasi_latt(a = a, s = s, nx = nx, ny = ny,
                       centers = centers, angles = angles)

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
curr_energy = calc_energy(distmap = distmap, centers = centers, mag = mag,
                          max_nn_dist = max_nn_dist, max_nn_num = max_nn_num)


#compute total number of temperature steps
