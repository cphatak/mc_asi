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
# the main function for running the MC simulations.

def afasi_mcrun(jobID = 'run1',
                a = 350, #lattice parameter
                s = 120, #island separation
                nx = 3, #repeat along x
                ny = 3, #repeat along y,
                mc_iters = 1000, #number of MC iterations
                eq_iters = 0, #number of equilibriation iterations
                start_temp = 1000, #Start temperature
                end_temp = 1, #end temperature
                red_fac = 0.95, #reduction factor
                save_file = 500, #save config data during MC runs
                verbose = True,
                display = False
                ):

    #Set the next nearest neghbors
    max_nn_num = 9
    max_nn_dist = 1.5 * a

    #Compute number of islands.
    #6 islands per motif.
    n_isl = nx * ny * 6

    #next we initialize the lattice.
    lattice = init_lattice(a = a, s = s, nx = nx, ny = ny,
                           centers = centers, angles = angles)


