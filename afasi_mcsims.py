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

jobID = 'run1'
a = 350.0 #lattice parameter
s = 120.0 #island separation
nx = 11 #repeat along x
ny = 11 #repeat along y,
mc_iters = 1000 #number of MC iterations
eq_iters = 0 #number of equilibriation iterations
start_temp = 1000 #Start temperature
end_temp = 1 #end temperature
n_temp = 100 #number of temperature steps
red_fac = 0.90 #reduction factor
save_file = 500 #save config data during MC runs
verbose = True
display = True
dir = '/Users/cphatak/work/af_test/'

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

dipolar_MC1.energy = dipolar_MC1.Calc_Energy()
print("Current Energy:",dipolar_MC1.energy)

#save the images if display option is set
if (display):
    
    fig1, ax1 = plt.subplots(figsize=(8,8))
    ax1.set_title('MC sims initial state')
    qq = ax1.quiver(dipolar_MC1.centers[0,:],dipolar_MC1.centers[1,:],dipolar_MC1.magx,dipolar_MC1.magy,pivot='mid')
    plt.draw()
    plt.savefig(dir+'initial_state.png',bbox_inches='tight')
    plt.close()

#saving the centers and initial mag data
f1 = open(dir+'MCrun_lattice_coords_'+jobID+'.txt','w+')
f1.write('# Centers of the islands for lattice run \n')
f1.write('# Num islands {0:3d} \n'.format(dipolar_MC1.n_isl))
for i in range(dipolar_MC1.n_isl):
    f1.write('{0:.3f}, {1:.3f} \n'.format(dipolar_MC1.centers[0,i],dipolar_MC1.centers[1,i]))

f1.close()

#save initial magnetization data
f2 = open(dir+'MCrun_mag_initial_'+jobID+'.txt', 'w+')
f2.write('# Num islands {0:3d} \n'.format(dipolar_MC1.n_isl))
for i in range(dipolar_MC1.n_isl):
    f2.write('{0:.3f}, {1:.3f} \n'.format(dipolar_MC1.magx[i],dipolar_MC1.magy[i]))

f2.close()


#MC iterations
dipolar_MC1.mc_iters = mc_iters
dipolar_MC1.eq_iters = eq_iters

#compute total number of temperature steps
temp_var = np.linspace(start_temp,end_temp,n_temp)
#variables to hold the values
latt_energy = np.zeros([n_temp])
latt_spheat = np.zeros([n_temp])
latt_susc = np.zeros([n_temp])
latt_mag = np.zeros([n_temp])

##open data file to save the runtime data
data_file = "Dipolar_MC1_"+jobID+".txt"
f = open(dir+data_file,"w+")
d = datetime.datetime.now()
f.write('# Runtime data file for variables.\n')
f.write('# Created: C. Phatak, ANL \n')
f.write('#{:%Y-%m-%d %H:%M:%S} \n'.format(d))
f.write('#\n')
#data dump of MC params
f.write('# a = {0:.1f}\n'.format(a))
f.write('# s = {0:.1f}\n'.format(s))
f.write('# nx = {0:2d}\n'.format(nx))
f.write('# ny = {0:2d}\n'.format(ny))
f.write('# Num isl = {0:4d}\n'.format(n_isl))
f.write('# MC iters = {0:5d}\n'.format(mc_iters))
f.write('# EQ_iters = {0:5d}\n'.format(eq_iters))
f.write('# St. temp = {0:4d}\n'.format(start_temp))
f.write('# End temp = {0:4d}\n'.format(end_temp))
f.write('# Num. temp = {0:3d}\n'.format(n_temp))
f.write('#\n')
f.write('# Temp     Energy     Mag     Sp.Heat    Susc.\n')
f.close()

#compute runtime
start = time.time()

#Start the sims
for i in range(n_temp):
    dipolar_MC1.temp = temp_var[i]
    dipolar_MC1.MC_move()
    f = open(dir+data_file,"a+")
    f.write('{0:.3f}, {1:.4e}, {2:.3f}, {3:.4e}, {4:.4e} \n'.format(dipolar_MC1.temp, dipolar_MC1.avgenergy, dipolar_MC1.netmag, dipolar_MC1.sp_heat, dipolar_MC1.suscep))
    f.close()
    print(dipolar_MC1.temp, dipolar_MC1.avgenergy, dipolar_MC1.netmag, dipolar_MC1.sp_heat, dipolar_MC1.suscep, dipolar_MC1.n_highaccept, dipolar_MC1.n_lowaccept, dipolar_MC1.n_noaccept)
    #save the draw lattice data
    if (display):
        
        fig, ax1 = plt.subplots(figsize=(8,8))
        ax1.set_title('Lattice State at {0:.3f}'.format(dipolar_MC1.temp))
        q1 = ax1.quiver(dipolar_MC1.centers[0,:],dipolar_MC1.centers[1,:],dipolar_MC1.magx,dipolar_MC1.magy,pivot='mid')
        plt.draw()
        plt.savefig(dir+'Lattice_state_'+str(i)+'.png',bbox_inches='tight')
        plt.close()
    
    #save magnetization data
    f2 = open(dir+'MCrun_mag_'+jobID+'_'+str(i)+'.txt','w+')
    f2.write('# Num islands {0:3d} \n'.format(dipolar_MC1.n_isl))
    for i in range(dipolar_MC1.n_isl):
        f2.write('{0:.3f}, {1:.3f} \n'.format(dipolar_MC1.magx[i],dipolar_MC1.magy[i]))

    f2.close()
        
end = time.time()
print('Finished the run in time:',end-start)
   

    
    
