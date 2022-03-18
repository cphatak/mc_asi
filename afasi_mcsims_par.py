#
#
# This set of programs and routines are to perform Monte Carlo
# simulations on the AF-ASI lattice using dipolar coupling as the
# the dominant energy interactions.
#
# Written, CD Phatak, ANL, 21.Feb.2019.
#
# Modified to use the Dipolar_MC Class file.
#
# Modified, CD Phatak, ANL, 10.Jul.2019 - to use multiprocessing for parallel running.

#import necessary modules
import numpy as np
import datetime
import time as time
import matplotlib as mpl
#determine if linux
import platform
    #if (platform.system() == 'Linux'):
mpl.use('Agg') # for linux operation.
from matplotlib import pyplot as plt
from dipolar_MC import Dipolar_MC
import os
import multiprocessing as mp

#------------------------------------------------------------------
#
# Printing out system information.
#
def print_sysinfo():
    
    print('\nPython version  :', platform.python_version())
    print('compiler        :', platform.python_compiler())
    
    print('\nsystem     :', platform.system())
    print('release    :', platform.release())
    print('machine    :', platform.machine())
    print('processor  :', platform.processor())
    print('CPU count  :', mp.cpu_count())
    print('interpreter:', platform.architecture()[0])
    print('\n\n')


#------------------------------------------------------------------
#
# the core function for running an individual MC simulation.
def run_MC(n_run, #job ID number
           a, #Lattice parameter
           s, #island separation
           nx, #num of islands along x
           ny, #num of islands along y
           th, #rotation angle for paired islands
           mc_iters, #number of MC iterations
           eq_iters, #number of equilibriation iterations
           start_temp, #Start temperature
           end_temp, #end temperature
           n_temp, #number of temperature steps
           red_fac, #reduction factor - 'Linear' or 'Geometric'
           load_file, #load magnetization config data for MC runs
           save_file, #to save files in intermediate MC steps. every N steps.
           file_name, #filename for loading the data.
           verbose,
           work_dir, #working directory
           display):
    
    #create the job ID
    jobID = 'run'+str(n_run)
    plt.ioff()
    
    #get current working directory
    cwd = os.getcwd()
    if not work_dir:
        dir = cwd+'/'+jobID+'/'
    else:
        dir = work_dir + '/'+jobID+'/'

    print('Saving to directory - '+dir)
    
    #create directory if it does not exist.
    if not os.path.exists(dir):
        os.makedirs(dir)

    #Set the next nearest neghbors
    max_nn_num = 9
    max_nn_dist = 1.5 * a
    
    #variable for pair flip
    pairflip = False

    #initialize random number
    np.random.seed()
    
    #Compute number of islands.
    #6 islands per motif.
    n_isl = nx * ny * 6
    print('Number of islands:',n_isl)

    #next we initialize the lattice.
    dipolar_MC1 = Dipolar_MC(a = a, s = s, nx = nx, ny = ny, th = th, max_nn_dist = max_nn_dist,
                            max_nn_num = max_nn_num, dir = dir, jobID = jobID, 
                            latt_type='rectangle', init_random = True)#, centers = centers, angles = angles, nn_inds = nn_inds)
    
    #check if we need to load previous data
    if load_file:
        mag_arr = np.genfromtxt(file_name, delimiter=',', skip_header=1)
        dipolar_MC1.magx = mag_arr[:,0]
        dipolar_MC1.magy = mag_arr[:,1]
      
    dipolar_MC1.energy = dipolar_MC1.Latt_Energy(debug=False)
    current_energy = dipolar_MC1.energy
    print("Current Energy:",dipolar_MC1.energy)

    ff = open(dir+'MC_nn_inds.txt','w+')
    ff.write(str(dipolar_MC1.nn_inds))
    ff.close()

    #save the images if display option is set
    if (display):
    
        fig1, ax1 = plt.subplots(figsize=(8,8))
        ax1.set_title('MC sims initial state')
        qq = ax1.quiver(dipolar_MC1.centers[0,:],dipolar_MC1.centers[1,:],dipolar_MC1.magx,dipolar_MC1.magy,pivot='mid')
        #plt.draw()
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
    if (red_fac == 'Linear'):
        temp_var = np.linspace(start_temp,end_temp,n_temp)

    if (red_fac == 'Geometric'):
        temp_var = np.geomspace(start_temp,end_temp,n_temp,endpoint=True)

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
    f.write('# Num Neighbors = {0:2d}\n'.format(max_nn_num))
    f.write('# Rotation angle = {0:.1f}\n'.format(th))
    f.write('# Num isl = {0:4d}\n'.format(n_isl))
    f.write('# MC iters = {0:5d}\n'.format(mc_iters))
    f.write('# EQ_iters = {0:5d}\n'.format(eq_iters))
    f.write('# St. temp = {0:.4e}\n'.format(start_temp))
    f.write('# End temp = {0:.4e}\n'.format(end_temp))
    f.write('# Num. temp = {0:3d}\n'.format(n_temp))
    f.write('#\n')
    f.write('Temp,     Energy,     Mag,    AvgSpin, Sp.Heat,    Susc., UL., ULspin, LowAccept,  HighAccept,   NoAccept\n')
    f.close()

    #compute runtime
    start = time.time()
    
    #break counter
    break_count = 0

    #Start the sims
    for i in range(n_temp):
        dipolar_MC1.temp = temp_var[i]
        if verbose:
            print(i,dipolar_MC1.temp)
        dipolar_MC1.MC_move(verbose=verbose, pairflip=pairflip, save_file=save_file, debug=False) #optional argument verbose.
        f = open(dir+data_file,"a+")
        f.write('{0:.4e}, {1:.4e}, {2:.3f}, {3:.3f}, {4:.4e}, {5:.4e}, {6:.3f}, {7:.3f}, {8:5d}, {9:5d}, {10:5d}\n'.format(
            dipolar_MC1.temp, dipolar_MC1.avgenergy, dipolar_MC1.netmag, dipolar_MC1.sumspin, dipolar_MC1.sp_heat, 
            dipolar_MC1.suscep, dipolar_MC1.ul, dipolar_MC1.ul2, dipolar_MC1.n_lowaccept, dipolar_MC1.n_highaccept, 
            dipolar_MC1.n_noaccept))
        f.close()
        if verbose:
            print(dipolar_MC1.temp, dipolar_MC1.avgenergy, dipolar_MC1.netmag, dipolar_MC1.sumspin, dipolar_MC1.sp_heat, dipolar_MC1.suscep, dipolar_MC1.ul, dipolar_MC1.ul2, dipolar_MC1.n_highaccept, dipolar_MC1.n_lowaccept, dipolar_MC1.n_noaccept)
        #save the draw lattice data
        if (display):
        
            fig, ax1 = plt.subplots(figsize=(8,8))
            ax1.set_title('Lattice State at {0:.3f}'.format(dipolar_MC1.temp))
            q1 = ax1.quiver(dipolar_MC1.centers[0,:],dipolar_MC1.centers[1,:],dipolar_MC1.magx,dipolar_MC1.magy,pivot='mid')
            #plt.draw()
            plt.savefig(dir+'Lattice_state_'+str(i)+'.png',bbox_inches='tight')
            plt.close()
    
        #save magnetization data
        f2 = open(dir+'MCrun_mag_'+jobID+'_'+str(i)+'.txt','w+')
        f2.write('# Num islands {0:3d} \n'.format(dipolar_MC1.n_isl))
        for im in range(dipolar_MC1.n_isl):
            f2.write('{0:.3f}, {1:.3f} \n'.format(dipolar_MC1.magx[im],dipolar_MC1.magy[im]))

        f2.close()
        
        #save spin correlation data
        f3 = open(dir+'MCrun_spincorr_'+jobID+'_'+str(i)+'.txt','w+')
        f3.write('# Num islands {0:3d} \n'.format(dipolar_MC1.n_isl))
        f3.write('# Max sp_cr islands {0:3d} \n'.format(dipolar_MC1.max_spcr_num))
        f3.write('# Temperature {0:.4e} \n'.format(dipolar_MC1.temp))
        for isc in range(dipolar_MC1.max_spcr_num):
            f3.write('{0:.3f}, {1:.3f}, {2:.3f}, {3:.3f}, {4:.3f} \n'.format(dipolar_MC1.sisj[0,isc],dipolar_MC1.sisj[1,isc],dipolar_MC1.sisj[2,isc],dipolar_MC1.sisj[3,isc],dipolar_MC1.sisj[4,isc]))
        
        f3.close()
        
        
        #Update current energy
        current_energy = dipolar_MC1.avgenergy
        
        #Now we are checking for activating pair flip.
        #If the total accepted changes are zero and net magnetization is zero, then we activate it.
        #if (((dipolar_MC1.n_highaccept + dipolar_MC1.n_lowaccept) == 0) and (dipolar_MC1.netmag == 0)):
        #    pairflip = True
        #    print('Activated Pair Flip at temp {0:.4e}'.format(dipolar_MC1.temp))
        
        #check if the MC simulation has reached final state.
        if (((dipolar_MC1.n_highaccept+dipolar_MC1.n_lowaccept) == 0) and ((dipolar_MC1.avgenergy - current_energy) == 0)):
            break_count += 1
            if (break_count > 5):
                break
        
        
    end = time.time()
    print('Finished the run in time:',end-start)

    return 1



#------------------------------------------------------------------
    
    
