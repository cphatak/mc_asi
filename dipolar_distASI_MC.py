#!/usr/bin/python
#
# Python Class file for Dipolar MC
#
# The purpose of this class file is to create an object
# for performing Monte Carlo Simulations that are based
# on dipolar energy interactions as the main energy term
#
# Written by CD Phatak, ANL, 23.Oct.2020.

#import necessary modules
import numpy as np
#import scipy.constants as physcon
from scipy.spatial import cKDTree as sp_cKDTree


class Dipolar_distASI_MC(object):

    def __init__(self,
                 a = 350, # lattice parameter
                 nx = 1, # repeat along x
                 ny = 1, # repeat along y
                 max_nn_dist = 500, # max. distance of nearest neighbors
                 max_nn_num = 9, # max. number of nearest neighbors
                 man_fname = 'quarry.txt', #filename for manual input of lattice info.
                 dir = '/', #folder name
                 jobID = 'run1', #job ID
                 init_random = True, #initial magnetization of the lattice
                 verbose = False):
        #This function is for initializing the dipolar_MC object
        #The above parameters can be set while initializing and others can
        #be set below.
        
        #parameters for describing lattice
        #multiplier for each motif depends on type of lattice.
        #we are reading the manual file
        data = np.genfromtxt(man_fname,delimiter=',',skip_header=3)
        n_isl, cols = data.shape
        self.n_isl = n_isl
        self.centers = np.zeros([2,self.n_isl])
        self.angles = np.zeros([self.n_isl])
        self.centers[0,:] = data[:,3]
        self.centers[1:,] = data[:,4]
        self.angles[:] = data[:,5]
    
        #set nearest neighbor parameters.
        self.max_nn_num = max_nn_num
        self.max_nn_dist = max_nn_dist
        
        #folder location.
        self.dir = dir
        self.jobID = jobID
        
        
        #now to use the cKDTree method
        comb_xy = self.centers.transpose()
        p_list = list(comb_xy)
        self.nn_inds = self.do_kdtree(comb_xy, p_list, max_nn_num+1, max_nn_dist)
        #MC simulation parameters - default values
        #they can be changed after definition.
        self.mc_iters = 1000 #total MC iters
        self.eq_iters = 0 #number of iters for equilbriation before computing

        #We will calculate the energy normalized by D where D = mu_0*mu^2/4/pi
        #The temperature scale is then in units of D/k_b. 
        self.temp = 1
        #self.mult_fac = physcon.mu_0*1e-9/physcon.k
        self.mult_fac = 1.0 # multiplication factor for partition function calculation. 
        
        #magnetization parameter
        self.magx = np.cos(np.deg2rad(self.angles))
        self.magy = np.sin(np.deg2rad(self.angles))

        #randomize the magnetization
        if init_random:
            for nn in range(self.n_isl):
                mult = 1
                if (np.random.random() <= 0.5):
                    mult = -1
                self.magx[nn] *= mult
                self.magy[nn] *= mult
        
        #compute and store the distance map.
        self.distmap = np.zeros([self.n_isl, self.n_isl])
        for ii in range(self.n_isl):
            for jj in range(self.n_isl):
                self.distmap[ii,jj] = np.sqrt((self.centers[0,ii]-self.centers[0,jj])**2 + (self.centers[1,ii]-self.centers[1,jj])**2)
        
        #parameters for storing various values
        self.n_lowaccept = 0
        self.n_highaccept = 0
        self.n_noaccept = 0
        self.energy = 0
        self.avgenergy = 0
        self.netmag = np.sqrt(np.sum(self.magx)**2 + np.sum(self.magy)**2)
        self.sp_heat = 0
        self.suscep = 0
        self.ul = 0
        
        #print output message
        print("Created the Dipolar_MC class.")
        print("Please run self.Latt_Energy to update self energy")

    #------------------------------------------------------------------
    #
    # Function using the kdtree algorithm to find the nearest neighbors.
    # USing the method described here - 
    # https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates
    #
    def do_kdtree(self, combined_x_y_arrays, points, max_nn, max_dd):
        mytree = sp_cKDTree(combined_x_y_arrays)
        dist, indexes = mytree.query(points, k=max_nn, distance_upper_bound=max_dd)
        return indexes
    
    
    
    #------------------------------------------------------------------
    #
    # Latt_Energy function for Rhombille Latt
    #
    # Computes the total energy of the entire lattice
    #
    def Latt_Energy(self, debug=False):
        
        # Energy variable
        tot_energy = 0
        count = 0
    
        # loop over each island and compute the neighbors
        for i in range(self.n_isl):
            for cnt in range(self.max_nn_num-1):
                j = self.nn_inds[i,cnt+1].astype('int')
                if ((i != j) and (j != self.n_isl)):
                    
                    si_sj = self.magx[i]*self.magx[j] + self.magy[i]*self.magy[j]
                    
                    r_ij = self.distmap[i,j]
                    
                    si_rij = (self.centers[0,i]-self.centers[0,j])*self.magx[i] 
                    + (self.centers[1,i]-self.centers[1,j])*self.magy[i]
                    
                    sj_rji = (self.centers[0,i]-self.centers[0,j])*self.magx[j]
                    + (self.centers[1,i]-self.centers[1,j])*self.magy[j]
                    
                    temp = (((si_sj)/r_ij**3) - ((3.0*si_rij*sj_rji)/r_ij**5))
                    tot_energy +=  temp
                    if debug:
                        print(i,j,r_ij,temp,tot_energy)
                        
                    count += 1
    
        #return total energy
        #self.energy = tot_energy/2.0
        if debug:
            print(count)
        return tot_energy/2.0

    #------------------------------------------------------------------
    #
    # Calc_del_Energy function for Rhombille Latt.
    #
    # Computes the energy of a given site in the lattice.
    #
    def Calc_del_Energy(self, site, pairflip = False, debug=False):
        
        # Energy variable
        site_energy = 0
        count = 0
    
        # compute the neighbors for the given site.
        if pairflip:
            site_indices = [site, self.nn_inds[site,1]]
            for i in site_indices:                
                for cnt in range(self.max_nn_num-1):
                    j = self.nn_inds[i,cnt+2].astype('int')
                    if ((i != j) and (j != self.n_isl)):
                        
                        si_sj = self.magx[i]*self.magx[j] + self.magy[i]*self.magy[j]
                        
                        r_ij = self.distmap[i,j]
                        
                        si_rij = (self.centers[0,i]-self.centers[0,j])*self.magx[i] 
                        + (self.centers[1,i]-self.centers[1,j])*self.magy[i]
                        
                        sj_rji = (self.centers[0,i]-self.centers[0,j])*self.magx[j]
                        + (self.centers[1,i]-self.centers[1,j])*self.magy[j]
                        
                        temp = (((si_sj)/r_ij**3) - ((3.0*si_rij*sj_rji)/r_ij**5))
                        site_energy +=  temp
                        if debug:
                            print(i,j,r_ij,temp,site_energy)
                            
                        count += 1
        else:
            i = site
            for cnt in range(self.max_nn_num-1):
                j = self.nn_inds[i,cnt+1].astype('int')
                if ((i != j) and (j != self.n_isl)):
                    
                    si_sj = self.magx[i]*self.magx[j] + self.magy[i]*self.magy[j]
                    
                    r_ij = self.distmap[i,j]
                    
                    si_rij = (self.centers[0,i]-self.centers[0,j])*self.magx[i] 
                    + (self.centers[1,i]-self.centers[1,j])*self.magy[i]
                    
                    sj_rji = (self.centers[0,i]-self.centers[0,j])*self.magx[j]
                    + (self.centers[1,i]-self.centers[1,j])*self.magy[j]
                    
                    temp = (((si_sj)/r_ij**3) - ((3.0*si_rij*sj_rji)/r_ij**5))
                    site_energy +=  temp
                    if debug:
                        print(i,j,r_ij,temp,site_energy)
                        
                    count += 1
        
        #return total energy
        #self.energy = tot_energy/2.0
        if debug:
            print(count)
        return site_energy


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
    def MC_move(self, pairflip = False, save_file = 1000,
                verbose=False,debug=False):
    
        #initialize arrays for holding various quantities
        avg_en = 0.0
        avg_en2 = 0.0
        avg_mag = 0.0
        avg_mag2 = 0.0
        avg_mag4 = 0.0
        
        #reset the counters for accepted values.
        self.n_lowaccept = 0
        self.n_highaccept = 0
        self.n_noaccept = 0
        
        #get current energy
        self.energy = self.Latt_Energy()
    
        for nn in range(self.mc_iters):
            if (verbose):
                if (np.mod(nn,100) == 0):
                    print(nn,',',end='')
            for ii in range(self.n_isl):
                # pick a random site in the lattice.
                site = np.random.randint(0,self.n_isl)
                if (debug):
                    print(site)
    
                #change the magnetization
                self.magx[site] *= (-1)
                self.magy[site] *= (-1)
                
                if pairflip:
                    pair_site = self.nn_inds[site,1]
                    self.magx[pair_site] *= (-1)
                    self.magy[pair_site] *= (-1)
    
                #calculate the change in energy
                dE = self.Calc_del_Energy(site, pairflip = pairflip)*2.0
                if (debug):
                    print(dE)
    
                #check if we should accept this energy or not
                if (dE < 0):
                    self.n_lowaccept += 1
                    self.energy += dE
                    if (debug):
                        print('Low accept')
                
                if (dE > 0):
                    #we check if we should accept the high energy change
                    rnum = np.random.random_sample()
                    part_fun = np.exp(-dE*self.mult_fac/self.temp)
                    if (debug):
                        print(rnum,part_fun)
                    
                    if (rnum < part_fun):
                        self.n_highaccept += 1
                        self.energy += dE
                        if (debug):
                            print('High accept')
                            
                    else:
                        #we do not accept the change
                        self.magx[site] *= (-1)
                        self.magy[site] *= (-1)
                        if pairflip:
                            self.magx[pair_site] *= (-1)
                            self.magy[pair_site] *= (-1)
                        self.n_noaccept += 1
                        if (debug):
                            print('No accept')
                    
                #Next we start computing various thermo. terms
                self.netmag = np.sqrt(np.sum(self.magx)**2 + np.sum(self.magy)**2)
                
                if (nn >= self.eq_iters):
                    avg_en += self.energy
                    avg_en2 += self.energy**2
                    avg_mag += self.netmag
                    avg_mag2 += self.netmag**2
                    avg_mag4 += self.netmag**4
                
                #Save the file if needed
                if (np.mod(nn,save_file) == 0):
                    f1 = open(self.dir+'mag_data_'+self.jobID+'_temp_'+str(self.temp)+'_MC_iter_'+str(nn)+'.txt','w+')
                    f1.write('# Num islands {0:3d} \n'.format(self.n_isl))
                    for i in range(self.n_isl):
                        f1.write('{0:.3f}, {1:.3f} \n'.format(self.magx[i],self.magy[i]))

                    f1.close()
                
        #we are out of the MC loop.
        #collect all the data.
        cn = 1.0/((self.mc_iters-self.eq_iters)*self.n_isl)
        self.avgenergy = avg_en*cn
        self.netmag = avg_mag*cn
        self.sp_heat = (avg_en2*cn - avg_en*avg_en*cn**2)/self.temp**2
        self.suscep = (avg_mag2*cn - avg_mag*avg_mag*cn**2)/self.temp
        self.ul = 1.0 - (avg_mag4*cn/(3.0*(avg_mag2*cn)**2))
        
        #print some output
        if verbose:
            print("\n MC iter complete at temp:",self.temp)
        return 1


## MAIN ##
if __name__ == '__main__':
    print(" Class definition for Dipolar_MC Class.")
