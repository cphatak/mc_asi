#!/usr/bin/python
#
# Python Class file for Dipolar MC
#
# The purpose of this class file is to create an object
# for performing Monte Carlo Simulations that are based
# on dipolar energy interactions as the main energy term
#
# Written by CD Phatak, ANL, 02.May.2019.

#import necessary modules
import numpy as np
import scipy.constants as physcon

class Dipolar_MC(object):

    def __init__(self,
                 centers = 0, #array of centers of each islands
                 angles = 0, #array of angles for each island
                 nn_inds = 0, #array of indices for nearest neighbors
                 max_nn_num = 1, #max. number of neighbors to consider
                 max_nn_dist = 1, #max. distance for nearest neighbor
                 verbose=False):
        #This function is for initializing the dipolar_MC object
        #The above parameters can be set while initializing and others can
        #be set below.
        
        #parameters for describing lattice
        self.centers = centers
        self.angles = angles
        self.nn_inds = nn_inds
        self.max_nn_num = max_nn_num
        self.max_nn_dist = max_nn_dist
        
        #MC simulation parameters - default values
        #they can be changed after definition.
        self.mc_iters = 1000 #total MC iters
        self.eq_iters = 0 #number of iters for equilbriation before computing
        self.temp = 10 #dimensionless temperature parameter
        self.mult_fac = physcon.mu_0*1e-9/physcon.k
        
        #other derived parameters
        temp, n_isl = centers.shape
        self.n_isl = n_isl #numner of islands
        #magnetization parameter
        self.magx = np.cos(np.deg2rad(angles))
        self.magy = np.sin(np.deg2rad(angles))
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
        
        #print output message
        print("Created the Dipolar_MC class.")
        print("Please run self.Calc_Energy to update self energy")


    #------------------------------------------------------------------
    #
    # Calc_Energy function for AFASI
    #
    def Calc_Energy(self, debug=False):
        
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
                    
                    sj_rji = (self.centers[0,j]-self.centers[0,i])*self.magx[j] 
                    + (self.centers[1,j]-self.centers[1,i])*self.magy[j]
                    
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
    # MC_move function to actual run the MC simulation.
    #
    # This function will take the input parameters for the MC simulation
    # then for number of MC iters, randomly select a spin site and flip its
    # direction, calculate energy and check whether to accept it or not.
    # It will also calculate the thermodynamic parameters such as sp. heat,
    # and susceptibility during the MC iters for a given temperature.
    #
    def MC_move(self,verbose=False,debug=False):
    
        #initialize arrays for holding various quantities
        avg_en = 0.0
        avg_en2 = 0.0
        avg_mag = 0.0
        avg_mag2 = 0.0
        
        #get current energy
        self.energy = self.Calc_Energy()
    
        for nn in range(self.mc_iters):
            if (verbose):
                if (np.mod(nn,10) == 0):
                    print(nn)
            for ii in range(self.n_isl):
                # pick a random site in the lattice.
                site = np.random.randint(0,self.n_isl)
                if (debug):
                    print(site)
    
                #change the magnetization
                self.magx[site] *= (-1)
                self.magy[site] *= (-1)
    
                #calculate the energy
                new_energy = self.Calc_Energy()
    
                #difference in energy
                dE = new_energy - self.energy
                if (debug):
                    print(dE)
    
                #check if we should accept this energy or not
                if (dE < 0):
                    self.n_lowaccept += 1
                    self.energy = new_energy
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
                        self.energy = new_energy
                        if (debug):
                            print('High accept')
                            
                    else:
                        #we do not accept the change
                        self.magx[site] *= (-1)
                        self.magy[site] *= (-1)
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
                
        #we are out of the MC loop.
        #collect all the data.
        cn = 1.0/((self.mc_iters-self.eq_iters)*self.n_isl)
        self.avgenergy = avg_en*cn
        self.netmag = avg_mag*cn
        self.sp_heat = (avg_en2*cn - avg_en*avg_en*cn**2)/self.temp**2
        self.suscep = (avg_mag2*cn - avg_mag*avg_mag*cn**2)/self.temp
        
        #print some output
        print("MC iter complete at temp:",self.temp)
        return 1


## MAIN ##
if __name__ == '__main__':
    print(" Class definition for Dipolar_MC Class.")
