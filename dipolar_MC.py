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
from numpy.core.numerictypes import maximum_sctype
import scipy.constants as physcon
from scipy.spatial import cKDTree as sp_cKDTree


class Dipolar_MC(object):

    def __init__(self,
                 a = 350, # lattice parameter
                 s = 120, # island separation
                 nx = 1, # repeat along x
                 ny = 1, # repeat along y
                 th = 30.0, # rotation angle of paired islands
                 max_nn_dist = 500, # max. distance of nearest neighbors
                 max_nn_num = 9, # max. number of nearest neighbors
                 latt_type = 'slanted', #slated or rectangular lattice.
                 dir = '/', #folder name
                 jobID = 'run1', #job ID
                 init_random = True, #initial magnetization of the lattice
                 verbose = False):
        #This function is for initializing the dipolar_MC object
        #The above parameters can be set while initializing and others can
        #be set below.
        
        #parameters for describing lattice
        self.n_isl = nx * ny * 6
        self.centers = np.zeros([2,self.n_isl])
        self.angles = np.zeros([self.n_isl])
        self.max_nn_num = max_nn_num
        self.max_nn_dist = max_nn_dist
        self.max_spcr_num = max_nn_num * 2
        
        #folder location.
        self.dir = dir
        self.jobID = jobID
        
        #initialize the lattice
        res = self.init_afasi_latt(a = a, s = s, nx = nx, ny = ny, th = th, latt_type = latt_type)
        
        #now to use the cKDTree method
        comb_xy = self.centers.transpose()
        #p_list = list(comb_xy)
        p_list = comb_xy
        self.nn_inds = self.do_kdtree(comb_xy[:], p_list[:], max_nn_num+1, max_nn_dist)
        self.nn_inds[self.nn_inds == self.n_isl] = 0
        #get indices for computing spin correlation for n=2*max_nn_num
        self.spcr_nn_inds = self.do_kdtree(comb_xy[:],p_list[:],max_nn_num*2+1, max_nn_dist*2)
        self.spcr_nn_inds[self.spcr_nn_inds == self.n_isl] = 0
        
        
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
        self.sumspin = np.sum(self.magx) + np.sum(self.magy)
        self.sisj = np.zeros([5,self.max_spcr_num])
        self.sp_heat = 0
        self.suscep = 0
        self.ul = 0
        self.ul2 = 0
        
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
    # Init lattice function for AFASI
    #
    # This function will initialize the lattice for given set of lattice parameters
    # and num of islands along x and y. The output will be an array centers consisting
    # of positions of each island (x,y), and an array angles with angle of each island
    # for magnetization, and an array nn_inds consisting of nearest neighbor indices
    # for each island to be considered for dipolar interactions.
    #
    # Updated to include a rotational angle for the paired islands.
    #
    def init_afasi_latt(self, a = 350, s = 120, nx = 1, ny = 1,
                        th = 30, latt_type= 'slanted',
                        verbose = False, debug = False):
    
        
        thr = np.deg2rad(th)
        r_mat = np.array([[np.cos(thr),-np.sin(thr)],[np.sin(thr),np.cos(thr)]])
    
        #loop over number of islands
        count = 0
        ext_cens = np.zeros([2,6])
        for nn in range(6):
            ext_cens[0,0] = s/2
            ext_cens[1,0] = 0
            ext_cens[0,1] = -s/2
            ext_cens[1,1] = 0
            ext_cens[0,2] = -s/4
            ext_cens[1,2] = np.sqrt(3)*s/4
            ext_cens[0,3] = s/4
            ext_cens[1,3] = -np.sqrt(3)*s/4
            ext_cens[0,4] = -s/4
            ext_cens[1,4] = -np.sqrt(3)*s/4
            ext_cens[0,5] = s/4
            ext_cens[1,5] = np.sqrt(3)*s/4
        
        ext_cens = np.matmul(r_mat,ext_cens)
            
        
        for i in range(nx):
            for j in range(ny):
                
                if (latt_type == 'rectangle'):
                            if (np.mod(j,2) != 0):
                                i -= 1
        
                self.centers[0,count] = 0 + ext_cens[0,0] + 2*a*i + a*j
                self.centers[1,count] = 0 + ext_cens[1,0] + np.sqrt(3)*a*j
                self.angles[count] = 90 + th
                
                self.centers[0,count+1] = 0 + ext_cens[0,1] + 2*a*i + a*j
                self.centers[1,count+1] = 0 + ext_cens[1,1] + np.sqrt(3)*a*j
                self.angles[count+1] = 90 + th
                
                self.centers[0,count+2] = a/2 + ext_cens[0,2] + 2*a*i + a*j
                self.centers[1,count+2] = np.sqrt(3)/2*a + ext_cens[1,2] + np.sqrt(3)*a*j
                self.angles[count+2] = 30 + th
                
                self.centers[0,count+3] = a/2 + ext_cens[0,3] + 2*a*i + a*j
                self.centers[1,count+3] = np.sqrt(3)/2*a + ext_cens[1,3] + np.sqrt(3)*a*j
                self.angles[count+3] = 30 + th
                
                self.centers[0,count+4] = (-a/2) + ext_cens[0,4] + 2*a*i + a*j
                self.centers[1,count+4] = np.sqrt(3)/2*a + ext_cens[1,4] + np.sqrt(3)*a*j
                self.angles[count+4] = 150 + th
                
                self.centers[0,count+5] = (-a/2) + ext_cens[0,5] + 2*a*i + a*j
                self.centers[1,count+5] = np.sqrt(3)/2*a + ext_cens[1,5] + np.sqrt(3)*a*j
                self.angles[count+5] = 150 + th
                
                count += 6
        
        
        # # Counter for tracking islands     
        # count = 0
        # for i in range(nx):
        #     for j in range(ny):
        #         if (latt_type == 'rectangle'):
        #             if (np.mod(j,2) != 0):
        #                 i -= 1
                
        #         if debug:
        #             print(i,j)
        #         #horizontal islands
        #         self.angles[count] = 0
        #         self.angles[count+1] = 0
        #         self.centers[0,count] = i*2*a - a/2 + j*a
        #         self.centers[1,count] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 + s/2
        #         self.centers[0,count+1] = i*2*a - a/2 + j*a
        #         self.centers[1,count+1] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 - s/2
        #         #first set of rotated islands
        #         self.angles[count+2] = 120
        #         self.angles[count+3] = 120
        #         self.centers[0,count+2] = i*2*a + a/2 + j*a + np.sqrt(3)/4*s
        #         self.centers[1,count+2] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 + s/4
        #         self.centers[0,count+3] = i*2*a + a/2 + j*a - np.sqrt(3)/4*s
        #         self.centers[1,count+3] = j*np.sqrt(3)*a - a*np.sqrt(3)/4 - s/4
        #         #second set of rotated islands
        #         self.angles[count+4] = 60
        #         self.angles[count+5] = 60
        #         self.centers[0,count+4] = i*2*a + j*a - np.sqrt(3)*s/4
        #         self.centers[1,count+4] = j*np.sqrt(3)*a + a*np.sqrt(3)/4 + s/4
        #         self.centers[0,count+5] = i*2*a + j*a + np.sqrt(3)*s/4
        #         self.centers[1,count+5] = j*np.sqrt(3)*a + a*np.sqrt(3)/4 - s/4
    
        #         #increment count
        #         count += 6
    
        return 1

    #------------------------------------------------------------------
    #
    # Latt_Energy function for AFASI
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
    # Calc_del_Energy function for AFASI
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
                        
                        sj_rji = (self.centers[0,j]-self.centers[0,i])*self.magx[j] 
                        + (self.centers[1,j]-self.centers[1,i])*self.magy[j]
                        
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
                    
                    sj_rji = (self.centers[0,j]-self.centers[0,i])*self.magx[j] 
                    + (self.centers[1,j]-self.centers[1,i])*self.magy[j]
                    
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
        avg_sumspin = 0.0
        avg_sumspin2 = 0.0
        avg_sumspin4 = 0.0
        #sisj_arr = np.zeros([4,self.max_spcr_num])
        sisj_distvals = np.zeros([self.max_spcr_num])
        sisj_corr = [[] for _ in range(self.max_spcr_num)]
        sp_corr_cnt = 0.0
        
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
                #if (debug):
                #    print(site)
    
                #change the magnetization
                self.magx[site] *= (-1)
                self.magy[site] *= (-1)
                
                if pairflip:
                    pair_site = self.nn_inds[site,1]
                    self.magx[pair_site] *= (-1)
                    self.magy[pair_site] *= (-1)
    
                #calculate the change in energy
                dE = self.Calc_del_Energy(site, pairflip = pairflip)*2.0
                #if (debug):
                #    print(dE)
    
                #check if we should accept this energy or not
                if (dE < 0):
                    self.n_lowaccept += 1
                    self.energy += dE
                #    if (debug):
                #        print('Low accept')
                
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
                self.sumspin = np.sum(self.magx) + np.sum(self.magy)
                
                if (nn >= self.eq_iters):
                    avg_en += self.energy
                    avg_en2 += self.energy**2
                    avg_mag += self.netmag
                    avg_mag2 += self.netmag**2
                    avg_mag4 += self.netmag**4
                    avg_sumspin += self.sumspin
                    avg_sumspin2 += self.sumspin**2
                    avg_sumspin4 += self.sumspin**4
                    
                    #compute the spin correlation
                    #for ii in range(self.n_isl):
                    #    for cc in range(self.max_spcr_num):
                    #    
                    #        jj = self.spcr_nn_inds[ii,cc].astype('int')
                    #        
                    #        if ((jj != self.n_isl)):
                    #            si_mod = np.sqrt(self.magx[ii]**2 + self.magy[ii]**2)
                    #            sj_mod = np.sqrt(self.magx[jj]**2 + self.magy[jj]**2)
                    #            si_sj = self.magx[ii]*self.magx[jj] + self.magy[ii]*self.magy[jj]
                    #            sisj_arr[0,cc] = self.distmap[ii,jj]
                    #            sisj_arr[1,cc] += (si_sj)
                    #            sisj_arr[2,cc] += (np.sqrt(self.magx[ii]**2 + self.magy[ii]**2))
                    #            sisj_arr[3,cc] += (np.sqrt(self.magx[jj]**2 + self.magy[jj]**2))
                    #            
                    #            sp_corr_cnt += 1.0
                    #            #if debug:
                    #            #    print('Spin Corr values:',si_sj, si_mod, sj_mod)
                    #for icor in range(self.n_isl):
                    for nth in range(self.max_spcr_num):
                        #jcor = self.spcr_nn_inds[:][nth]
                        #jcor[jcor >= self.n_isl] = self.n_isl-1
                        #if (jcor != self.n_isl):
                        sisj_distvals[nth] = self.distmap[self.n_isl//2,(self.spcr_nn_inds[self.n_isl//2][nth])]
                        c_temp = self.magx[:]*self.magx[self.spcr_nn_inds[:,nth]] + self.magy[:]*self.magy[self.spcr_nn_inds[:,nth]]
                        
                        #sisj_distvals[nth] = self.distmap[0,jcor[nth]]
                        #c_temp = self.magx[:]*self.magx[jcor] + self.magy[:]*self.magy[jcor]
                        
                        sisj_corr[nth].append(c_temp[:])
                       

                                                                
                            
                
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
        self.sumspin = avg_sumspin*cn
        self.sp_heat = (avg_en2*cn - avg_en*avg_en*cn**2)/self.temp**2
        self.suscep = (avg_mag2*cn - avg_mag*avg_mag*cn**2)/self.temp
        self.ul = 1.0 - (avg_mag4*cn/(3.0*(avg_mag2*cn)**2))
        self.ul2 = 1.0 - (avg_sumspin4*cn/(3.0*(avg_sumspin2*cn)**2))
        #self.sisj[0,:] = sisj_arr[0,:]
        #self.sisj[1,:] = sisj_arr[1,:]*cn/self.max_spcr_num
        #self.sisj[2,:] = sisj_arr[2,:]*cn/self.max_spcr_num
        #self.sisj[3,:] = sisj_arr[3,:]*cn/self.max_spcr_num
        #self.sisj[4,:] = sisj_arr[1,:]*cn/self.max_spcr_num - sisj_arr[2,:]*sisj_arr[3,:]*cn**2/self.max_spcr_num**2
        #print(c_temp.shape)
        #print(sum(sum(sisj_corr[0]))/len(sisj_corr[0])/self.n_isl)
        for nth2 in range(self.max_spcr_num):
            self.sisj[0,nth2] = sisj_distvals[nth2]
            self.sisj[1,nth2] = sum(sum(sisj_corr[nth2]))/len(sisj_corr[nth2])/self.n_isl
        
        #self.sisj[4,:] = self.sisj[1,:] - self.sisj[2,:]*self.sisj[3,:]
        
        
        #print some output
        if verbose:
            print("\n MC iter complete at temp:",self.temp)

        return 1
#------------------------------------------------------------------
    #
    # spin_corr function.
    #
    # This function will take the input parameters of the mag array, indices
    # of nearest neighbors and which NN corelation we want to compute, and return
    #

## MAIN ##
if __name__ == '__main__':
    print(" Class definition for Dipolar_MC Class.")
