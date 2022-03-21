"""
Python script to generate the distorted square ASI
lattice with option of setting the rotation angle 
for each island to convert to pin-wheel like lattice

@author, CD Phatak, 05.Mar.2022

"""

#import necessary modules
import numpy as np
import scipy.constants as physcon
import matplotlib.pyplot as plt
import subprocess
import datetime
import time
import math

#arguments - can be passed later on and edited using arg_parse.
#verbose = False
#work_dir = "/home/cphatak/micromagnetic_sims/dist_asi/"
#jobID = '001'
#mx3_fname = 'gen_dist_ASI_latt'
#latt_info_name = 'nondist_ASI_latt_info_'
#create_mx3file = False
#max_row_isl = 10

"""
Function to generate the lattice information and optionally
save the mx3 file.

@author CD Phatak, ANL, 19.Mar.2022.
"""

def gen_dist_asi(l = 200, # nm - length of rectangle part of island
                 w = 100, # nm - width of the island
                 thk = 10, # nm - thickness of islands
                 s1 = 250, # nm - separation parameter 1
                 s2 = 300, # nm - separation parameter 2
                 isl_rot = 0, # rotation angle of island in deg.
                 d_period = 1, # distortion period
                 max_row_isl = 10, # max. no. of islands along x/y
                 work_dir = '/home/cphatak', #working directory
                 latt_info_name = 'dist_ASI_latt_', #basefilename
                 verbose = False,
                 save_fig = True,
                 create_mx3file = False,
                 jobID = '001'):
    
                 
                 
                


    #spin ice geometric parameters
    diam = w #diameter of the end circle
    radius = diam/2

    #lattice size
    xisl_ny = max_row_isl
    xisl_nx = xisl_ny-1
    yisl_nx = xisl_ny
    yisl_ny = xisl_nx

    # number of islands
    num_isl = xisl_nx*xisl_ny + yisl_nx*yisl_ny

    #define the centers
    centers = np.zeros([num_isl,2])
    angles = np.zeros([num_isl])
    mag = centers

    # loop through x and y to create islands
    # initially though all x-islands
    shift_arr = l + (s1+s2)/2
    count = 0
    for iy in range(xisl_ny):
        if np.mod(iy//d_period,2)==0:
            xs = -(l+s1)/2
        else:
            xs = -(l+s2)/2
        xcen = xs
        ycen = iy*shift_arr
        for ix in range(xisl_nx):
            if np.mod(ix//d_period,2)==0:
                xcen += (l+s1)
            else:
                xcen += (l+s2)
            centers[count,:] = [xcen,ycen]
            angles[count] = 0+isl_rot
            count+=1
            if verbose:
                print('Horizontal island: ',count,ix,iy,xcen,ycen)
        

    #next all vertical islands
    for ix in range(yisl_nx):
        if np.mod(ix//d_period,2)==0:
            ys = -(l+s2)/2
        else:
            ys = -(l+s1)/2
        ycen = ys
        xcen = ix*shift_arr
        for iy in range(yisl_ny):
            if np.mod(iy//d_period,2)==0:
                ycen += (l+s2)
            else:
                ycen += (l+s1)
            centers[count,:] = [xcen,ycen]
            angles[count] = 90+isl_rot
            count+=1
            if verbose:
                print('Vertical island: ',count,ix,iy,xcen,ycen)

    #save display
    if save_fig:

        fig,ax = plt.subplots()
        ax.quiver(centers[:,0],centers[:,1],np.cos(np.deg2rad(angles[:])),np.sin(np.deg2rad(angles[:])),pivot='mid')
        ax.set_aspect('equal')
        plt.savefig(work_dir+latt_info_name+jobID+'.png')

    #save lattice coordinates

    f1 = open(work_dir+latt_info_name+'.txt','w+')
    f1.write('# Lattice information: \n')
    f1.write('# Num islands {0:3d} \n'.format(num_isl))
    f1.write('# L (nm), s1 (nm), s2 (nm), cen_x (nm), cen_y (nm), ang, mag_x, mag_y\n')
    for i in range(num_isl):
        f1.write('{0:.3f}, {1:.3f}, {2:.3f}, {3:.3f}, {4:.3f}, {5:.3f}, {6:.3f}, {7:.3f} \n'.format(l, s1, s2, centers[i,0], centers[i,1], angles[i], np.cos(np.deg2rad(angles[i])), np.sin(np.deg2rad(angles[i]))))

    f1.close()

    if create_mx3file:

        #length multiplier to be in m.
        len_mult = 1.0e-9

        #set Permalloy parameters
        Msat = 8e+5
        Aex = 13e-12
        alpha= 0.05
        lex=np.sqrt(Aex/(0.5*physcon.mu_0*Msat**2)) #exchange length
        print("Exchange length: ",lex)

        #world size definition
        max_cent = (l+diam+max(s1,s2))*len_mult
        gridx = max(xisl_nx,xisl_ny,yisl_nx,yisl_ny) * max_cent
        gridy = gridx
        gridz = thk 
        lmax = 0.52*lex
        cell = 2e-9
        cell_size = min(cell,lmax)
        print(cell_size)
        #determine the simulation size closest to power of 2
        Nx = 2**math.ceil(math.log2(gridx//cell_size))
        Ny = 2**math.ceil(math.log2(gridy//cell_size))
        Nz = 2**math.ceil(math.log2(gridz//cell_size))
        #for now use just the closest round value.
        Nx = math.ceil(gridx/cell_size)
        Ny = math.ceil(gridy/cell_size)
        Nz = math.ceil(gridz/cell_size)
        print(Nx, Ny, Nz)
        shift_cenx = shift_arr/2*(xisl_ny-1)*len_mult#gridx/2#*cell_size
        shift_ceny = shift_arr/2*(yisl_nx-1)*len_mult#gridy/2#*cell_size

        # #write the file
        d = datetime.datetime.now()
        f = open(work_dir+latt_info_name+'.mx3','w')
        f.write('// Mumax runtime file for generating distorted sq ASI.\n')
        f.write('// Part of JOB:'+jobID+'\n')
        f.write('// Created: C. Phatak, ANL, ')
        f.write('{:%Y-%m-%d %H:%M:%S}'.format(d))
        f.write('\n')
        f.write('setgridsize({0:3.0f}, {1:3.0f}, {2:3.0f})\n'.format(Nx,Ny,Nz))
        f.write('setcellsize({0:.2e}, {1:.2e}, {2:.2e})\n'.format(cell_size,cell_size,cell_size))
        f.write('\n')
        f.write('\n// Magnetization parameters\n')
        f.write('Msat = {:.2e}\n'.format(Msat))
        f.write('Aex = {:.2e}\n'.format(Aex))
        f.write('\n')
        f.write('// Define island geometry parameters\n')
        f.write('diam:={:.3e}\n'.format(diam*len_mult))
        f.write('l:={:.3e}\n'.format(l*len_mult))
        f.write('w:={:.3e}\n'.format(w*len_mult))
        f.write('thk:={:.3e}\n'.format(thk*len_mult))
        f.write('s1:={:.3e}\n'.format(s1*len_mult))
        f.write('s2:={:.3e}\n'.format(s2*len_mult))
        f.write('isl_rot:={:3.0f}\n'.format(np.deg2rad(isl_rot)))
        stadium_define = """
        circleleft:=circle(diam).transl(-l/2,0,0)
        circleright:=circle(diam).transl(l/2,0,0)
        stadium:=rect(l,w).add(circleleft).add(circleright)
        """
        f.write(stadium_define+'\n')
        f.write('// Build lattice\n')
        for i in range(num_isl):
            f.write('isl{0:02d}:=stadium.rotz({1:.3e}).transl({2:.3e},{3:.3e},0.0)\n'.format(i,np.deg2rad(angles[i]),centers[i,0]*len_mult-shift_cenx,centers[i,1]*len_mult-shift_ceny))
        f.write('\nlatt:=isl00')
        for i in range(1,num_isl):
            f.write('.add(isl{:02d})'.format(i))
        f.write('\n')
        f.write('setGeom(latt)\n')
        f.write('save(geom)\n')
        f.close()

        mumax_run_command = 'mumax3 '+work_dir+latt_info_name+'.mx3'
        subprocess.run(mumax_run_command.split())
        mumax_convert_command = 'mumax3-convert -png '+work_dir+latt_info_name+'.out/*.ovf'
        subprocess.run(mumax_convert_command.split())
    
    return 1

""""
End of funciton
"""