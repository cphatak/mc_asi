#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 12:46:22 2019

@author: cphatak
"""

import numpy as np
from matplotlib import pyplot as plt
from skimage import draw as sk_draw
from skimage.transform import rotate as sk_rot
from comp_phase import mansPhi
from skimage import io as sk_io
from microscopes import Microscope
 
def draw_lattice(microscope,
                 cen_fname = 'MCrun_lattice_coords_run1.txt', #centers
                 mag_fname = 'MCrun_mag_initial_run1.txt', #mag file name
                 dim = 400, #size of images
                 del_px = 10, #pixel resolution for images. (nm/px)
                 Lx = 300, #lenght of island along horizontal (nm)
                 Ly = 100, #length of island alog vertical (nm)
                 thk = 10, #thickness for phase computations
                 save_tfs = True, #save mag phase, ephase and TFs files
                 save_lattice = True):#save the lattice image with mag indicatio

#dir = '/Users/cphatak/work/af_test/run1/'
#microscope = Microscope(Cs = 200.0e3, theta_c = 0.02e-3, def_spr = 80.0)
#cen_fname = dir+'MCrun_lattice_coords_run1' #centers
#mag_fname = dir+'MCrun_mag_initial_run1' #mag file name
#dim = 400 #size of images
#del_px = 10 #pixel resolution for images. (nm/px)
#Lx = 300 #lenght of island along horizontal (nm)
#Ly = 100 #length of island alog vertical (nm)
#thk = 10 #thickness for phase computations
#save_tfs = True #save mag phase, ephase and TFs files
#save_lattice = True#save the lattice image with mag indicatio    

    #some pre-dfined parameters
    #microscope defocus
    def_val = 1000000.0

    #offset values for shifing the islands of lattice
    offx = 1000/del_px#im_sz/2#/del_px
    offy = 800/del_px
    
    #material parameters for phase and amp. computation.
    #lattice values
    latt_V0 = 20.0 #V
    latt_xip0 = 200 #nm
    #Magnetization parameters
    b0 = 1.0e4 #Gauss
    phi0 = 20.7e6 #Gauss.nm^2
    cb = b0/phi0*del_px**2 #1/px^2
    
    #Support membrane
    mem_thk = 50.0 #nm
    mem_V0 = 10.0 #V
    mem_xip0 = 800.0 #nm
    
    # read the data
    centers = np.genfromtxt(cen_fname+'.txt',delimiter=',',skip_header=2)
    mag = np.genfromtxt(mag_fname+'.txt', delimiter=',', skip_header=1)
    
    #get number of islands
    n_isl, coord = centers.shape
    
    #convert to pixel coordinates
    new_cens = centers/del_px
    max_coord = np.amax(new_cens)
    
    im_sz = dim
    if max_coord >= im_sz:
        im_sz *= 2
    
    isl_img = np.zeros([im_sz,im_sz])
    
    #create a single horizontal island.
    Lx /= del_px
    Ly /= del_px
    thk/= del_px
    
    #create the rectangle
    r_start = (im_sz/2-(Ly)/2,im_sz/2-(Lx-Ly)/2)
    rr,cc = sk_draw.rectangle(r_start, extent = (Ly,(Lx-Ly)))
    #rr,cc = sk_draw.rectangle(r_start, extent = (Ly,Lx))
    isl_img[rr.astype('int'),cc.astype('int')] = 1
    
    #add the circles
    rr,cc = sk_draw.circle(im_sz/2,im_sz/2-(Lx-Ly)/2,Ly/2)
    isl_img[rr.astype('int'),cc.astype('int')] = 1
    rr,cc = sk_draw.circle(im_sz/2,im_sz/2+(Lx-Ly)/2,Ly/2)
    isl_img[rr.astype('int'),cc.astype('int')] = 1
    
    
    #arrays for holding variables
    latt = np.zeros([im_sz,im_sz])
    magx = np.zeros([im_sz,im_sz])
    magy = np.zeros([im_sz,im_sz])
    
    for i in range(n_isl):
        ang = np.rad2deg(np.arctan2(mag[i,1],mag[i,0]))
        temp = sk_rot(isl_img,-ang)
        xs = new_cens[i,0]-offx
        ys = new_cens[i,1]-offy
        temp2 = np.roll(np.roll(temp,ys.astype('int'),axis=0),xs.astype('int'),axis=1)
        magx += temp2*mag[i,0]
        magy += temp2*mag[i,1]
        latt += temp2
        
    if save_tfs:
        
        #compute the magnetic phase shift
        mag_phi = mansPhi(bx = magx, by = magy, thick = thk)*cb*np.pi
        #compute the ephase shift
        latt_ephi = microscope.sigma * latt_V0 * thk * latt * del_px
        #back ground phase
        mem_phi = microscope.sigma * mem_V0 * mem_thk * np.random.uniform(
                low = -np.pi/32, high = np.pi/32, size=latt.shape)
        #total phase 
        Tphi = mag_phi + latt_ephi #+ mem_phi
        
        #amplitude
        Amp = np.exp((-np.ones(latt.shape) * mem_thk / mem_xip0) - (thk / latt_xip0 * latt))
        #Object Wave
        ObjWave = Amp * (np.cos(Tphi) + 1j * np.sin(Tphi))
        
        #coordinates for image simulations.
        #Dimensions and coordinates
        d2=dim/2
        line = np.arange(dim)-float(d2)
        [X,Y] = np.meshgrid(line,line)
        th = np.arctan2(Y,X)
        qq = np.sqrt(X**2 + Y**2) / float(dim)
        
        #simulate images
        microscope.defocus = 0.0
        full_im_in = microscope.getImage(ObjWave,qq,del_px)
        microscope.defocus = -def_val
        full_im_un = microscope.getImage(ObjWave,qq,del_px)
        microscope.defocus = def_val
        full_im_ov = microscope.getImage(ObjWave,qq,del_px)
        
        
        sk_io.imsave(mag_fname+'_Phi.tiff',Tphi.astype('float32'))
        sk_io.imsave(mag_fname+'_Im_In.tiff',full_im_in.astype('float32'))
        sk_io.imsave(mag_fname+'_Im_Un.tiff',full_im_un.astype('float32'))
        sk_io.imsave(mag_fname+'_Im_Ov.tiff',full_im_ov.astype('float32'))
    
