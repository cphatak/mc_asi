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
    
    #some pre-dfined parameters
    
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
    cb = b0/phi0#*del_px**2 #1/px^2

    #Support membrane
    mem_thk = 50.0 #nm
    mem_V0 = 10.0 #V
    mem_xip0 = 800.0 #nm
    
    #get number of islands
    n_isl, coord = centers.shape
    
    #convert to pixel coordinates
    new_cens = centers/del_px
    max_coord = np.amax(new_cens)
    
    im_sz = dim
    if max_coord >= im_sz:
        im_sz = 512
    
    isl_img = np.zeros([im_sz,im_sz])
    
    #create a single horizontal island.
    Lx /= del_px
    Ly /= del_px
    
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
        mag_phi = mansPhi(bx = magx, by = magy, thick = thk)*cb
        #compute the ephase shift
        latt_ephi = microscope.sigma * latt_V0 * thk
        #back ground phase
        mem_phi = microscope.sigma * mem_V0 * mem_thk * np.random.uniform(
                low = -np.pi, high = np.pi, size=disc.shape)
        #total phase 
        Tphi = mag_phi + latt_ephi + mem_phi
        
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
        
    
    
    
    plt.imshow(latt)
    plt.show()
