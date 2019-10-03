#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 12:46:22 2019

@author: cphatak
"""

import numpy as np
#from matplotlib import pyplot as plt
from skimage.transform import rotate as sk_rot
from comp_phase import mansPhi
from skimage import io as sk_io
 
def draw_lattice(microscope,
                 cen_fname = 'MCrun_lattice_coords_run1', #centers
                 mag_fname = 'MCrun_mag_initial_run1', #mag file name
                 dim = 400, #size of images
                 del_px = 10, #pixel resolution for images. (nm/px)
                 Lx = 300, #lenght of island along horizontal (nm)
                 Ly = 100, #length of island alog vertical (nm)
                 thk = 10, #thickness for phase computations
                 save_tfs = True, #save mag phase, ephase and TFs files
                 defocus_step = 1000000.0, #defocus for TFS images.
                 n_tfs = 3, #number of defocus images.
                 s_type = 'quadratic', #type of defocus series (linear or quadratic)
                 random_phase = False, #bckground random phase.
                 save_lattice = True, #save the lattice image with mag indication
                 buff_offset = 100, #offset for lattice to fit in the view.
                 resize_im = False): #Resizing of FOV to fit the lattice.


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
    
    #lattice offset buffer
    buff = buff_offset
    
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
    Lx /= del_px
    Ly /= del_px
    thk/= del_px
    new_cen = centers/del_px
    max_coord = np.amax(new_cen)
    
    #check the size of image.
    im_sz = dim
    if resize_im:
        if max_coord >= im_sz:
            im_sz *= 2
    
    #create single horizontal island
    dim = np.int(Lx*2)
    d2 = dim/2
    line = np.arange(dim)-d2
    X,Y = np.meshgrid(line,line)
    
    #create a single horizontal island.
    isl_img = np.zeros([dim,dim])
    
    #first the rectangular part
    isl_img[int(d2-Ly/2):int(d2+Ly/2-1),int(d2-(Lx-Ly)/2):int(d2+(Lx-Ly)/2-1)] = 1.0

    #then the circle part
    cc = np.sqrt(X**2 + Y**2)
    temp = np.zeros([dim,dim])
    temp[np.where(cc <= Ly/2)] = 1.0
    cc = np.roll(np.roll(temp,int(-(Lx-Ly)/2),axis=1),-1,axis=0)
    isl_img[np.where(cc == 1.0)] = 1.0
    cc = np.roll(np.roll(temp,int((Lx-Ly)/2),axis=1),-1,axis=0)
    isl_img[np.where(cc == 1.0)] = 1.0
    
    
    #arrays for holding variables
    latt = np.zeros([im_sz,im_sz])
    magx = np.zeros([im_sz,im_sz])
    magy = np.zeros([im_sz,im_sz])

    #offsets for lattice
    #xoff = d2-np.abs(np.amin(new_cen[:,1])) + buff
    #yoff = d2-np.abs(np.amin(new_cen[:,0])) + buff
    xoff = buff
    yoff = buff
    

    for i in range(n_isl):
        ang = np.rad2deg(np.arctan2(mag[i,1],mag[i,0]))
        temp = sk_rot(isl_img,-ang)
        xs = im_sz//2+new_cen[i,0]-xoff-d2
        xe = im_sz//2+new_cen[i,0]-xoff+d2
        ys = im_sz//2+new_cen[i,1]-yoff-d2
        ye = im_sz//2+new_cen[i,1]-yoff+d2
        #temp2 = np.roll(np.roll(temp,ys.astype('int'),axis=0),xs.astype('int'),axis=1)
        latt[int(ys):int(ye),int(xs):int(xe)] += temp
        magx[int(ys):int(ye),int(xs):int(xe)] += temp*mag[i,0]
        magy[int(ys):int(ye),int(xs):int(xe)] += temp*mag[i,1]
        #magx += temp2*mag[i,0]
        #magy += temp2*mag[i,1]
        
    if save_tfs:
        
        #compute the magnetic phase shift
        mag_phi = mansPhi(bx = magx, by = magy, thick = thk)*cb*np.pi
        #compute the ephase shift
        latt_ephi = microscope.sigma * latt_V0 * thk * latt * del_px
        #back ground phase
        if random_phase:
            mem_phi = microscope.sigma * mem_V0 * mem_thk * np.random.uniform(low = -np.pi/64, high = np.pi/64, size=latt.shape)
        else:
            mem_phi = microscope.sigma * mem_V0 * mem_thk

        #total phase 
        Tphi = mag_phi + latt_ephi + mem_phi

        #Save phase image
        sk_io.imsave(mag_fname+'_Phi.tiff',Tphi.astype('float32'))
        
        #amplitude
        Amp = np.exp((-np.ones(latt.shape) * mem_thk / mem_xip0) - (thk * del_px / latt_xip0 * latt))

        #Save amp image
        sk_io.imsave(mag_fname+'_Amp.tiff',Amp.astype('float32'))

        #Object Wave
        ObjWave = Amp * (np.cos(Tphi) + 1j * np.sin(Tphi))
        
        #coordinates for image simulations.
        #Dimensions and coordinates
        dim = im_sz
        d2=dim/2
        line = np.arange(dim)-float(d2)
        [X,Y] = np.meshgrid(line,line)
        qq = np.sqrt(X**2 + Y**2) / float(dim)
        
        #simulate images
        #im_stack = np.zeros([dim,dim,n_tfs])
        #calculate all the defocus values
        if (s_type == 'quadratic'):
        
            aa = np.arange(-(n_tfs-1)/2,(n_tfs-1)/2+1)
            def_vals = (2**(np.abs(aa)-1))*np.sign(aa) * defocus_step
        
        else:
            aa = np.arange(-(n_tfs-1)/2,(n_tfs-1)/2+1)
            def_vals = aa * defocus_step

        for nntfs in range(n_tfs):
            microscope.defocus = def_vals[nntfs]
            full_im = microscope.getImage(ObjWave,qq,del_px)
            #im_stack[:,:,nntfs] = full_im_in
            # save image
            sk_io.imsave(mag_fname+'_def_'+str(def_vals[nntfs])+'.tiff',full_im.astype('float32'))

    return 1

    
