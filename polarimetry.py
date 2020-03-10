# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 12:06:49 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as unit
from photutils import DAOStarFinder, aperture_photometry, CircularAperture
from astropy.stats import mad_std
#import aplpy


fits_dir = "C:\\Users\\premv\\Documents\\acad\\gradschool\\term-3\\Astrotech-I\\Experiments\\Exp-1\\Data\\"
#fits_subdir = "Polarized\\"
fits_subdir = "Unpolarized\\"


gain = 2.5
fwhm = 26.
r = 20

fits_filenames = {}
R = {}
No = {}
Ne = {}
Io = {}
Ie = {}
bkgnd = {}
sig_R = {}
K = {}
q = {}
u = {}
p = {}
theta = {}
sig_q = {}
sig_u = {}
sig_p = {}
sig_theta = {}



for sn in range(1,6):
    fits_filenames[sn] = {}
    R[sn] = {}
    No[sn] = {}
    Ne[sn] = {}
    Io[sn] = {}
    Ie[sn] = {}
    bkgnd[sn] = {}
    
    sig_R[sn] = {}
    
    
    for alpha in range(0,900,900//4):        
        fits_filenames[sn][alpha] = "Set{0}_{1}.FIT".format(sn, alpha)
        fits_filepath = fits_dir + fits_subdir + fits_filenames[sn][alpha]
        #gc = aplpy.FITSFigure(fits_file)
        #gc.show_grayscale()
        
        hdul = fits.open(fits_filepath)
        image = np.float64(hdul[0].data)
        
        x_pixels = hdul[0].header['NAXIS1'] * unit.pix
        T = hdul[0].header['EXPTIME']
        hdul.close()

        bkgnd[sn][alpha] = np.median(image)
        image -= bkgnd[sn][alpha]
        bkg_sigma = mad_std(image)
        daofind = DAOStarFinder(fwhm=fwhm, threshold=5.*bkg_sigma)  
        sources = daofind(image)  
#        for col in sources.colnames:  
#            sources[col].info.format = '%.8g'  # for consistent table output
#        print(sources)  
        
        positions = np.transpose((sources['xcentroid'],sources['ycentroid']))
        
        apertures = CircularAperture(positions,r=r)
        
        plt.imshow(image)
        apertures.plot()
        
        photo_apertures = aperture_photometry(image, apertures)
#        print(photo_apertures)
        
        Ne[sn][alpha] = np.max(photo_apertures['aperture_sum'][photo_apertures['xcenter']<x_pixels/2]) / gain
        No[sn][alpha] = np.max(photo_apertures['aperture_sum'][photo_apertures['xcenter']>x_pixels/2]) / gain
        
        Ie[sn][alpha] = Ne[sn][alpha] / T
        Io[sn][alpha] = No[sn][alpha] / T
        
    
    K[sn] = (Io[sn][0]*Io[sn][225]*Io[sn][450]*Io[sn][675])**.25 / (Ie[sn][0]*Ie[sn][225]*Ie[sn][450]*Ie[sn][675])**.25
    
    for alpha in range(0,900,900//4):
        Ie[sn][alpha] *= K[sn]
        R[sn][alpha] = (Io[sn][alpha]-Ie[sn][alpha])/(Io[sn][alpha]+Ie[sn][alpha])
        
        No[sn][alpha] += bkgnd[sn][alpha] * np.pi * r**2 / gain
        Ne[sn][alpha] += bkgnd[sn][alpha] * np.pi * r**2 / gain
        sig_R[sn][alpha] = ( (4* No[sn][alpha] * Ne[sn][alpha]) / (No[sn][alpha] + Ne[sn][alpha])**3 )**.5
        
    q[sn] = R[sn][0]    
    u[sn] = R[sn][225]
    
    sig_q[sn] = sig_R[sn][0]    
    sig_u[sn] = sig_R[sn][225]
    
    p[sn] = (q[sn]**2 + u[sn]**2)**.5
    theta[sn] = np.arctan2(u[sn],q[sn])/2 * 180 / np.pi
    if theta[sn] < 0:
        theta[sn] += 180
    
    sig_p[sn] = ((q[sn]**2* sig_q[sn]**2 + u[sn]**2 * sig_u[sn]**2) / (q[sn]**2 + u[sn]**2))**.5
    sig_theta[sn] = ( (q[sn]**2* sig_u[sn]**2 + u[sn]**2 * sig_q[sn]**2) / (q[sn]**2 + u[sn]**2)**2 )**.5 / 2 * 180 / np.pi
    
    
    print ("Degree of polarisation in set {0}, p = {1:.3} +/- {2:.3}".format(sn,p[sn],sig_p[sn]) )
    print ("Angle of polarisation in set {0}, theta = {1:.4} +/- {2:.3}".format(sn,theta[sn],sig_theta[sn]) )




p_bar = 0
varinv_p_bar = 0

theta_bar = 0
varinv_theta_bar = 0

for sn in range(1,6):
    p_bar += p[sn] / sig_p[sn]**2
    varinv_p_bar += 1 / sig_p[sn]**2
    
    theta_bar += theta[sn] / sig_theta[sn]**2
    varinv_theta_bar += 1 / sig_theta[sn]**2

p_bar /= varinv_p_bar
theta_bar /= varinv_theta_bar

sig_p_bar = varinv_p_bar**-.5
sig_theta_bar = varinv_theta_bar**-.5


print("Degree of polarisation, p = {1:.3} +/- {2:.3}".format(sn,p_bar,sig_p_bar) )

print("Angle of polarisation, theta = {1:.4} +/- {2:.3} deg".format(sn,theta_bar,sig_theta_bar) )










