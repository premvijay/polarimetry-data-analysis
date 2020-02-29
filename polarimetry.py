# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 12:06:49 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
#import photutils.aperture
from photutils import DAOStarFinder, aperture_photometry, CircularAperture
from astropy.stats import mad_std
#import aplpy

fits_file = r"C:\Users\premv\Documents\acad\gradschool\term-3\Astrotech-I\Experiments\Exp-1\Data\Polarized\Set1_0.FIT"

#gc = aplpy.FITSFigure(fits_file)

#gc.show_grayscale()

hdul = fits.open(fits_file)
image = np.float64(hdul[0].data)
x_pixels = image.shape[0] * u.pix
fwhm = 16.

image -= np.median(image)
bkg_sigma = mad_std(image)
daofind = DAOStarFinder(fwhm=fwhm, threshold=10.*bkg_sigma)  
sources = daofind(image)  
for col in sources.colnames:  
    sources[col].info.format = '%.8g'  # for consistent table output
print(sources)  

positions = np.transpose((sources['xcentroid'],sources['ycentroid']))

apertures = CircularAperture(positions,r=fwhm)



plt.imshow(image)
apertures.plot()

photo_apertures = aperture_photometry(image, apertures)

print(photo_apertures)





N0 = np.max(photo_apertures['aperture_sum'][photo_apertures['xcenter']<x_pixels/2])
N1 = np.max(photo_apertures['aperture_sum'][photo_apertures['xcenter']>x_pixels/2])

R = (N0-N1)/(N0+N1)







#photo_aperture_array = np.array((photo_apertures['xcenter'],photo_apertures['ycenter'],photo_apertures['aperture_sum'])).T
#
#for photo_aperture in photo_apertures:
#    print(photo_aperture)
#    photo_aperture







