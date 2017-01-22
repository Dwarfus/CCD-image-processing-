import os
import numpy as np
import numpy.ma as ma
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

os.chdir("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic")

from astropy.io import fits
hdulist = fits.open("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic/mosaic.fits")
header = hdulist[0].header

data = hdulist[0].data

a, b = 3100, 1400 # These are the coordinates of the center where a is the y coordinate and b is x for y,x defined in ds9

r = 290 # radius, obviously,DOH :D

y,x = np.ogrid[-a:4611-a, -b:2570-b] # I dont really know what ogrid does. Should read later
mask = x*x + y*y <= r*r # Check the circle equation, YaY!!

data[mask] = 0


masked = fits.PrimaryHDU(data) # these two lines are all you need. 
masked.writeto('masked9.fits')
