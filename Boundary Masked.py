import os
import numpy as np
import numpy.ma as ma


os.chdir("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic")

from astropy.io import fits
hdulist = fits.open("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic/mosaic.fits")
header = hdulist[0].header

data = hdulist[0].data
# If you want to crop the image, use the one below
#c=data[119:4513,115:2470]

# If you want to just set the boundaries to zero, use this one below
data[:119,:]=0
data[4513:,:] = 0
data[:,:115] =0
data[:,2470:] = 0
masked = fits.PrimaryHDU(data) # these two lines are all you need. 
masked.writeto('BoundaryMasked.fits')
