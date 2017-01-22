import os
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.optimize import curve_fit
#os.chdir("/Mordor/16-17/labs/astronomical")
os.chdir("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic")

from astropy.io import fits

class Gaussian:
    def function(x, x0, a, sigma):
        return a*np.e**(-(x-x0)**2/(2*sigma**2))


hdulist = fits.open("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic/BoundaryMasked.fits")
header = hdulist[0].header

data = hdulist[0].data

x = ma.compressed(data)

hcounts, intens = np.histogram(x, bins = 1000, range=(3000,4000))
peakint = np.argmax(hcounts)
background = intens[peakint]

jnk2 = intens[:-1]
popt, pcov = curve_fit(Gaussian.function, jnk2,hcounts, p0=[background, 5000,20] )

plt.figure(1)
plt.title("The histogram of counts")
plt.ylabel("Number of occurances")
plt.xlabel("Number of counts")

ordata, = plt.plot(jnk2,hcounts, label='Original data')
gaus, = plt.plot(jnk2, Gaussian.function(jnk2, *popt), label = 'Fitted Gaussian')

plt.legend(handles=[ordata, gaus])
plt.show

print ("The peak of the histogram is at:",background)
print("THe fitted gaussian mean, norm.constant and stand dev is:", popt)



