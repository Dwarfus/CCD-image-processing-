# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 10:41:40 2017

@author: Pavel
The Main code for catalogization of the galaxies
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from astropy.io import fits

#os.chdir("/Mordor/16-17/labs/astronomical")
os.chdir("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic")

class Main:
    def __init__(self):
        """
        Loads the appropriate image and gets the data as an array from it. """
        
        
        self.hdulist2 = fits.open("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic/BoundaryMasked.fits")
        self.data = self.hdulist2[0].data
        self.mean = 3418
        self.sigma = 11.87
        self.treshold = 40000 #self.mean + 5*self.sigma
        self.totalgalaxies = 0
        self.call()
#        self.x = ma.compressed(self.data)

    def getheader(self):
        """
        This one will get the header from the original file and the appropriate values from it
        """
        self.hdulist = fits.open("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic/mosaic.fits")
        self.header = self.hdulist[0].header
        self.zeropoint = self.header["MAGZPT"]
        self.erzeropoint = self.header["MAGZRR"]

    def findmax(self):
        """
        Finds the x and y coordinates of the point with highest count.
        Prints all of the coordinates, returns only the first one. 
        """
        self.posmax= np.where(self.data==np.max(self.data))
        self.count = np.max(self.data)
        #self.posmaxy = np.ndarray.argmax(self.data)
        print(self.posmax, self.count)
        return self.posmax[0][0], self.posmax[1][0], self.count

    def apperture(self):
        """
        This function calls for position of maximum, draws an apperture sphere
        around of radius 6 to measure the total magnitude. Around that sphere
        it draws a shell of radius 2 from which the background radiation is
        calculated. For the complete set of equations used see lab book. 
        """
        ymax, xmax, count = self.findmax()
        if (count> self.treshold):
            self.totalgalaxies +=1
            #a, b = 2300, 700 # These are the coordinates of the center where a is the y coordinate and b is x for y,x defined in ds9

            r = 8 # radius

            y,x = np.ogrid[-ymax:4611-ymax, -xmax:2570-xmax] # I dont really know what ogrid does. Should read later
            mask1 = x*x + y*y <= r*r # Check the circle equation, YaY!!
            self.totalcount = np.sum(self.data[mask1])
            self.totalpixels = len(self.data[mask1])
            
            r=10
            #y,x = np.ogrid[-ymax:4611-ymax, -xmax:2570-xmax] # I dont really know what ogrid does. Should read later
            mask2 = x*x + y*y <= r*r # Check the circle equation, YaY!!
            self.totalbackcount = np.sum(self.data[mask2])
            self.totalbackpixels = len(self.data[mask2])
            #The below is the local background 
            self.avbackground = (self.totalbackcount-self.totalcount)/(self.totalbackpixels-self.totalpixels)
            # This is the total count of the galaxy
            self.nobackcount = self.totalcount -self.totalpixels*self.avbackground # We need to subtract the background 
            print (self.nobackcount)
            
            self.data[mask1] = 0  
            self.apperture()
        else:
            print ("No more points above treshold")
            masked = fits.PrimaryHDU(self.data) # these two lines are all you need. 
            masked.writeto('After.fits')

    def call(self):
        self.getheader()
        self.apperture()
        
a=Main()
