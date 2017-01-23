# -*- coding: utf-8 -*-
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
        Loads the appropriate image and gets the data as an array from it.
        Initiate all the necessary variables"""
        
        
        self.hdulist2 = fits.open("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic/img_test.fits")
        self.data = self.hdulist2[0].data # this is the data from fits file
        self.mean = 3418 # mean calculated from histogram
        self.sigma = 11.87 # Standard deviation from histogram
        self.treshold =  self.mean + 5*self.sigma # THe treshold with user picked number of sigma
        self.totalgalaxies = 0 # This is the total count of galaxies, I guess we can just use len(catalog) but whatever
        self.overlap = 0 # This will count the number of occurences where we discard a point as a residual, e.g. the apperture average is lower than treshold
        self.catalog = [] # The catalog values, y, x coordinates of the centre, magnitude
        self.call()

    def getheader(self):
        """
        This one will get the header from the original file and the appropriate values from it.
        Technically can be changed for manual input of the values, but it might be useful
        to have the header anyway. The header cant be loaded from the masked file
        as there the header doesn't contain the info
        """
        self.hdulist = fits.open("/Users/Pavel/Documents/aaaaaaaa/Anglie/imperial/2016-2017/labs/astronomy/mosaic/mosaic.fits")
        self.header = self.hdulist[0].header
        self.zeropoint = self.header["MAGZPT"]
        self.erzeropoint = self.header["MAGZRR"]

    def findmax(self):
        """
        Finds the x and y coordinates of the point with highest count.
        Returns the first pair of coordinates, so only one point is considered
        """
        self.posmax= np.where(self.data==np.max(self.data))
        self.count = np.max(self.data)
        
        return self.posmax[0][0], self.posmax[1][0], self.count

    def apperture(self):
        """
        This function calls for position of maximum.
        It checks that the count value is above treshold. If yes, proceed with
        algorithm, if not it exports the file.
        For values above treshold, draws two circles, if average count of the
        inner circle above treshold, it proceeds (this is a check for counting
        the residual from other masked galaxies as a new galaxy). If not, it 
        masks the point and add one to residual count.
        If proceed it calculates the contribution from non zero points,
        calculates the background from non zero outer shell
        gets the total count and calculates the magnitude of brightness of whatever
        saves y,x, mag into catalog, masks the inner radius sphere
        """
        ymax, xmax, count = self.findmax() #Gets the position of the highest value
        #Treshold condition
        if (count> self.treshold):
                                   
            r = 8 # radius of inner sphere

            y,x = np.ogrid[-ymax:4611-ymax, -xmax:2570-xmax] # I dont really know what ogrid does. Should read later
            mask1 = x*x + y*y <= r*r # Do a circle mask of given radius            
            # Check that the average of the mask is above treshold... e.g. it is not a residual
            if (np.mean(self.data[mask1])>self.treshold):

                b = ma.masked_where(self.data[mask1]!=0,self.data[mask1]) # find all zeros values
                circapp=self.data[mask1] # Creates an array with only the values within a circle
                nonzeroapp = circapp[b.mask]  # Creates an array with non zero values within the circle
                self.totalcount = np.sum(nonzeroapp) # Total counts from the circle
                self.totalpixels = len(nonzeroapp) # Total number of non zero pixels in the circle
                
                r=10 # The radius of the outer circle
                

                mask2 = x*x + y*y <= r*r # Second circle mask
                b2 = ma.masked_where(self.data[mask2]!=0, self.data[mask2]) # Again check for zero values
                circapp2 = self.data[mask2] # creates an array of just values in outer circle
                nonzeroapp2 = circapp2[b2.mask] # Creates an array of non zero values in outer circle
                self.totalbackcount = np.sum(nonzeroapp2) # Total count
                self.totalbackpixels = len(nonzeroapp2) # TOtal number of non zero pixels
                self.avbackground = (self.totalbackcount-self.totalcount)/(self.totalbackpixels-self.totalpixels) # Background per pixel
                
                self.nobackcount = self.totalcount -self.totalpixels*self.avbackground # Total number of counts without the background
                
                
                mag = self.zeropoint -2.5*np.log10(self.nobackcount) # calculating the calibrated magnitude 
            
                self.catalog.append([ymax,xmax,mag]) # Adding to the catalog
                self.totalgalaxies +=1 # non necessary, can take the len of catalog
                self.data[mask1] = 0 #mask the data that was just evaluated.

            # If the point is a residual, it will mask it and add one to overlap count so we can optimize the radius choice
            else: 
                self.data[ymax,xmax]=0
                self.overlap +=1
            
            #self.apperture()
        # When there are no more points above the treshold, it will terminate and export the image. 
        else:
            print ("No more points above treshold")
            #masked = fits.PrimaryHDU(self.data) # these two lines are all you need. 
            #masked.writeto('After.fits')

            

    def call(self):
        """
        Just a call function that calls the appropriate method.
        I like to have the order of methods at one place instead of calling one method fromother"""
        self.getheader()
        i=0
        # THis loop will not be here. It will be looped by having self.apperture() in the apperture function to run as long as points above treshold.
        # For now, we want to run it just some number of times so this is done.
        while i<200:
            self.apperture()
            i+=1
            
        masked = fits.PrimaryHDU(self.data) # these two lines are all you need. 
        masked.writeto('After4.fits')
        
a=Main()
