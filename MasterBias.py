# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 00:25:36 2023

@author: James Acton
Description:
    prog to create a masterbias image of the bias images taken for the 1.2m
    telescope at OHP.
"""


import ccdproc
from astropy.io import fits
from astropy.nddata import CCDData
import glob

# Defining the directories of the bias frames for WASP-33 and WASP-11
bias_dir_33 = 'F:/FYP Stuff/DBF night1/Bias/'
bias_dir_11 = 'F:/FYP Stuff/OHP Exoplanet files 2019/Biases/'

# Creating a list of the bias file names for each star
biasNames_33 = sorted(glob.glob(bias_dir_33+'/*FITS'))
biasNames_11 = sorted(glob.glob(bias_dir_11+'/*.FIT'))

# Creating empty lists to be filled with the bias frames
biases_33 = []
biases_11 = []

# Retrieving and reading in the bias frames
for i in biasNames_33:
    biases_33.append(CCDData(fits.getdata(i), unit = 'adu'))
    
for i in biasNames_11:
    biases_11.append(CCDData(fits.getdata(i), unit = 'adu'))
    
# Combine the bias frames for each star into a masterbias frame
Masterbias_33 = ccdproc.combine(biases_33, bias_dir_33+'Masterbias.FITS', 'median')
Masterbias_11 = ccdproc.combine(biases_11, bias_dir_11+'Masterbias.FIT', 'median')