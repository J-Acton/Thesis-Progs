# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 19:30:04 2023

@author: James Acton
Description:
    A prog to combine all flats taken with the visible light filter on the 1.2m
    telescope OHP into a masterflat file and then subtract the MasterBias file
    from the masterflat file
"""

import ccdproc
from astropy.io import fits
from astropy.nddata import CCDData
import glob


# Creating a masterflat file

# Defining directories
flat_dir33 = 'F:/FYP Stuff/DBF night1/Flat/Vis/'
flat_dir11 = 'F:/FYP Stuff/OHP Exoplanet files 2019/Flats/'

# creating lists of flat file names
flatNames33 = sorted(glob.glob(flat_dir33+'/*.FITS'))
flatNames11 = sorted(glob.glob(flat_dir11+'/*.FIT'))

# creating empty 
flats33 = []
flats11 = []

# Retrieving FITS files from list of file names
for i in flatNames33:
    flats33.append(CCDData(fits.getdata(i), unit = 'adu'))

for i in flatNames11:
    flats11.append(CCDData(fits.getdata(i), unit = 'adu'))
  
# combining flat frames into masterflat
Masterflat33 = ccdproc.combine(flats33)
Masterflat11 = ccdproc.combine(flats11)

# retrieving masterbias files
Mbias33 = CCDData(fits.getdata('F:/FYP Stuff/ DBF night1/Bias/MasterBias/Masterbias.FITS'), unit='adu')
Mbias11 = CCDData(fits.getdata('F:/FYP Stuff/OHP Exoplanet files 2019/Biases/Masterbias.FIT'), unit='adu')

# Subtracting masterbias from masterflat
Masterflat_subMBias33 = Masterflat33.subtract(Mbias33)
Masterflat_subMBias33.write('F:/FYP Stuff/DBF night1/Flat/Vis/Masterflat_subMBias.FITS')

Masterflat_subMBias11 = Masterflat11.subtract(Mbias11)
Masterflat_subMBias11.write('F:/FYP Stuff/OHP Exoplanet files 2019/Masterflat_subMBias.FIT')

