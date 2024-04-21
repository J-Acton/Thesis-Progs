# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 22:10:27 2023

@author: James Acton
Description:
    a prog to subtract the MasterBias.Fits from all images taken of WASP-33
    using the 1.2m at OHP
"""

import ccdproc
from astropy.io import fits
from astropy.nddata import CCDData
import glob


# Retrieving the MasterBias frames for WASP-33 and WASP-11
Mbias33 = CCDData(fits.getdata('F:/FYP Stuff/DBF night1/Bias/MasterBias/Masterbias.FITS'), unit='adu')
Mbias11 = CCDData(fits.getdata('F:/FYP Stuff/OHP Exoplanet files 2019/Biases/Masterbias.FIT'), unit='adu')

# Defining directories for each set of data
Tstart_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-start/'
Tmid_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-middle/'
Tend_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-end/'

WASP11_dir = 'F:/FYP Stuff/OHP Exoplanet files 2019/Raw_Data/'

def Mbias_sub(direc, Mbias, folder):
    Picnames = (glob.glob(direc+'/*.FITS')) # creating list of FITS file to be bias subtracted
    tot = len(Picnames)
    
    
    for i in range(len(Picnames)):
        progress = (i/tot)
        arrow = '=' * int(30 * progress) #progress bar for the loop
        spaces = ' ' * (30 - len(arrow)) #^^^
        print(f'\r[{arrow}{spaces}] {int(progress * 100)}%', end='', flush=True) #^^^
        subB = CCDData(fits.getdata(Picnames[i]), unit='adu').subtract(Mbias) #subtracting Masterbias from FITS files
        subB.write(direc+folder+str(i+1) + '.FITS') # saving FITS files to folder
        

# Subtracting WASP-33's masterbias from all of the raw images of WASP-33
Mbias_sub(Tstart_dir, Mbias33, 'T-start-subB/Bsub')
Mbias_sub(Tmid_dir, Mbias33, 'T-mid_subB/mid_Bsub')
Mbias_sub(Tend_dir, Mbias33, 'T-end_subB/end_Bsub')


# Subtracting WASP-11's masterbias from all of the raw images of WASP-11, these were FIT files not FITS file so couldn't use Mbias_sub function
WASP11_Picnames = (glob.glob(WASP11_dir+'/*.FIT')) 
tot = len(WASP11_Picnames)


for i in range(len(WASP11_Picnames)):
    progress = (i/tot)
    arrow = '=' * int(30 * progress) #progress bar for the loop
    spaces = ' ' * (30 - len(arrow)) #^^^
    print(f'\r[{arrow}{spaces}] {int(progress * 100)}%', end='', flush=True) #^^^
    WASP11_subB = CCDData(fits.getdata(WASP11_Picnames[i]), unit='adu').subtract(Mbias11)
    WASP11_subB.write('F:/FYP Stuff/OHP Exoplanet files 2019/Data_subB/Bsub'+str(i+1) + '.FIT')



#-----------------------------------------------------------------------------#
