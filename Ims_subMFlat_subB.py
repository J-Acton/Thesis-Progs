# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:16:21 2023

Description:
    a prog to normalise Mflat_subBias and then divide all images, that have had 
    the masterbias subtracted from it, by the normalised Mflat_subBias.


@author: James Acton (20325303)
"""

import ccdproc
from astropy.io import fits
from astropy.nddata import CCDData
import glob
from skimage import io, img_as_float
import numpy as np


MflatsubB33 = img_as_float(fits.getdata('F:/FYP Stuff/DBF night1/Flat/Vis/Masterflat_subMBias.FITS'))
MflatsubB11 = img_as_float(fits.getdata('F:/FYP Stuff/OHP Exoplanet files 2019/Masterflat_subMBias.FIT'))

# Defining directories for the Bias subtracted images
Tstart_Bsub_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-start/T-start_subB/'
Tmid_Bsub_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-middle/T-mid_subB/'
Tend_Bsub_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-end/T-end_subB/'

# Defining directories for the fully processed images
Tstart_final_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/Tstart/T-start_finals/start_final'
Tmid_final_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-middle/T-mid_finals/mid_final'
Tend_final_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-end/T-end_finals/end_final'

# function to normalise the masterflat frame
def flat_correct(flat):
    #Normalize the flat
    mean = np.mean(flat)
    flat_normed = flat / mean
    
    return flat_normed


NormMasterFlatsubB33 = CCDData(flat_correct(MflatsubB33), unit='adu')
NormMasterFlatsubB11 = CCDData(flat_correct(MflatsubB11), unit='adu')

#NormMasterFlatsubB.write('F:/FYP Stuff/ DBF night1/Flat/Vis/NormMasterFlatsubB.FITS)
#NormMasterFlatsubB.write('F:/FYP Stuff/OHP Exoplanet files 2019/NormMasterFlatsubB.FIT')

# function to divide bias subtracted images by the normalised masterflat frame
def Flat_div(Bsub_dir, final_dir):
    pics = sorted(glob.glob(Bsub_dir+'/*.FITS')) #creating list of names of bias subtracted images
    
    tot = len(pics)
    
    for i in range(len(pics)):
        progress = (i/tot) #progress bar stuff
        arrow = '=' * int(30 * progress) #^^^
        spaces = ' ' * (30 - len(arrow)) #^^^
        print(f'\r[{arrow}{spaces}] {int(progress * 100)}%', end='', flush=True) #^^^
        
        final = CCDData(fits.getdata(pics[i]), unit='adu').divide(NormMasterFlatsubB33) # dividing each image by normalised masterflat
        final.write(final_dir+str(i+1)+'.FITS') #saving fully processed images


        
Flat_div(Tstart_Bsub_dir, Tstart_final_dir)
Flat_div(Tmid_Bsub_dir, Tmid_final_dir)
Flat_div(Tend_Bsub_dir, Tend_final_dir)

    

    
# Same process for WASP-11 b but slightly altered as the images of this star are FIT files not FITS files
WASP11_dir = 'F:/FYP Stuff/OHP Exoplanet files 2019/Data_subB/Bsub'
WASP11_pics=[]

    
for i in range(1, 218): 
    name = glob.glob(WASP11_dir+str(i)+'.FIT')
    name = str(name)
    name = name[2:len(name)-2]
    WASP11_pics.append(str(name))
tot = len(WASP11_pics)


for i in range(len(WASP11_pics)):
    progress = (i/tot)
    arrow = '=' * int(30*progress)
    spaces = ' ' * (30 - len(arrow))
    print(f'\r[{arrow}{spaces}] {int(progress * 100)}%', end='', flush=True)
    T_start_final = CCDData(fits.getdata(WASP11_pics[i]), unit='adu').divide(NormMasterFlatsubB11)
    T_start_final.write('F:/FYP Stuff/OHP Exoplanet files 2019/Processed_Data/final'+str(i+1) + '.FIT')

