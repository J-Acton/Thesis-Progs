# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 15:28:49 2023

Description:


@author: James Acton (20325303)
"""

from photutils.aperture import CircularAperture as cap
from photutils.aperture import aperture_photometry as ap
from astropy.io import fits
import glob
import numpy as np


#Defining the directories for the fully processed fits files
Tstart_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-start/T-start_finals/'
Tmid_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-middle/T-mid_finals/'
Tend_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-end/T-end_finals/'

WASP11_dir = 'F:/FYP Stuff/OHP Exoplanet files 2019/Processed_Data/'

#Making list of T-start, middle and end fits files
startpics = (glob.glob(Tstart_dir+'/*.FITS'))
midpics = (glob.glob(Tmid_dir+'/*.FITS'))
endpics = (glob.glob(Tend_dir+'/*.FITS'))

#Making a list of pre-processed start, middle, and end pics
prepro_sp = (glob.glob('F:/FYP Stuff/SCJA WASP33b/N1/T-start/*.FITS'))
prepro_mp = (glob.glob('F:/FYP Stuff/SCJA WASP33b/N1/T-middle/*.FITS'))
prepro_ep = (glob.glob('F:/FYP Stuff/SCJA WASP33b/N1/T-end/*.FITS'))

WASP11_pics = (glob.glob(WASP11_dir+'/*FIT'))
#-----------------------------------------------------------------------------#
r = 25 #8 #25 #radius of circular aperture in pixels, r = 8 for WASP33, r = 25 for WASP11

# function for retrieving the mean background pixel values from a text file
def mean_retriev(file):
    file = open('F:/FYP Stuff/Photometry data/'+file, 'r')
    bkg_mean = file.readlines()
    file.close()
    
    means = []
    
    for i in range(len(bkg_mean)):
        means.append(float(bkg_mean[i].split('\n')[0]))
        
    means=np.array(means)
    
    return means

#-----------------------------------------------------------------------------#

# function for retrieving the centre coords of a star from a text file
def centre_retriev(file):
    coords = open('F:/FYP Stuff/Photometry data/'+file, 'r')
    
    x , y, vals = [], [], []
    for line in coords:
        vals.append(line.split(' centre at: '))
        
    for n in vals:
        n[-1] = n[-1].split(' ')
        
    coords.close()
    
    for i in range(len(vals)):
        x.append(int(vals[i][1][0]))
        y.append(int(vals[i][1][1]))
        
    x = np.array(x)
    y = np.array(y)
    
    return x, y
    
#-----------------------------------------------------------------------------#

# function to generate a circular aperature around a star given the centre coordinates
def Flux_finder(pics, file, x, y, means):
    data = []
    for i in range(len(pics)):
        data.append(fits.getdata(pics[i]))
        
    flux_txt = open('F:/FYP Stuff/Photometry data/'+file, 'w')
    pos = (np.array([x, y])).T
    
    for i in range(len(data)):
        aper = cap(pos[i], r) #draws circular aperture on each fits file
        bkg_sub = data[i] - means[i] #subtracts background count from fits file
        flux_val = ap(bkg_sub, aper) #finds flux value from bkg subbed fits files
        flux_txt.write('\n'+str(flux_val[0]['aperture_sum']))
    
    flux_txt.close()
    return aper
    
        
#-----------------------------------------------------------------------------#

#---------------------------retrieving means----------------------------------#

#retrieving means for T-Start
st_means = mean_retriev('T-Start/st_bkg_mean.txt')
BD_st_means = mean_retriev('T-Start/BDst_bkg_mean.txt')
TYC_st_means = mean_retriev('T-Start/TYCst_bkg_mean.txt')
TYC2_st_means = mean_retriev('T-Start/TYC2st_bkg_mean.txt')
TYC3_st_means = mean_retriev('T-Start/TYC3st_bkg_mean.txt')

#retrieving means for T-Middle
m_means = mean_retriev('T-Middle/m_bkg_mean.txt')
BD_m_means = mean_retriev('T-Middle/BDm_bkg_mean.txt')
TYC_m_means = mean_retriev('T-Middle/TYCm_bkg_mean.txt')
TYC2_m_means = mean_retriev('T-Middle/TYC2m_bkg_mean.txt')
TYC3_m_means = mean_retriev('T-Middle/TYC3m_bkg_mean.txt')

#retrieving means for T-End
e_means = mean_retriev('T-End/e_bkg_mean.txt')
BD_e_means = mean_retriev('T-End/BDe_bkg_mean.txt')
TYC_e_means = mean_retriev('T-End/TYCe_bkg_mean.txt')
TYC2_e_means = mean_retriev('T-End/TYC2e_bkg_mean.txt')
TYC3_e_means = mean_retriev('T-End/TYC3e_bkg_mean.txt')

WASP11_means = mean_retriev('WASP11/bkg_mean.txt')
LSPM_means = mean_retriev('WASP11/LSPM_bkg_mean.txt')


#retrieving means for pre-processed images
prepro_sbm = mean_retriev('T-Start/prepro_sbm.txt')
prepro_mbm = mean_retriev('T-Middle/prepro_mbm.txt')
prepro_ebm = mean_retriev('T-End/prepro_ebm.txt')

#-----------------------------------------------------------------------------#

#-----------------------retrieving centre coords------------------------------#

#retrieving centre coords for T-start
xs, ys = centre_retriev('T-Start/Tstart_maxs.txt')
BD_xs, BD_ys = centre_retriev('T-Start/stbkg_maxs.txt')
TYC_xs, TYC_ys = centre_retriev('T-Start/stbkg_maxs2.txt')
TYC2_xs, TYC2_ys = centre_retriev('T-Start/stbkg_maxs3.txt')
TYC3_xs, TYC3_ys = centre_retriev('T-Start/stbkg_maxs4.txt')

#retrieving centre coords for T-middle
xm, ym = centre_retriev('T-Middle/Tmid_maxs.txt')
BD_xm, BD_ym = centre_retriev('T-Middle/mbkg_maxs.txt')
TYC_xm, TYC_ym = centre_retriev('T-Middle/mbkg_maxs2.txt')
TYC2_xm, TYC2_ym = centre_retriev('T-Middle/mbkg_maxs3.txt')
TYC3_xm, TYC3_ym = centre_retriev('T-Middle/mbkg_maxs4.txt')

#retrieving centre coords for T-end
xe, ye = centre_retriev('T-End/Tend_maxs.txt')
BD_xe, BD_ye = centre_retriev('T-End/ebkg_maxs.txt')
TYC_xe, TYC_ye = centre_retriev('T-End/ebkg_maxs2.txt')
TYC2_xe, TYC2_ye = centre_retriev('T-End/ebkg_maxs3.txt')
TYC3_xe, TYC3_ye = centre_retriev('T-End/ebkg_maxs4.txt')

#retrieving centre coords for WASP-11 and reference star
WASP11_x, WASP11_y = centre_retriev('WASP11/maxs.txt')
LSPM_x, LSPM_y = centre_retriev('WASP11/LSPM_maxs.txt')
#-----------------------------------------------------------------------------#

#--------------------Aperture Photometry (finally!)---------------------------#
#code commented out so as to not overwrite existing text files
"""
#T-Start
Flux_finder(startpics, 'T-Start/st_flux.txt', xs, ys, st_means)
Flux_finder(startpics, 'T-Start/BD_st_flux.txt', BD_xs, BD_ys, BD_st_means)
Flux_finder(startpics, 'T-Start/TYC_st_flux.txt', TYC_xs, TYC_ys, TYC_st_means)
Flux_finder(startpics, 'T-Start/TYC2_st_flux.txt', TYC2_xs, TYC2_ys, TYC2_st_means)
Flux_finder(startpics, 'T-Start/TYC3_st_flux.txt', TYC3_xs, TYC3_ys, TYC3_st_means)


#T-Middle
Flux_finder(midpics, 'T-Middle/mid_flux.txt', xm, ym, m_means)
Flux_finder(midpics, 'T-Middle/BD_mid_flux.txt', BD_xm, BD_ym, BD_m_means)
Flux_finder(midpics, 'T-Middle/TYC_mid_flux.txt', TYC_xm, TYC_ym, TYC_m_means)
Flux_finder(midpics, 'T-Middle/TYC2_mid_flux.txt', TYC2_xm, TYC2_ym, TYC2_m_means)
Flux_finder(midpics, 'T-Middle/TYC3_mid_flux.txt', TYC3_xm, TYC3_ym, TYC3_m_means)

#T-End
Flux_finder(endpics, 'T-End/end_flux.txt', xe, ye, e_means)
Flux_finder(endpics, 'T-End/BD_end_flux.txt', BD_xe, BD_ye, BD_e_means)
Flux_finder(endpics, 'T-End/TYC_end_flux.txt', TYC_xe, TYC_ye, TYC_e_means)
Flux_finder(endpics, 'T-End/TYC2_end_flux.txt', TYC2_xe, TYC2_ye, TYC2_e_means)
Flux_finder(endpics, 'T-End/TYC3_end_flux.txt', TYC3_xe, TYC3_ye, TYC3_e_means)

#Pre-processed images
Flux_finder(prepro_sp, 'T-Start/prepro_sf.txt', xs, ys, prepro_sbm)
Flux_finder(prepro_mp, 'T-Middle/prepro_mf.txt', xm, ym, prepro_mbm)
Flux_finder(prepro_ep, 'T-End/prepro_ef.txt', xe, ye, prepro_ebm)
"""


"""
#WASP11
Flux_finder(WASP11_pics, 'WASP11/flux.txt', WASP11_x, WASP11_y, WASP11_means)
 
#WASP-11 Reference Star - LSPM J0309+3042
Flux_finder(WASP11_pics, 'WASP11/LSPM_flux.txt', LSPM_x, LSPM_y, LSPM_means)
"""

# small bit of code used to measure the number of pixels within the aperture generated
"""
st_pic1 = fits.getdata(startpics[0])
# finding number of pixels in aperture of radius 8  
test_aper = Flux_finder([startpics[0]], 'T-Start/aper_test.txt', xs, ys, st_means)
#test_aper = Flux_finder([WASP11_pics[0]], 'T-Start/aper_test.txt', WASP11_x, WASP11_y, WASP11_means)
test_mask = test_aper.to_mask()

tot = 0
tot_tru = 0

for i in test_mask.data:
    for j in i:
        tot += j
print(tot)
print(tot_tru)
"""