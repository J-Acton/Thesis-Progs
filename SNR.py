# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:51:30 2024

Description:
    A prog to calculate the signal-to-noise ratio of each image, as well as
    the average SNR.

@author: James Acton (20325303)
"""

import numpy as np
from astropy.io import fits
import glob


#Defining the directories for the fully processed fits files
Tstart_dir = 'F:/FYP Stuff/Photometry data/T-Start/'
Tmid_dir = 'F:/FYP Stuff/Photometry data/T-Middle/'
Tend_dir = 'F:/FYP Stuff/Photometry data/T-End/'

WASP11_dir = 'F:/FYP Stuff/Photometry data/WASP11/'

Mbias_33 = fits.getdata('F:/FYP Stuff/DBF night1/Bias/MasterBias/Masterbias.FITS')
Mbias_11 = fits.getdata('F:/FYP Stuff/OHP Exoplanet files 2019/Biases/Masterbias.FIT')

#function for retrieving flux values from text file
def flux_retriev(directory, file):
    flux_file = open(directory+file, 'r')
    vals = []
    for line in flux_file:
        vals.append(line.split('\n'))
        
    flux = []
    for i in range(len(vals)):
        flux.append(float(vals[i][0]))
        
    flux_file.close()
    
    return np.array(flux)

def mean_retriev(file):
    file = open('F:/FYP Stuff/Photometry data/'+file, 'r')
    bkg_mean = file.readlines()
    file.close()
    
    means = []
    
    for i in range(len(bkg_mean)):
        means.append(float(bkg_mean[i].split('\n')[0]))
        
    means=np.array(means)
    
    return means


def SNR(flux_vals, mean_vals, Mbias, ap_pixels, ann_pixels, bins):
    n_ap = ap_pixels #number of pixels in star aperture - see Flux_finder.py
    C_ap = flux_vals #Sum of pixel values in star aperture
    n_ann = ann_pixels #number of pixels in annulus - see Flux_finder.py
    
    g = 0.94*(bins**2) #gain from OHP website
    
    sig_q = 0.29 #ADUs of quantization noise from digitization of CCDs analog output, use 0.29 from handbook of astro image processing
    C_sky = mean_vals #ADUs of sky background present in every pixel
    C_dark = 0 #ADUs of dark current added to every pixel
    sig_ron = np.std(Mbias) #readout noise in electrons rms added to every pixel
    
    S_star = g*(C_ap - n_ap*(C_ap/n_ann)) #Star signal 
    
    SNR = S_star/(np.sqrt(S_star + n_ap*(1 + n_ap/n_ann)*(g*C_sky + g*C_dark +
                            sig_ron**2 + (g*sig_q)**2)))
    return SNR


#-----------------------retrieving flux of stars------------------------------#
#WASP33
st_flux = flux_retriev(Tstart_dir, 'st_flux.txt')
mid_flux = flux_retriev(Tmid_dir, 'mid_flux.txt')
end_flux = flux_retriev(Tend_dir, 'end_flux.txt')

#Pre-processed WASP33
prepro_sf = flux_retriev(Tstart_dir, 'prepro_sf.txt')
prepro_mf = flux_retriev(Tmid_dir, 'prepro_mf.txt')
prepro_ef = flux_retriev(Tend_dir, 'prepro_ef.txt')

#WASP11
WASP11_flux = flux_retriev(WASP11_dir, 'flux.txt')
#-----------------------------------------------------------------------------#

#-----------retrieving mean background counts around each star----------------#
#WASP33
st_mean = mean_retriev('T-Start/st_bkg_mean.txt')
mid_mean = mean_retriev('T-Middle/m_bkg_mean.txt')
end_mean = mean_retriev('T-End/e_bkg_mean.txt')

BD_st_mean = mean_retriev('T-Start/BDst_bkg_mean.txt')
BD_m_mean = mean_retriev('T-Middle/BDm_bkg_mean.txt')
BD_e_mean = mean_retriev('T-End/BDe_bkg_mean.txt')

#TYC 2335-839-1
TYC_st_mean = mean_retriev('T-Start/TYCst_bkg_mean.txt')
TYC_m_mean = mean_retriev('T-Middle/TYCm_bkg_mean.txt')
TYC_e_mean = mean_retriev('T-End/TYCe_bkg_mean.txt')

#Pre-processed WASP33
prepro_sm = mean_retriev('T-Start/prepro_sbm.txt')
prepro_mm = mean_retriev('T-Middle/prepro_mbm.txt')
prepro_em = mean_retriev('T-End/prepro_ebm.txt')

#WASP11
WASP11_mean = mean_retriev('WASP11/bkg_mean.txt')
#-----------------------------------------------------------------------------#

#------------------getting the SNR values for each image----------------------#
#WASP33
st_SNR = SNR(st_flux, st_mean, Mbias_33, 201, 455, 2)
m_SNR = SNR(mid_flux, mid_mean, Mbias_33, 201, 455, 2)
e_SNR = SNR(end_flux, end_mean, Mbias_33, 201, 455, 2)


#Pre-processed WASP33
prepro_sSNR = SNR(prepro_sf, prepro_sm, Mbias_33, 201, 455, 2)
prepro_mSNR = SNR(prepro_mf, prepro_mm, Mbias_33, 201, 455, 2)
prepro_eSNR = SNR(prepro_ef, prepro_em, Mbias_33, 201, 455, 2)

#WASP11
SNR_11 = SNR(WASP11_flux, WASP11_mean, Mbias_11, 1963, 2827, 1)
avg_SNR_11 = np.average(SNR_11)
#-----------------------------------------------------------------------------#

#-----getting the average SNR value for start, middle and, end of transit-----#
#WASP33
st_SNR_avg = np.average(st_SNR)
m_SNR_avg = np.average(m_SNR)
e_SNR_avg = np.average(e_SNR)

#Pre-processed WASP33
prepro_sSNR_avg = np.average(prepro_sSNR)
prepro_mSNR_avg = np.average(prepro_mSNR)
prepro_eSNR_avg = np.average(prepro_eSNR)
#-----------------------------------------------------------------------------#

Tstart_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-start/T-start_finals/'
startpics = (glob.glob(Tstart_dir+'/*.FITS'))

test_st_SNR = []
for i in range(len(startpics)):
    med_bkg = np.median(fits.getdata(startpics[i]))
    test_st_SNR.append(st_flux[i]/med_bkg)


print('\n')
print('SNR at the start of the transit (post-processing): ', st_SNR_avg)
print('SNR at the start of the transit (pre-processing): ', prepro_sSNR_avg)
print('\n')
print('SNR during the transit (post-processing): ', m_SNR_avg)
print('SNR during the transit (pre-processing): ', prepro_mSNR_avg)
print('\n')
print('SNR at the end of the transit (post-processing): ', e_SNR_avg)
print('SNR at the end of the transit (pre-processing): ', prepro_eSNR_avg)
print('\n')
print('SNR for WASP-11:', avg_SNR_11)
