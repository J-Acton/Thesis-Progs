# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 18:31:22 2023

Description:
    A prog to perform photometry on the fits file of WASP33a to determine its
    brightness.


@author: James Acton (20325303)
"""

from photutils.aperture import CircularAnnulus as cann
#from photutils.aperture import CircularAperture as cap
#from photutils.aperture import aperture_photometry as ap
from photutils.aperture import ApertureStats as apstats
from astropy.io import fits
import glob

# function for retrieving the centre coords of a star from a text file
def coord_retriev(file):
    direc = 'F:/FYP Stuff/Photometry data/'
    coords_txt = open(direc+file, 'r')
    x = []
    y = []
    vals = []
    for line in coords_txt:
        vals.append(line.split(' centre at: ')) #removing unnecessary text from text file
        
    for n in vals:
        n[-1] = n[-1].split(' ')
        
    coords_txt.close()
    
    #[file number][0 = file name | 1 = coordinates][if previous index = 1 => 0 = x | 1 = y]
    
    for i in range(len(vals)): #seperating out x and y coords
        x.append(int(vals[i][1][0]))
        y.append(int(vals[i][1][1]))
        
    return x, y

# function for calculating the mean background pixel value and writing it to a text file
def bkg_mean(x, y, pics, file): 
    #r = 8 #25 #circular aperture radius
    r_in = 12 # 40 #inner annulus radius
    r_out = 17 # 50 #outer annulus radius
    pos=[]
    for i in range(len(x)):
        pos.append([x[i], y[i]])
        
    #aps = cap(pos, r)
    anns = cann(pos, r_in, r_out)
    
    data = []
    for i in range(len(pics)):
        data.append(fits.getdata(pics[i]))
        
    stats = []
    bkg_mean = []
    
    file = open('F:/FYP Stuff/Photometry data/'+file, 'w')
    for i in range(len(pics)):
        stats.append(apstats(fits.getdata(pics[i]), anns[i]))
        bkg_mean.append(stats[i].mean)
        
    for i in range(len(bkg_mean)):
        file.write('\n'+str(bkg_mean[i]))
    file.close()
    
# retrieving centre coords for main stars and reference stars

#WASP33
xs, ys = coord_retriev('T-Start/Tstart_maxs.txt')
xm, ym = coord_retriev('T-Middle/Tmid_maxs.txt')
xe, ye = coord_retriev('T-End/Tend_maxs.txt')

#BD+36 488
xs_BD, ys_BD = coord_retriev('T-Start/stbkg_maxs.txt')
xm_BD, ym_BD = coord_retriev('T-Middle/mbkg_maxs.txt')
xe_BD, ye_BD = coord_retriev('T-End/ebkg_maxs.txt')

#TYC 2335-839-1
xs_TYC, ys_TYC = coord_retriev('T-Start/stbkg_maxs2.txt')
xm_TYC, ym_TYC = coord_retriev('T-Middle/mbkg_maxs2.txt')
xe_TYC, ye_TYC = coord_retriev('T-End/ebkg_maxs2.txt')

#TYC 2831-1138-1
xs_TYC2, ys_TYC2 = coord_retriev('T-Start/stbkg_maxs3.txt')
xm_TYC2, ym_TYC2 = coord_retriev('T-Middle/mbkg_maxs3.txt')
xe_TYC2, ye_TYC2 = coord_retriev('T-End/ebkg_maxs3.txt')

#TYC 2831-1356-1
xs_TYC3, ys_TYC3 = coord_retriev('T-Start/stbkg_maxs4.txt')
xm_TYC3, ym_TYC3 = coord_retriev('T-Middle/mbkg_maxs4.txt')
xe_TYC3, ye_TYC3 = coord_retriev('T-End/ebkg_maxs4.txt')

#WASP11
xs_11, ys_11 = coord_retriev('WASP11/maxs.txt')
LSPM_xs, LSPM_ys = coord_retriev('WASP11/LSPM_maxs.txt')

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

WASP11_pics = (glob.glob(WASP11_dir+'/*.FIT'))

#Commented out calling function so as to not overwrite existing text files
"""
#Background mean for start of transit
#bkg_mean(xs, ys, startpics, 'T-Start/st_bkg_mean.txt') #WASP33
#bkg_mean(xs_BD, ys_BD, startpics, 'T-Start/BDst_bkg_mean.txt') #BD+36 488
#bkg_mean(xs_TYC, ys_TYC, startpics, 'T-Start/TYCst_bkg_mean.txt') #TYC 2335-839-1
#bkg_mean(xs_TYC2, ys_TYC2, startpics, 'T-Start/TYC2st_bkg_mean.txt') #TYC 2831-1138-1
bkg_mean(xs_TYC3, ys_TYC3, startpics, 'T-Start/TYC3st_bkg_mean.txt') #TYC 2831-1356-1


#Background mean for middle of transit
#bkg_mean(xm, ym, midpics, 'T-Middle/m_bkg_mean.txt') #WASP33
#bkg_mean(xm_BD, ym_BD, midpics, 'T-Middle/BDm_bkg_mean.txt') #BD+36 488
#bkg_mean(xm_TYC, ym_TYC, midpics, 'T-Middle/TYCm_bkg_mean.txt') #TYC 2335-839-1
#bkg_mean(xm_TYC2, ym_TYC2, midpics, 'T-Middle/TYC2m_bkg_mean.txt') #TYC 2831-1138-1
bkg_mean(xm_TYC3, ym_TYC3, midpics, 'T-Middle/TYC3m_bkg_mean.txt') #TYC 2831-1356-1


#Background mean for end of transit
#bkg_mean(xe, ye, endpics, 'T-End/e_bkg_mean.txt') #WASP33
#bkg_mean(xe_BD, ye_BD, endpics, 'T-End/BDe_bkg_mean.txt') #BD+36 488
#bkg_mean(xe_TYC, ye_TYC, endpics, 'T-End/TYCe_bkg_mean.txt') #TYC 2335-839-1
#bkg_mean(xe_TYC2, ye_TYC2, endpics, 'T-End/TYC2e_bkg_mean.txt') #TYC 2831-1138-1
bkg_mean(xe_TYC3, ye_TYC3, endpics, 'T-End/TYC3e_bkg_mean.txt') #TYC 2831-1356-1

#Background mean for WASP11 and reference star
bkg_mean(xs_11, ys_11, WASP11_pics, 'WASP11/bkg_mean.txt')
bkg_mean(LSPM_xs, LSPM_ys, WASP11_pics, 'WASP11/LSPM_bkg_mean.txt')


#----------------Measuring background before image processing-----------------#
bkg_mean(xs, ys, prepro_sp, 'T-Start/prepro_sbm.txt')
bkg_mean(xm, ym, prepro_mp, 'T-Middle/prepro_mbm.txt')
bkg_mean(xe, ye, prepro_ep, 'T-End/prepro_ebm.txt')
"""
#-----------------------------------------------------------------------------#