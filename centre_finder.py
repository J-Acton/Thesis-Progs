# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 09:22:08 2023

Description:
    A prog to determine the coords for the centre of WASP33a in each of the 
    fits files taken.

@author: James Acton (20325303)
"""

from astropy.io import fits
import glob


#Defining the directories for the fully processed fits files
Tstart_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-start/T-start_finals/'
Tmid_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-middle/T-mid_finals/'
Tend_dir = 'F:/FYP Stuff/SCJA WASP33b/N1/T-end/T-end_finals/'

WASP11_dir = 'F:/FYP Stuff/OHP Exoplanet files 2019/Processed_Data/'

#Making list of T-start, middle and end fits files
startpics = (glob.glob(Tstart_dir+'/*.FITS'))
midpics = (glob.glob(Tmid_dir+'/*.FITS'))
endpics = (glob.glob(Tend_dir+'/*.FITS'))

WASP11_pics = (glob.glob(WASP11_dir+'/*.FIT'))

#-------------------------Progress bar stuff----------------------------------#
sttot = len(startpics)
midtot = len(midpics)
endtot = len(endpics)

WASP11_tot = len(WASP11_pics)

#-----------------------------------------------------------------------------#

#----------------Locating centre of annulus for each image--------------------#
def centre_finder(xll, yll, xul, yul, pics, file, tot):
    Phot_direc = 'F:/FYP Stuff/Photometry data/'
    
    file = open(Phot_direc+file, 'w') # opening text file to write centre coords
    for k in range(len(pics)): # for loop to find centre coords in each image
        max_val = 0
        progress = (k/tot)
        arrow = '=' * int(30*progress)
        spaces = ' ' * (30 - len(arrow))
        print(f'\r[{arrow}{spaces}] {int(progress * 100)}%', end='', flush=True)
        data = fits.getdata(pics[k])
        for j in range(len(data)):
            for i in range(len(data[j])):
                if i > xll and i < xul and j > yll and j < yul:
                    if data[j][i] > max_val:
                        max_val = data[j][i]
                        xval = i
                        yval = j
        file.write(str('\n'+str(pics[k])+' centre at: '+str(xval)+' '+str(yval))) # writing centre coords to text file
    file.close()
#-----------------------------------------------------------------------------#

#---------------------------Box limit values----------------------------------#
#T-start
xll_WASP_st, yll_WASP_st, xul_WASP_st, yul_WASP_st = 500, 370, 560, 440
xll_BD_st, yll_BD_st, xul_BD_st, yul_BD_st = 730, 166, 810, 216
xll_TYC_st, yll_TYC_st, xul_TYC_st, yul_TYC_st = 176, 733, 232, 778
xll_TYC2_st, yll_TYC2_st, xul_TYC2_st, yul_TYC2_st = 480, 430, 570, 490
xll_TYC3_st, yll_TYC3_st, xul_TYC3_st, yul_TYC3_st = 820, 400, 930, 500

#T-mid
xll_WASP_m, yll_WASP_m, xul_WASP_m, yul_WASP_m = 391, 413, 491, 513
xll_BD_m, yll_BD_m, xul_BD_m, yul_BD_m = 660, 235, 694, 265
xll_TYC_m, yll_TYC_m, xul_TYC_m, yul_TYC_m = 82, 797, 142, 843
xll_TYC2_m, yll_TYC2_m, xul_TYC2_m, yul_TYC2_m = 370, 490, 480, 560
xll_TYC3_m, yll_TYC3_m, xul_TYC3_m, yul_TYC3_m = 730, 460, 840, 560

#T-end
xll_WASP_e, yll_WASP_e, xul_WASP_e, yul_WASP_e = 387, 416, 487, 516
xll_BD_e, yll_BD_e, xul_BD_e, yul_BD_e = 653, 220, 705, 270
xll_TYC_e, yll_TYC_e, xul_TYC_e, yul_TYC_e = 85, 790, 134, 840
xll_TYC2_e, yll_TYC2_e, xul_TYC2_e, yul_TYC2_e = 370, 480, 480, 550
xll_TYC3_e, yll_TYC3_e, xul_TYC3_e, yul_TYC3_e = 710, 450, 830, 560

#WASP11
xll_WASP11, yll_WASP11, xul_WASP11, yul_WASP11 = 990, 800, 1208, 940
xll_LSPM, yll_LSPM, xul_LSPM, yul_LSPM = 0, 477, 75, 563
#-----------------------------------------------------------------------------#

#---------------------------Start of transit----------------------------------#
#centre_finder(xll_WASP_st, yll_WASP_st, xul_WASP_st, yul_WASP_st, startpics, 
#              'T-start/Tstart_maxs.txt', len(startpics)) #WASP33

#centre_finder(xll_BD_st, yll_BD_st, xul_BD_st, yul_BD_st, startpics, 
#             'T-start/stbkg_maxs.txt', len(startpics)) #BD+36 488

#centre_finder(xll_TYC_st, yll_TYC_st, xul_TYC_st, yul_TYC_st, startpics, 
#              'T-start/stbkg_maxs2.txt', len(startpics)) #TYC 2335-839-1

#centre_finder(xll_TYC2_st, yll_TYC2_st, xul_TYC2_st, yul_TYC2_st, startpics,
 #             'T-start/stbkg_maxs3.txt', len(startpics)) #TYC 2831-1138-1

#centre_finder(xll_TYC3_st, yll_TYC3_st, xul_TYC3_st, yul_TYC3_st, startpics,
 #             'T-start/stbkg_maxs4.txt', len(startpics)) #TYC 2831-1356-1

#-----------------------------------------------------------------------------#


#-------------------------Middle of transit-----------------------------------#
#centre_finder(xll_WASP_m, yll_WASP_m, xul_WASP_m, yul_WASP_m, midpics,
 #             'T-middle/Tmid_maxs.txt', len(midpics)) #WASP33

#centre_finder(xll_BD_m, yll_BD_m, xul_BD_m, yul_BD_m, midpics,
 #             'T-middle/mbkg_maxs.txt', len(midpics)) #BD+36 488

#centre_finder(xll_TYC_m, yll_TYC_m, xul_TYC_m, yul_TYC_m, midpics,
 #             'T-middle/mbkg_maxs2.txt', len(midpics)) #TYC 2335-839-1
 
#centre_finder(xll_TYC2_m, yll_TYC2_m, xul_TYC2_m, yul_TYC2_m, midpics, 
 #             'T-middle/mbkg_maxs3.txt', len(midpics)) #TYC 2831-1138-1

#centre_finder(xll_TYC3_m, yll_TYC3_m, xul_TYC3_m, yul_TYC3_m, midpics, 
 #             'T-middle/mbkg_maxs4.txt', len(midpics)) #TYC 2831-1356-1
#-----------------------------------------------------------------------------#


#---------------------------End of transit------------------------------------#
#centre_finder(xll_WASP_e, yll_WASP_e, xul_WASP_e, yul_WASP_e, endpics,
 #             'T-end/Tend_maxs.txt', len(endpics)) #WASP33

#centre_finder(xll_BD_e, yll_BD_e, xul_BD_e, yul_BD_e, endpics,
 #             'T-end/ebkg_maxs', len(endpics)) #BD+36 488

#centre_finder(xll_TYC_e, yll_TYC_e, xul_TYC_e, yul_TYC_e, endpics,
 #             'T-end/ebkg_maxs2.txt', len(endpics)) #TYC 2335-839-1

#centre_finder(xll_TYC2_e, yll_TYC2_e, xul_TYC2_e, yul_TYC2_e, endpics,
 #             'T-end/ebkg_maxs3.txt', len(endpics)) #TYC 2831-1138-1

#centre_finder(xll_TYC3_e, yll_TYC3_e, xul_TYC3_e, yul_TYC3_e, endpics,
 #             'T-end/ebkg_maxs4.txt', len(endpics)) #TYC 2831-1356-1
#-----------------------------------------------------------------------------#


#----------------------------------WASP11-------------------------------------#
#centre_finder(xll_WASP11, yll_WASP11, xul_WASP11, yul_WASP11, WASP11_pics, 'WASP11/maxs.txt', WASP11_tot)
centre_finder(xll_LSPM, yll_LSPM, xul_LSPM, yul_LSPM, WASP11_pics, 'WASP11/LSPM_maxs.txt', WASP11_tot)
#-----------------------------------------------------------------------------#
