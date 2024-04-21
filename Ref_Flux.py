# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:17:23 2024

Description:


@author: James Acton (20325303)
"""

import numpy as np
import matplotlib.pyplot as plt

#Defining the directories for the fully processed fits files
Tstart_dir = 'F:/FYP Stuff/Photometry data/T-Start/'
Tmid_dir = 'F:/FYP Stuff/Photometry data/T-Middle/'
Tend_dir = 'F:/FYP Stuff/Photometry data/T-End/'

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

def plotting(y1, y2, y3, star, c):
    
    
    xs = np.linspace(-(86/60), -(3/4), 1063, endpoint=True)
    xm = np.linspace((13/60), (28/60), 394, endpoint=True)
    xe = np.linspace((73/60), (94/60), 572, endpoint=True)
    
    y = np.concatenate((y1, y2, y3))
    
    plt.scatter(xs, y1/(max(y)), marker='o', color=c, s=5, alpha=0.35)
    plt.scatter(xm, y2/(max(y)), marker='o', color=c, s=5, alpha=0.35)
    plt.scatter(xe, y3/(max(y)), marker='o', color=c, s=5, alpha=0.35)
    
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.title('Normalised Flux Measured of '+star, fontsize=26)
    plt.ylabel('Normalised Flux', fontsize=24)
    plt.xlabel('Hours From Centre of Transit', fontsize=24)
    plt.grid()

    


#retrieving flux for BD+36 488 (Mag 9.37)
BD_st_flux = flux_retriev(Tstart_dir, 'BD_st_flux.txt')
BD_mid_flux = flux_retriev(Tmid_dir, 'BD_mid_flux.txt')
BD_end_flux = flux_retriev(Tend_dir, 'BD_end_flux.txt')

#retrieving flux for TYC 2335-839-1 (Mag 11.08)
TYC_st_flux = flux_retriev(Tstart_dir, 'TYC_st_flux.txt')
TYC_mid_flux = flux_retriev(Tmid_dir, 'TYC_mid_flux.txt')
TYC_end_flux = flux_retriev(Tend_dir, 'TYC_end_flux.txt')


#retrieving flux for TYC 2831-1138-1 (Mag 11.27)
TYC2_st_flux = flux_retriev(Tstart_dir, 'TYC2_st_flux.txt')
TYC2_mid_flux = flux_retriev(Tmid_dir, 'TYC2_mid_flux.txt')
TYC2_end_flux = flux_retriev(Tend_dir, 'TYC2_end_flux.txt')


#retrieving flux for TYC 2831-1356-1 (Mag 11.48)
TYC3_st_flux = flux_retriev(Tstart_dir, 'TYC3_st_flux.txt')
TYC3_mid_flux = flux_retriev(Tmid_dir, 'TYC3_mid_flux.txt')
TYC3_end_flux = flux_retriev(Tend_dir, 'TYC3_end_flux.txt')


#-----------------------------------------------------------------------------#

fig1 = plt.figure()
plotting(BD_st_flux, BD_mid_flux, BD_end_flux, 'BD+36 488', 'r')


fig2 = plt.figure()
plotting(TYC_st_flux, TYC_mid_flux, TYC_end_flux, 'TYC 2335-839-1', 'b')


fig3 = plt.figure()
plotting(TYC2_st_flux, TYC2_mid_flux, TYC2_end_flux, 'TYC 2831-1138-1', 'g')


fig4 = plt.figure()
plotting(TYC3_st_flux, TYC3_mid_flux, TYC3_end_flux, 'TYC 2831-1356-1', 'k')
