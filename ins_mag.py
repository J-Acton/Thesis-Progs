# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 16:19:23 2024

Description:


@author: James Acton (20325303)
"""
import numpy as np
import matplotlib.pyplot as plt


#Defining the directories for the fully processed fits files
Tstart_dir = 'F:/FYP Stuff/Photometry data/T-Start/'
Tmid_dir = 'F:/FYP Stuff/Photometry data/T-Middle/'
Tend_dir = 'F:/FYP Stuff/Photometry data/T-End/'


def Ins_mag(flux, Z, n_ap, n_ann, t):
    C_ap = flux
    
    m = -2.5*np.log10((C_ap - n_ap*(C_ap/n_ann))/t) + Z
    
    return m

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
    
"""
#--------------------------------WASP33---------------------------------------#
st_flux = flux_retriev(Tstart_dir, 'st_flux.txt')
mid_flux = flux_retriev(Tmid_dir, 'mid_flux.txt')
end_flux = flux_retriev(Tend_dir, 'end_flux.txt')

#WASP-33 b transit start and end times
T_start_33 = -81/60
T_end_33 = 82/60

st_mag = Ins_mag(st_flux, 21.525, 201, 455, 0.5)
m_mag = Ins_mag(mid_flux, 21.525, 201, 455, 0.5)
e_mag = Ins_mag(end_flux, 21.525, 201, 455, 0.5)

xs = np.linspace(-(86/60), -(3/4), 1063, endpoint=True)
xm = np.linspace((13/60), (28/60), 394, endpoint=True)
xe = np.linspace((73/60), (94/60), 572, endpoint=True)

plt.title('WASP-33 - Raw Instrumental Magnitude', fontsize=26)
plt.scatter(xs, st_mag, marker='o', color='#00128c', s=10, alpha=0.6)
plt.scatter(xm, m_mag, marker='o', color='#00128c', s=10, alpha=0.6)
plt.scatter(xe, e_mag, marker='o', color='#00128c', s=10, alpha=0.6)
plt.axvline(T_start_33, linestyle='-', color = 'g', 
            label='Transit Start & End')
plt.axvline(T_end_33, linestyle='-', color = 'g')
plt.ylabel('Raw Instrumental Magnitude', fontsize=24)
plt.xlabel('Hours From Centre of Transit', fontsize=24)
plt.gca().invert_yaxis()
plt.xlim(-2,2)
plt.xticks(size=20)
plt.yticks(size=20)
plt.grid()
plt.legend(loc=1, fontsize=22)
#-----------------------------------------------------------------------------#
"""
#---------------------------------WASP11--------------------------------------#
WASP11_dir = 'F:/FYP Stuff/Photometry data/WASP11/'

T_start_11 = -(80/60)

WASP11_flux = flux_retriev(WASP11_dir, 'flux.txt')

mag_11 = Ins_mag(WASP11_flux, 22.4, 1963, 2827, 40)

x_11 = np.linspace(-(108/60), (58/60), 217, endpoint=True)

plt.title('WASP-11 - Raw Instrumental Magnitude', fontsize=26)
plt.scatter(x_11, mag_11, marker='o', color='#a40000', s=10, alpha=0.6)
plt.axvline(T_start_11, linestyle='-', color = 'b', 
            label='Transit Start')
plt.ylabel('Raw Instrumental Magnitude', fontsize=24)
plt.xlabel('Hours From Centre of Transit', fontsize=24)
plt.xlim(-2,2)
plt.xticks(size=20)
plt.yticks(size=20)
plt.gca().invert_yaxis()
plt.grid()
plt.legend(loc=1, fontsize=22)




