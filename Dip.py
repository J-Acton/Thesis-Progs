# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 15:29:21 2024

Description:
    A prog to calculate the magnitude of target star using the reference star 
    BD+36 488

@author: James Acton (20325303)
"""

import numpy as np
import matplotlib.pyplot as plt


#Defining the directories for the fully processed fits files
Tstart_dir = 'F:/FYP Stuff/Photometry data/T-Start/'
Tmid_dir = 'F:/FYP Stuff/Photometry data/T-Middle/'
Tend_dir = 'F:/FYP Stuff/Photometry data/T-End/'

WASP11_dir = 'F:/FYP Stuff/Photometry data/WASP11/'

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

#function for calculating the difference in magnitude between two stars given their flux
def mag_diff(F1, F2):
    del_m = -2.5*np.log10(F1/F2)
    return del_m

#function for estimating the instrumental magnitude of a star
def Ins_mag_calculator(ref_mag, mag_diff):
    Ins_mag = mag_diff + ref_mag
    return Ins_mag

"""
def recorder(directory, file, data):
    file = open(directory+file, 'w')
    file.write('\n'+str(data))
    file.close()
"""
#----------------------retrieving flux values---------------------------------#

#retrieving flux for WASP33 (Mag 8.14)
start_flux = flux_retriev(Tstart_dir, 'st_flux.txt')
mid_flux = flux_retriev(Tmid_dir, 'mid_flux.txt')
end_flux = flux_retriev(Tend_dir, 'end_flux.txt')
WASP33_flux = np.concatenate((start_flux, mid_flux, end_flux))

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


#retrieving flux for WASP11 (Mag 11.57)
WASP11_flux = flux_retriev(WASP11_dir, 'flux.txt')
LSPM_flux = flux_retriev(WASP11_dir, 'LSPM_flux.txt')
#-----------------------------------------------------------------------------#

#-------Setting catalogue values for the V-band magnitudes of each star-------#
Mc_BD = 9.37 #BD+36 488
Mc_TYC = 11.08 #TYC 2335-839-1
Mc_TYC2 = 11.27 #TYC 2831-1138-1
Mc_TYC3 = 11.48 # TYC 2831-1356-1

Mc_WASP33 = 8.3 #WASP33
#-----------------------------------------------------------------------------#
#------------------------x arrays for plotting--------------------------------#

#2023
xs = np.linspace(-(86/60), -(3/4), 1063, endpoint=True)
xm = np.linspace((13/60), (28/60), 394, endpoint=True)
xe = np.linspace((73/60), (94/60), 572, endpoint=True)

x = np.concatenate((xs, xm, xe))

#WASP-33 b transit start and end times
T_start = -81/60
T_end = 82/60

#2019
x_11 = np.linspace(-(108/60), (58/60), 217, endpoint=True)
#-----------------------------------------------------------------------------#

#-----Calculating differences in magnitude between WASP-33 and references-----#
#BD+36 488
ys_BD = mag_diff(start_flux, BD_st_flux)
ym_BD = mag_diff(mid_flux, BD_mid_flux)
ye_BD = mag_diff(end_flux, BD_end_flux)

#TYC 2335-839-1
ys_TYC = mag_diff(start_flux, TYC_st_flux)
ym_TYC = mag_diff(mid_flux, TYC_mid_flux)
ye_TYC = mag_diff(end_flux, TYC_end_flux)

#TYC 2831-1138-1
ys_TYC2 = mag_diff(start_flux, TYC2_st_flux)
ym_TYC2 = mag_diff(mid_flux, TYC2_mid_flux)
ye_TYC2 = mag_diff(end_flux, TYC2_end_flux)

#TYC 2831-1356-1
ys_TYC3 = mag_diff(start_flux, TYC3_st_flux)
ym_TYC3 = mag_diff(mid_flux, TYC3_mid_flux)
ye_TYC3 = mag_diff(end_flux, TYC3_end_flux)
#-----------------------------------------------------------------------------#

#--------comparing the magnitude differences between reference stars----------#
ys_BDvsTYC = mag_diff(BD_st_flux, TYC_st_flux)
ym_BDvsTYC = mag_diff(BD_mid_flux, TYC_mid_flux)
ye_BDvsTYC = mag_diff(BD_end_flux, TYC_end_flux)
BDvsTYC = np.concatenate((ys_BDvsTYC, ym_BDvsTYC, ye_BDvsTYC))

ys_BDvsTYC2 = mag_diff(BD_st_flux, TYC2_st_flux)
ym_BDvsTYC2 = mag_diff(BD_mid_flux, TYC2_mid_flux)
ye_BDvsTYC2 = mag_diff(BD_end_flux, TYC2_end_flux)
BDvsTYC2 = np.concatenate((ys_BDvsTYC2, ym_BDvsTYC2, ye_BDvsTYC2))

ys_BDvsTYC3 = mag_diff(BD_st_flux, TYC3_st_flux)
ym_BDvsTYC3 = mag_diff(BD_mid_flux, TYC3_mid_flux)
ye_BDvsTYC3 = mag_diff(BD_end_flux, TYC3_end_flux)
BDvsTYC3 = np.concatenate((ys_BDvsTYC3, ym_BDvsTYC3, ye_BDvsTYC3))

ys_TYCvsTYC2 = mag_diff(TYC_st_flux, TYC2_st_flux)
ym_TYCvsTYC2 = mag_diff(TYC_mid_flux, TYC2_mid_flux)
ye_TYCvsTYC2 = mag_diff(TYC_end_flux, TYC2_end_flux)
TYCvsTYC2 = np.concatenate((ys_TYCvsTYC2, ym_TYCvsTYC2, ye_TYCvsTYC2))

ys_TYCvsTYC3 = mag_diff(TYC_st_flux, TYC3_st_flux)
ym_TYCvsTYC3 = mag_diff(TYC_mid_flux, TYC3_mid_flux)
ye_TYCvsTYC3 = mag_diff(TYC_end_flux, TYC3_end_flux)
TYCvsTYC3 = np.concatenate((ys_TYCvsTYC3, ym_TYCvsTYC3, ye_TYCvsTYC3))

ys_TYC2vsTYC3 = mag_diff(TYC2_st_flux, TYC3_st_flux)
ym_TYC2vsTYC3 = mag_diff(TYC2_mid_flux, TYC3_mid_flux)
ye_TYC2vsTYC3 = mag_diff(TYC2_end_flux, TYC3_end_flux)
TYC2vsTYC3 = np.concatenate((ys_TYC2vsTYC3, ym_TYC2vsTYC3, ye_TYC2vsTYC3))
#-----------------------------------------------------------------------------#

#-------------------Instrumental magnitude estimates--------------------------#
Imag_st_1 = Ins_mag_calculator(Mc_BD, ys_BD)
Imag_st_2 = Ins_mag_calculator(Mc_TYC, ys_TYC)
Imag_st_3 = Ins_mag_calculator(Mc_TYC2, ys_TYC2)
Imag_st_4 = Ins_mag_calculator(Mc_TYC3, ys_TYC3)
Imag_st_avg = (Imag_st_1 + Imag_st_2 + Imag_st_3 + Imag_st_4)/4

#recorder(Tstart_dir, 'st_mag.txt', Imag_st_avg)


Imag_m_1 = Ins_mag_calculator(Mc_BD, ym_BD)
Imag_m_2 = Ins_mag_calculator(Mc_TYC, ym_TYC)
Imag_m_3 = Ins_mag_calculator(Mc_TYC2, ym_TYC2)
Imag_m_4 = Ins_mag_calculator(Mc_TYC3, ym_TYC3)
Imag_m_avg = (Imag_m_1 + Imag_m_2 + Imag_m_3 + Imag_m_4)/4

#recorder(Tmid_dir, 'm_mag.txt', Imag_m_avg)


Imag_e_1 = Ins_mag_calculator(Mc_BD, ye_BD)
Imag_e_2 = Ins_mag_calculator(Mc_TYC, ye_TYC)
Imag_e_3 = Ins_mag_calculator(Mc_TYC2, ye_TYC2)
Imag_e_4 = Ins_mag_calculator(Mc_TYC3, ye_TYC3)
Imag_e_avg = (Imag_e_1 + Imag_e_2 + Imag_e_3 + Imag_e_4)/4

#recorder(Tend_dir, 'e_mag.txt', Imag_e_avg)

#-----------------------------------------------------------------------------#

#--------Average y value for plots before, during, and after transit----------#
#Outside - start
def start_avg(st_data):
    st_sum = 0
    for i in range(len(st_data)):
        if i < (5/41)*len(st_data):
            st_sum += st_data[i]
    st_avg = st_sum/((5/41)*len(st_data))
    return st_avg

BDavg_st = start_avg(ys_BD)
TYCavg_st = start_avg(ys_TYC)
TYC2avg_st = start_avg(ys_TYC2)
TYC3avg_st = start_avg(ys_TYC3)
   
Ins_mag_st_sum = 0
for i in range(len(Imag_st_avg)):
    if i < (5/41)*len(Imag_st_avg):
        Ins_mag_st_sum += (Imag_st_avg[i])
Ins_mag_st_avg = Ins_mag_st_sum/((5/41)*len(Imag_st_avg))


#middle
BDavg_m = np.average(ym_BD)
TYCavg_m = np.average(ym_TYC)
TYC2avg_m = np.average(ym_TYC2)
TYC3avg_m = np.average(ym_TYC3)


#outside - end
def end_avg(e_data):
    e_sum = 0
    for i in range(len(e_data)):
        if i > (3/7)*len(e_data):
            e_sum += e_data[i]
    e_avg = e_sum/((4/7)*len(e_data))
    return e_avg

BDavg_e = end_avg(ye_BD)
TYCavg_e = end_avg(ye_TYC)
TYC2avg_e = end_avg(ye_TYC2)
TYC3avg_e = end_avg(ye_TYC3)


Ins_mag_e_sum = 0
for i in range(len(Imag_e_avg)):
    if i > (3/7)*len(Imag_e_avg):
        Ins_mag_e_sum += (Imag_e_avg[i])
Ins_mag_e_avg = Ins_mag_e_sum/((4/7)*len(Imag_e_avg))
#-----------------------------------------------------------------------------#

#getting the average mag diff values outside the transit
BDavg = (BDavg_st+BDavg_e)/2
TYCavg = (TYCavg_st+TYCavg_e)/2
TYC2avg = (TYC2avg_st+TYC2avg_e)/2
TYC3avg = (TYC3avg_st+TYC3avg_e)/2
#avg_comp = (np.average(ys_BDvsTYC)+np.average(ye_BDvsTYC))/2
Ins_mag_avg = (Ins_mag_st_avg + Ins_mag_e_avg)/2
#-----------------------------------------------------------------------------#

#--------------Plotting mag diffs of target with reference stars--------------#
fig, ax = plt.subplots()

# Setting up subplots
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

# start of transit scatters
ax1.scatter(xs, ys_BD, marker='o', color='b', s=1, alpha=0.15)
ax2.scatter(xs, ys_TYC, marker='o', color='b', s=1, alpha=0.15)
ax3.scatter(xs, ys_TYC2, marker='o', color='b', s=1, alpha=0.15)
ax4.scatter(xs, ys_TYC3, marker='o', color='b', s=1, alpha=0.15)

# middle of transit scatters
ax1.scatter(xm, ym_BD, marker='o', color='b', s=1, alpha=0.15)
ax2.scatter(xm, ym_TYC, marker='o', color='b', s=1, alpha=0.15)
ax3.scatter(xm, ym_TYC2, marker='o', color='b', s=1, alpha=0.15)
ax4.scatter(xm, ym_TYC3, marker='o', color='b', s=1, alpha=0.15)

# end of transit scatters
ax1.scatter(xe, ye_BD, marker='o', color='b', s=1, alpha=0.15)
ax2.scatter(xe, ye_TYC, marker='o', color='b', s=1, alpha=0.15)
ax3.scatter(xe, ye_TYC2, marker='o', color='b', s=1, alpha=0.15)
ax4.scatter(xe, ye_TYC3, marker='o', color='b', s=1, alpha=0.15)

#--------------------formatting the different graphs--------------------------#


ax1.axvline(T_start, linestyle='-', color = 'g', 
            label='Transit Start & End')
ax1.axvline(T_end, linestyle='-', color = 'g')
ax1.set_title('Difference in Apparent Magnitude Between WASP-33 & BD+36 488')
ax1.set_ylim(BDavg - 0.07, BDavg + 0.07)
ax1.axhline(BDavg, linestyle='--', color='k', 
            label = 'Avg Mag Diff Outside Transit')
ax1.axhline(BDavg_m, linestyle='--', color='r', 
            label='Avg Mag Diff During Transit')
ax1.set_ylabel('Magnitude Difference')
#ax1.set_xlabel('Hours From Centre of Transit')
ax1.invert_yaxis()
ax1.legend(loc=1)
ax1.grid()
#ax1.set_xlim((0, max(xe)))
ax1.set_xlim(-2,2)


ax2.axvline(T_start, linestyle='-', color = 'g')
ax2.axvline(T_end, linestyle='-', color = 'g')
ax2.set_title('Difference in apparent magnitude between WASP-33 & TYC 2335-839-1')
ax2.set_ylim(TYCavg - 0.1, TYCavg + 0.1)
ax2.axhline(TYCavg, linestyle='--', color='k')
ax2.axhline(TYCavg_m, linestyle='--', color='r')
ax2.invert_yaxis()
ax2.grid()

ax2.set_xlim(-2,2)

ax3.axvline(T_start, linestyle='-', color = 'g')
ax3.axvline(T_end, linestyle='-', color = 'g')
ax3.set_title('Difference in apparent magnitude between WASP-33 and TYC 2831-1138-1')
ax3.set_ylim(TYC2avg - 0.1, TYC2avg + 0.1)
ax3.axhline(TYC2avg, linestyle='--', color='k')
ax3.axhline(TYC2avg_m, linestyle='--', color='r')
ax3.invert_yaxis()
ax3.set_ylabel('Magnitude difference')
ax3.set_xlabel('Hours from centre of Transit')
#ax3.set_xlabel('"Time"')
ax3.grid()
ax3.set_xlim(-2,2)

ax4.axvline(T_start, linestyle='-', color = 'g')
ax4.axvline(T_end, linestyle='-', color = 'g')
ax4.set_title('Difference in apparent magnitude between WASP-33 & TYC 2831-1356-1')
ax4.axhline(TYC3avg, linestyle='--', color='k')
ax4.axhline(TYC3avg_m, linestyle='--', color='r')
ax4.set_ylabel('Magnitude difference')
ax4.set_ylim(TYC3avg - 0.1, TYC3avg + 0.1)
ax4.invert_yaxis()
#ax4.set_xlabel('"Time"')
ax4.set_xlabel('Hours from centre of Transit')
ax4.grid()
#ax4.set_xlim((0, max(xe)))
ax4.set_xlim(-2, 2)
#-----------------------------------------------------------------------------#

#-------------------plotting mag diffs of reference stars---------------------#
fig2, ax = plt.subplots()

ax5 = plt.subplot(231)
ax6 = plt.subplot(232)
ax7 = plt.subplot(233)
ax8 = plt.subplot(234)
ax9 = plt.subplot(235)
ax10 = plt.subplot(236)

ax5.scatter(x, BDvsTYC, marker='o', color='g', s=1, alpha=0.25)
ax6.scatter(x, BDvsTYC2, marker='o', color='g', s=1, alpha=0.25)
ax7.scatter(x, BDvsTYC3, marker='o', color='g', s=1, alpha=0.25)
ax8.scatter(x, TYCvsTYC2, marker='o', color='g', s=1, alpha=0.25)
ax9.scatter(x, TYCvsTYC3, marker='o', color='g', s=1, alpha=0.25)
ax10.scatter(x, TYC2vsTYC3, marker='o', color='g', s=1, alpha=0.25)

    
ax5.set_title('Mag diff - BD+36 488 vs TYC 2335-839-1')
#ax5.set_xlim(0, max(xe))
ax5.set_xlim(-2,2)
ax5.set_ylim(-1.65, -1.35)
ax5.grid()
ax5.set_ylabel('Magnitude Difference')
ax5.invert_yaxis()

ax6.set_title('Mag diff - BD+36 488 vs TYC 2831-1138-1')
#ax6.set_xlim(0, max(xe))
ax6.set_xlim(-2,2)
#ax6.set_ylim(-2.3, -2.1)
ax6.grid()
#ax6.set_ylabel('Magnitude Difference')
ax6.invert_yaxis

ax7.set_title('Mag diff - BD+36 488 vs TYC 2831-1356-1')
#ax7.set_xlim(0, max(xe))
ax7.set_xlim(-2,2)
ax7.set_ylim(-2.35, -2.05)
ax7.grid()
ax7.invert_yaxis

ax8.set_title('Mag diff - TYC 2335-839-1 vs TYC 2831-1138-1')
#ax8.set_xlim(0, max(xe))
ax8.set_xlim(-2,2)
ax8.set_ylim(-0.8, -0.5)
ax8.set_ylabel('Magnitude Difference')
ax8.set_xlabel('Hours from centre of Transit')
ax8.grid()
ax8.invert_yaxis

ax9.set_title('Mag diff - TYC 2335-839-1 vs TYC 2831-1356-1')
#ax9.set_xlim(0, max(xe))
ax9.set_xlim(-2,2)
ax9.set_ylim(-0.85, -0.55)
ax9.set_xlabel('Hours from centre of Transit')
ax9.grid()
ax9.invert_yaxis()

ax10.set_title('Mag_diff - TYC 2831-1138-1 vs TYC 2831-1356-1')
#ax10.set_xlim(0, max(xe))
ax10.set_xlabel('Hours from centre of Transit')
ax10.set_xlim(-2,2)
ax10.set_ylim(-0.15, 0.15)
ax10.grid()
ax10.invert_yaxis
#-----------------------------------------------------------------------------#

#-----------plotting average instrumental mag of target star------------------#
fig3 = plt.figure()

plt.scatter(xs, Imag_st_avg, marker='o', color='b', s=5, alpha=0.35)
plt.scatter(xm, Imag_m_avg, marker='o', color='b', s=5, alpha=0.35)
plt.scatter(xe, Imag_e_avg, marker='o', color='b', s=5, alpha=0.35)


plt.xticks(size=20)
plt.yticks(size=20)
plt.title('Instrumental Magnitude Estimate', fontsize=26)
plt.axhline(Ins_mag_avg, linestyle='--', color='k', label='Avg Mag Outside Transit')
plt.axhline(np.average(Imag_m_avg), linestyle='--', color='r', label='Avg Mag During Transit')
plt.axvline(T_start, linestyle='-', color = 'g', 
            label='Transit Start & End')
plt.axvline(T_end, linestyle='-', color = 'g')
plt.ylabel('Instrumental Magnitude', fontsize=24)
#plt.xlabel('"Time"')
plt.xlabel('Hours From Centre of Transit', fontsize=24)
#plt.xlim(0, max(xe))
plt.xlim(-2,2)
plt.ylim(Ins_mag_avg - 0.02, Ins_mag_avg + 0.02)
#plt.ylim(Ins_mag_avg - 0.1, Ins_mag_avg + 0.1)
plt.grid()
plt.gca().invert_yaxis()
plt.legend(loc=1, fontsize=22)
#-----------------------------------------------------------------------------#

#-----------------------plotting flux of WASP33-------------------------------#
fig4 = plt.figure()

plt.xticks(size=20)
plt.yticks(size=20)
plt.scatter(xs, (start_flux/(max(WASP33_flux))), marker='o', color='b', s=10, alpha=1)
plt.scatter(xm, (mid_flux/(max(WASP33_flux))), marker='o', color='b', s=10, alpha=1)
plt.scatter(xe, (end_flux/(max(WASP33_flux))), marker='o', color='b', s=10, alpha=1)
plt.title('Measured Flux of WASP-33', fontsize=26)
plt.axvline(T_start, linestyle='-', color='g', label='Transit Start and End')
plt.axvline(T_end, linestyle='-', color='g')
#plt.axhline(np.average(WASP33_flux)/max(WASP33_flux))
plt.ylabel('Normalised Flux', fontsize=24)
plt.xlabel('Hours From Centre of Transit', fontsize=24)
plt.grid()
plt.legend(fontsize=22)
#-----------------------------------------------------------------------------#

#------------------------plotting flux of WASP11------------------------------#
fig5 = plt.figure()

plt.xticks(size=20)
plt.yticks(size=20)
plt.scatter(x_11, (WASP11_flux/(max(WASP11_flux))), marker='o', color='r', s=10, alpha=1)
plt.title('Measured Flux of WASP-11 - 2019 Archival Data', fontsize=26)
plt.axvline(-(80/60), linestyle='-', color='b', label='transit start')
#plt.axhline(np.average(WASP11_flux)/max(WASP11_flux))
#plt.axhline(np.average(WASP11_st), linestyle='--', color='k', label='avg outside transit')
#plt.axhline(np.average(WASP11_m), linestyle='--', color='grey', label='avg inside transit')
plt.ylabel('Normalised Flux', fontsize=24)
plt.xlabel('Hours From Centre of Transit', fontsize=24)
#plt.xlim(0, max(x_11))
#plt.ylim(0.975, 1.0)
plt.grid()
plt.legend(fontsize=22)


#-----------------------------------------------------------------------------#

#--------------------plotting flux of LSPM J0309+3042-------------------------#
fig6 = plt.figure()

plt.xticks(size=20)
plt.yticks(size=20)
plt.scatter(x_11, (LSPM_flux/(max(LSPM_flux))), marker='o', color='m', s=10, alpha=1)
plt.title('Measured Flux of LSPM J0309+3042 - 2019 Reference Star', fontsize=26)
plt.ylabel('Normalised Flux', fontsize=24)
plt.xlabel('Hours From Centre of Transit', fontsize=24)
#plt.ylim(0.75, 1.0)
plt.grid()
#-----------------------------------------------------------------------------#

#-----------------\Delta m between 11 and LSPM--------------------------------#
Delm_2019 = mag_diff(WASP11_flux, LSPM_flux)

Delm_2019_st = []
for i in range(37):
    Delm_2019_st.append(Delm_2019[i])
    
Delm_2019_m = []
for i in range(37, 217):
    Delm_2019_m.append(Delm_2019[i])
    
Delm_2019_avg = (np.average(Delm_2019_st) - np.average(Delm_2019_m))
print(Delm_2019_avg)

fig7 = plt.figure()

plt.xticks(size=20)
plt.yticks(size=20)
plt.scatter(x_11, Delm_2019, marker='o', color='#D18700', s=10, alpha=1)
plt.axhline(np.average(Delm_2019_st), linestyle='--', color='k', label='Avg Mag Diff Before Transit')
plt.axhline(np.average(Delm_2019_m), linestyle='--', color='r', label='Avg Mag Diff During Transit')
plt.title('Magnitude Difference Between WASP-11 and LSPM J0309+3042', fontsize=26)
plt.ylabel('$\Delta m$', fontsize=24)
plt.xlabel('Hours From Centre of Transit', fontsize=24)
plt.axvline(-(80/60), linestyle='-', color='b', label='transit start')
plt.ylim(-4.6, -3.6)
plt.grid()
plt.gca().invert_yaxis()
plt.legend(fontsize=22)

#-----------------------------------------------------------------------------#


def mag_diff_print(Ref_star, avg, transit_avg):
    statement = print('\nThe average change in mag diff between WASP33 and',
                      'reference star ', Ref_star, 'is: ', 
                      '{:.5f}'.format(avg-transit_avg))
    return statement

mag_diff_print('BD+36 488', BDavg, BDavg_m)
mag_diff_print('TYC 2335-839-1', TYCavg, TYCavg_m)
mag_diff_print('TYC 2831-1138-1', TYC2avg, TYC2avg_m)
mag_diff_print('TYC 2831-1356-1', TYC3avg, TYC3avg_m)

print(start_flux[68])
print(mid_flux[207])
print(end_flux[568])
