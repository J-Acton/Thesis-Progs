# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 17:37:34 2024

Description:


@author: James Acton (20325303)
"""

import numpy as np
import matplotlib.pyplot as plt

from pytransit import QuadraticModel

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

npv = 50
#seed(0)

times_sc = np.linspace(-1/12, 1/12, 6223)  # Short cadence time stamps
times_lc = np.linspace(1/12, 1/12,  100)  # Long cadence time stamps

#-------------------------Transit parameters----------------------------------#
k = 0.110811506 # ratio of planet radius to star radius
t0 = 0.0 # zero epoch
p = 1.21986967 # period of transit (days)
a = 3.558 # semi-major axis (in terms of planet radii)
i = 1.530654 # orbital inclination
e = 0.0 # eccentricity
w  = 0.5*np.pi # argument of periastron, 90 degrees for e=0
ldc = np.array([0.3, 0.21, 0.1]) #limb darkening coeficients
#-----------------------------------------------------------------------------#

normal = np.random.normal
uniform = np.random.uniform

ks    = normal(k, 0.002, npv)
t0s   = normal(t0, 0.004, npv)
ps    = normal(p, 0.01, npv)
smas  = normal(a, 0.1, npv)
incs  = uniform(0.49*np.pi, 0.5*np.pi, npv)
es    = uniform(0, 0.25, npv)
ws    = uniform(0, 2*np.pi, npv)
ldc   = uniform(0.1, 0.5, size=(npv, 2))

tm = QuadraticModel()

tm.set_data(times_sc)

def plot_transits(tm, ldc, fmt='g'):
    fig, axs = plt.subplots(1, 1 , figsize = (13,3), constrained_layout=True, sharey=True)

    # Evaluate for scalar parameters
    flux = tm.evaluate(k, ldc[1], t0, p, a, i, e, w)
    axs.plot(tm.time, flux, fmt)
    axs.tick_params(axis='x', labelsize=20)
    axs.tick_params(axis='y', labelsize=20)
    axs.scatter(xs, start_flux_norm, s=5, color='b', alpha=0.35)
    axs.scatter(xm, mid_flux_norm, s=5, color='b', alpha=0.35)
    axs.scatter(xe, end_flux_norm, s=5, color='b', alpha=0.35)
    axs.set_title('WASP-33 Transit Model', fontsize=26)
    axs.set_xlim(min(tm.time), max(tm.time))
    #axs.set_ylim(0.7, 1.05)
    axs.grid()
    axs.set_ylabel('Normalised Flux', fontsize=24)
    axs.set_xlabel('Days From Centre of Transit', fontsize=24)

    # Evaluate for npv parameters in parallel
    #flux = tm.evaluate(ks, ldc, t0s, ps, smas, incs, es, ws)
    #axs[1].plot(tm.time, flux.T, 'k', alpha=0.2)
    #axs[1].scatter(xs, start_flux_norm, s=1, color='b', alpha=0.25)
    #axs[1].scatter(xm, mid_flux_norm, s=1, color='b', alpha=0.25)
    #axs[1].scatter(xe, end_flux_norm, s=1, color='b', alpha=0.25)
    #axs[1].set_title('Parameter vector')

    plt.setp(axs, ylabel='Normalised flux')
    plt.setp(axs, xlabel='Time [days from centre of transit]', xlim=(-1/12,1/12))#, ylim=(0.85, 1.02))
   # plt.grid()
    
    tm.set_data(times_sc)
    

#Defining the directories for the fully processed fits files
Tstart_dir = 'F:/FYP Stuff/Photometry data/T-Start/'
Tmid_dir = 'F:/FYP Stuff/Photometry data/T-Middle/'
Tend_dir = 'F:/FYP Stuff/Photometry data/T-End/'    

#retrieving flux for WASP33 (Mag 8.14)
start_flux = flux_retriev(Tstart_dir, 'st_flux.txt')
mid_flux = flux_retriev(Tmid_dir, 'mid_flux.txt')
end_flux = flux_retriev(Tend_dir, 'end_flux.txt')

flux_total = np.concatenate((start_flux, mid_flux, end_flux))
max_flux = max(flux_total)

start_flux_norm = start_flux/max_flux
mid_flux_norm = mid_flux/max_flux
end_flux_norm = end_flux/max_flux

xs = np.linspace(-(86/(60*24)), -(3/(4*(24))), len(start_flux), endpoint=True)
xm = np.linspace((13/(60*24)), (28/(60*24)), len(mid_flux), endpoint=True)
xe = np.linspace((73/(60*24)), (94/(60*24)), len(end_flux), endpoint=True)

mod_xs = np.linspace(xs[0], xs[-1], len(start_flux), endpoint=True)
mod_xm = np.linspace(xm[0], xm[-1], len(mid_flux), endpoint=True)
mod_xe = np.linspace(xe[0], xe[-1], len(end_flux), endpoint=True)

half_length = 0.08333333

mod_flux_st = (np.around(6223*((half_length+mod_xs)/(2*(half_length))))).astype(int)
mod_flux_m = (np.around(6223*((half_length+mod_xm)/(2*(half_length))))).astype(int)
mod_flux_e = (np.around(6223*((half_length+mod_xe)/(2*(half_length))))).astype(int)

plot_transits(tm, ldc)

#calculating percentage dip of flux according to transit
model_flux = tm.evaluate(k, ldc[1], t0, p, a, i, e, w)
print((1-min(model_flux))*100)
