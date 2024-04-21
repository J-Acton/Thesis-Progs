# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 16:36:17 2024

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

times_sc = np.linspace(-1/12, 1/12, 1000)  # Short cadence time stamps
times_lc = np.linspace(1/12, 1/12,  100)  # Long cadence time stamps

#-------------------------Transit parameters----------------------------------#
k = 0.129588823 # ratio of planet radius to star radius
t0 = 0.0 # zero epoch
p = 3.7 # period of transit (days)
a = 11.649 # semi-major axis (in terms of planet radii)
i = (898/1800)*np.pi # orbital inclination
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

def plot_transits(tm, ldc, fmt='b'):
    fig, axs = plt.subplots(1, 1 , figsize = (13,3), constrained_layout=True, sharey=True)

    # Evaluate for scalar parameters
    flux = tm.evaluate(k, ldc[1], t0, p, a, i, e, w)
    axs.plot(tm.time, flux, fmt)
    axs.tick_params(axis='x', labelsize=20)
    axs.tick_params(axis='y', labelsize=20)
    axs.scatter(xs, WASP11_flux_norm, s=5, color='r', alpha=1)
    axs.set_title('WASP-11 Transit Model', fontsize=26)
    axs.set_ylim(0.92, 1.02)
    axs.set_ylabel('Normalised Flux', fontsize=24)
    axs.set_xlabel('Time [days from centre of transit]', fontsize=24)
    axs.set_xlim(min(tm.time), max(tm.time))
    axs.grid()

    # Evaluate for npv parameters in parallel
    #flux = tm.evaluate(ks, ldc, t0s, ps, smas, incs, es, ws)
    #axs[1].plot(tm.time, flux.T, 'k', alpha=0.2)
    #axs[1].scatter(xs, WASP11_flux_norm, s=1, color='b', alpha=0.25)
    #axs[1].set_title('Parameter vector')

    plt.setp(axs, ylabel='Normalised flux')
    #plt.setp(axs, xlabel='Time [days from centre of transit]', xlim=(-1/12,1/12), ylim=(0.85, 1.02))
    #plt.grid()
    
    tm.set_data(times_sc)
    
    
WASP11_dir = 'F:/FYP Stuff/Photometry data/WASP11/'

#retrieving flux for WASP11 (Mag 11.57)
WASP11_flux = flux_retriev(WASP11_dir, 'flux.txt')
WASP11_flux_norm = WASP11_flux/(max(WASP11_flux))

xs = np.linspace(-3/40, 19/480, len(WASP11_flux), endpoint=True)

plot_transits(tm, ldc)

model_flux = tm.evaluate(k, ldc[1], t0, p, a, i, e, w)
print((1-min(model_flux))*100)