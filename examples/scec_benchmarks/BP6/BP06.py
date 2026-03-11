#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 12:22:43 2025

@author: pierre
"""


# Imports
import numpy as np
from matplotlib import pyplot as plt

from qdyn import qdyn

# Initialise the wrapperclear
p = qdyn()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Set BP06 benchmark
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set up size of the fault
ds = 10 # Discretisation size
N = 2**12 # Number of elements
N_half = N/2
L = ds*N # Size of the fault


# Set up initial theta
theta_ini = 4e-3/1e-6*np.exp(7.0/5.0*np.log(2.0*1e-6/1e-12*np.sinh(29.2e6/(0.007*50e6)))-0.6/0.005)

# General simulation parameters
set_dict = {
    "FRICTION_MODEL": "RSF", # Rate and state friction
    "MESHDIM": 1,   # 1D antiplane fault in 2D medium
    "FINITE": 1, # Make it a finite fault loaded by constant velocity on the side
    "SOLVER": 2, # Runge-Kutta-Fehlberg solver
    "V_PL": 0.0, # The fault is block on some part
    "L": L , # Length of the fault 
    "SIGMA":50.0e6,# Normal traction along the fault
    "VS":3464.0, # Shear Velocity
    "MU":32.04e9, # Shear modulus
    "FEAT_VAR_K": 0, # Coupling permeability
    "FEAT_FLUID_DIFF": 1, # Allowing for fault diffusion
    "FEAT_LOCALISATION": 0,
    "ACC": 1e-7, # Accuracy of the time solver
    "N": N, # Number of the elements (must be a power of 2)
    "NSTOP": 0, # Stop at t= TMAX
    "TMAX": 2*365.25*86400 # Finish time of the simulation
}
   
# RSF parameters
set_dict_RSF = {
    "RNS_LAW": 0,
    "THETA_LAW": 1, # Aging law
    "A": 0.007,
    "B": 0.005,
    "DC": 0.004,
    "MU_SS": 0.6,
    "V_SS": 1e-6,
    "V_0":1e-12,
    "TH_0": theta_ini, # Initial state,
}

 # Fluid diffusion model 
set_dict_fluid_diff = {
     "RHOF": 1e3,						# Density Fluid
     "BETA": 1e-8,						# Bulk compressibility (fluid + pore) [1/Pa]
     "ETA": 1e-3,						# Fluid dynamic viscosity [Pa s]
     "PHI": 0.1,                       # Porosity 
     "PERMEABILITY": 1e-13,            # Permeability
     "nb_source": 1,                   # Number of source (only one for now)
     "rate_injection":  [1.25e-6],       # Injection rate (m/s)
     "index_injection": [N_half],          # Indec of the fault where injection is performed (for FORTRAN NOT PYTHON)
     "t_injection_beg": [0.0],           # Beginning time of injection
     "t_injection_end": [100.0*86400.0],  # End time of injection
     "P_a":0.0
 }
 
 # Variable permeability
set_dict_var_k = {
    "Lk": 1.0,
    "Tk": 1000,
    "kmin": 1e-15,      
    "kmax": 1e-13,
    "Snk": 1e6,
 }


# Merge dictionnaries
set_dict["SET_DICT_RSF"] = set_dict_RSF
set_dict["SET_DICT_FLUID_DIFF"] = set_dict_fluid_diff
set_dict["SET_DICT_VAR_K"] = set_dict_var_k
set_dict["NTOUT_OX"] = 10        # Save output every N steps

# Set the simulation
p.settings(set_dict)

# Make the mesh
p.render_mesh()

# Write input
p.write_input()

# Run the simulation
p.run()

#%% Regorganize based on the number of step
# Read output
p.read_output()

x_unique = p.ox["x"].unique()
sort_inds = np.argsort(x_unique)
x_unique = x_unique[sort_inds]
t_vals = np.sort(p.ox["t"].unique())
v = p.ox['v']
v=np.reshape(v,(len(t_vals),len(x_unique)))
theta = np.reshape(p.ox['theta'],(len(t_vals),len(x_unique)))
P = np.reshape(p.ox['P'],(len(t_vals),len(x_unique)))

plt.figure()
plt.semilogy(x_unique,v[1,:])
plt.show()

plt.figure()
plt.semilogy(t_vals/(86400*365.25),v[:,np.int16(N_half-1)])
plt.xlim((0,1))
plt.xlabel('Time (year)')
plt.ylabel('Slip velocity (m/s)')
plt.title('Slip Velocity at injection point')

plt.figure()
plt.plot(t_vals/(86400*365.25),np.log10(theta[:,np.int16(N_half-1)]))
plt.xlim((0,1))
plt.xlabel('Time (year)')
plt.ylabel('Theta (s)')
plt.title('Slip Velocity at injection point')
