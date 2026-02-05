#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 13:42:57 2025

@author: pierre
"""


# Imports
import sys
from matplotlib import pyplot as plt
import numpy as np
import scipy
import loadandprocessdata

# Path to the python wrapper
#path_wrapper = '/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/qdyn-release-3.0.0/qdyn/'
#sys.path.append(path_wrapper)

from pyqdyn import *
import plot_functions as qdyn_plot


# Initialise the wrapperclearclear
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
    "SOLVER": 1, # Runge-Kutta-Fehlberg solver
    "V_PL": 0.0, # The fault is block on some part
    "L": L , # Length of the fault 
    "SIGMA":50.0e6,# Normal traction along the fault
    "VS":3464.0, # Shear Velocity
    "MU":32.04e9, # Shear modulus
    "FEAT_VAR_K": 0, # Coupling permeability
    "FEAT_FLUID_DIFF": 1, # Allowing for fault diffusion
    "FEAT_LOCALISATION": 0,
    "ACC": 1e-6, # Accuracy of the time solver
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
#p.run()



#%
p.read_output()





#%% Regorganize based on the number of step


def calculate_P(z,time):
    # Define all the variable 
    alpha = 0.1
    toff = 100*86400.0

    # Calculate pressure
    if time<=toff:
        P = 1.25e-6/(1e-8*0.1*np.sqrt(alpha))*(Gfunction(z,time,alpha))
    else:
        P = 1.25e-6/(1e-8*0.1*np.sqrt(alpha))*(Gfunction(z,time,alpha)-Gfunction(z,time-toff,alpha))
    

    # Return P
    return P 


def Gfunction(z,t,alpha):
    G = np.sqrt(t)*(np.exp(-np.power(z,2)/(4*alpha*t))/np.sqrt(np.pi)-np.abs(z)/np.sqrt(4*alpha*t)*scipy.special.erfc(np.abs(z)/np.sqrt(4*alpha*t)))
    return G




x_unique = p.ox["x"].unique()
sort_inds = np.argsort(x_unique)
x_unique = x_unique[sort_inds]
t_vals = np.sort(p.ox["t"].unique())
v = p.ox['v']
v=np.reshape(v,(len(t_vals),len(x_unique)))
theta = np.reshape(p.ox['theta'],(len(t_vals),len(x_unique)))
P = np.reshape(p.ox['P'],(len(t_vals),len(x_unique)))




plt.figure()

plt.semilogy(x_unique,np.transpose(v[0:-1:5,:]))



plt.figure()
plt.semilogy(t_vals/(86400*365.25),v[:,np.int16(N_half-1)])
plt.xlim((0,1))
plt.show()

#%%
plt.figure(figsize=(15, 10))
#time_id = 10


# 50 days
time_id = np.argmin(np.abs(t_vals-50*86400))
plt.plot(x_unique/1000,P[time_id,:])
P_th = calculate_P(x_unique,t_vals[time_id])
plt.plot(x_unique/1000,P_th,'--')

# 10 days
time_id = np.argmin(np.abs(t_vals-10*86400))
plt.plot(x_unique/1000,P[time_id,:])
P_th = calculate_P(x_unique,t_vals[time_id])
plt.plot(x_unique/1000,P_th,'--')
print('Time ')
#%%
# 100 days
time_id = np.argmin(np.abs(t_vals-100*86400))
plt.plot(x_unique/1000,P[time_id,:])
P_th = calculate_P(x_unique,t_vals[time_id])
plt.plot(x_unique/1000,P_th,'--')



# 200 days
time_id = np.argmin(np.abs(t_vals-200*86400))
plt.plot(x_unique/1000,P[time_id,:])
P_th = calculate_P(x_unique,t_vals[time_id])
plt.plot(x_unique/1000,P_th,'--')


plt.xlim([0,3])

plt.legend(['Calc','Th'])

plt.savefig('qdyn_diff_BP6.pdf')



#%%
x_unique = p.ox["x"].unique()
sort_inds = np.argsort(x_unique)
x_unique = x_unique[sort_inds]
t_vals = np.sort(p.ox["t"].unique())
v = p.ox['v']
v=np.reshape(v,(len(t_vals),len(x_unique)))
theta = np.reshape(p.ox['theta'],(len(t_vals),len(x_unique)))
P = np.reshape(p.ox['P'],(len(t_vals),len(x_unique)))



# Définir les chemins des fichiers
tandem_0_path = '/Users/pierre/Dropbox/PapersinPreparation/romanet2024_DFN/data/Tandem_0.0.txt'
tandem_25_path = '/Users/pierre/Dropbox/PapersinPreparation/romanet2024_DFN/data/Tandem_2.5.txt'
dublanchet_0_path = '/Users/pierre/Dropbox/PapersinPreparation/romanet2024_DFN/data/Dublanchet_0.0.txt'
dublanchet_25_path = '/Users/pierre/Dropbox/PapersinPreparation/romanet2024_DFN/data/Dublanchet_2.5.txt'

# Importer les données en ignorant les lignes d'en-tête (comme 'importdata' avec un offset)
tandem_0 = np.loadtxt(tandem_0_path, skiprows=21)
tandem_25 = np.loadtxt(tandem_25_path, skiprows=21)
dublanchet_0 = np.loadtxt(dublanchet_0_path, skiprows=23)
dublanchet_25 = np.loadtxt(dublanchet_25_path, skiprows=23)


# plt.figure()
# plt.semilogy(x_unique,v[1,:])
# plt.show()




# plt.figure()
# plt.semilogy(t_vals/(86400*365.25),v[:,np.int16(N_half-1)])
# plt.xlim((0,1))
# plt.xlabel('Time (year)')
# plt.ylabel('Slip velocity (m/s)')
# plt.title('Slip Velocity at injection point')


# plt.figure()
# plt.plot(t_vals/(86400*365.25),np.log10(theta[:,np.int16(N_half-1)]))
# plt.xlim((0,1))
# plt.xlabel('Time (year)')
# plt.ylabel('Theta (s)')
# plt.title('Slip Velocity at injection point')
#%%
skip =40

index = np.argmin(np.abs(x_unique - 2500))


plt.figure(figsize=(15, 10))
plt.subplot(2,2,1)

# Qdyn data
plt.semilogy(t_vals/(86400*365.25),v[:,np.int16(N_half-1)])

# Tandem
plt.semilogy(
    tandem_0[::skip, 0] / (86400 * 365.25),
    10 ** (tandem_0[::skip, 2]),
    'k+',
    label='Tandem 0.0'
)

# Dublanchet
plt.semilogy(
    dublanchet_0[::skip, 0] / (86400 * 365.25),
    10 ** (dublanchet_0[::skip, 2]),
    'ko',fillstyle='none',
    linewidth=1,
    label='Dublanchet 0.0'
)




plt.semilogy([100/365.25,100/365.25],[1e-12,1e-6],'k')

plt.xlim((0,1))
plt.xlabel('Time (year)')
plt.ylabel('Slip velocity (m/s)')
plt.title('Slip Velocity at injection point')



plt.subplot(2,2,2)
# Qdyn
plt.plot(t_vals/(86400*365.25),np.log10(theta[:,np.int16(N_half-1)]))



# Premier tracé : Tandem 0.0 ('k+')
plt.plot(
    tandem_0[::skip, 0] / (86400 * 365.25),
    tandem_0[::skip, 6],        # colonne 7 en MATLAB → index 6 en Python
    'k+',
    label='Tandem 0.0'
)

# Deuxième tracé : Dublanchet 0.0 ('ko')
plt.plot(
    dublanchet_0[::skip, 0] / (86400 * 365.25),
    dublanchet_0[::skip, 6],
    'ko',fillstyle='none',
    linewidth=1,
    label='Dublanchet 0.0'
)


plt.xlim((0,1))
plt.xlabel('Time (year)')
plt.ylabel('log[Theta] (log[s])')
plt.title('Theta at injection point')


plt.subplot(2,2,3)

# Qdyn data
plt.semilogy(t_vals/(86400*365.25),v[:,np.int16(index)])

# Tandem
plt.semilogy(
    tandem_25[::skip, 0] / (86400 * 365.25),
    10 ** (tandem_25[::skip, 2]),
    'k+',
    label='Tandem 0.0'
)

# Dublanchet
plt.semilogy(
    dublanchet_25[::skip, 0] / (86400 * 365.25),
    10 ** (dublanchet_25[::skip, 2]),
    'ko',fillstyle='none',
    linewidth=1,
    label='Dublanchet 0.0'
)




plt.semilogy([100/365.25,100/365.25],[1e-12,1e-6],'k')

plt.xlim((0,1))
plt.xlabel('Time (year)')
plt.ylabel('Slip velocity (m/s)')
plt.title('Slip Velocity at injection point')



plt.subplot(2,2,4)
# Qdyn
plt.plot(t_vals/(86400*365.25),np.log10(theta[:,index]))



# Premier tracé : Tandem 0.0 ('k+')
plt.plot(
    tandem_25[::skip, 0] / (86400 * 365.25),
    tandem_25[::skip, 6],        # colonne 7 en MATLAB → index 6 en Python
    'k+',
    label='Tandem 0.0'
)

# Deuxième tracé : Dublanchet 0.0 ('ko')
plt.plot(
    dublanchet_25[::skip, 0] / (86400 * 365.25),
    dublanchet_25[::skip, 6],
    'ko',fillstyle='none',
    linewidth=1,
    label='Dublanchet 0.0'
)


plt.xlim((0,1))
plt.xlabel('Time (year)')
plt.ylabel('log[Theta] (log[s])')
plt.title('Theta at injection point')


plt.savefig('qdyn_BP6.pdf')
# plt.show()




