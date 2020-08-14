'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    units.py: this file contains useful information about units and conversion
'''

import math
import numpy as np
from scipy import interpolate

def kpc2km():
	kpc2km = 3.086e+16
	return kpc2km

def km2kpc():
	return 1.0 / kpc2km()

def s2Myr():
	return 1.0 / Myr2s()

def Myr2s():
	Myr2s = 3.154e+13
	return Myr2s

def Myr2z(time):
	a_value = Myr2a(time)
	z_value = (1.0 / a_value) - 1.0
	return z_value

def Myr2a(time):
	# Do some interpolation
	a_fname = 'data/output_list_5Myr.txt' 
	a_file = open(a_fname, 'r')
	t_step = 5.0	# Mega-Years

	all_a = a_file.readlines()
	n_a = len(all_a)
	a_interp = np.zeros((n_a))
	t_interp = np.zeros((n_a))

	for iz in range(0, n_a):
		this_a = float(all_a[iz].rstrip()) 
		this_t = (1 + iz) * t_step
		a_interp[iz] = this_a
		t_interp[iz] = this_t

	f_t = interpolate.interp1d(t_interp, a_interp)

	a_value = f_t(time)

	return a_value

def z2Myr(z_factor):
	a_value = 1.0 / (z_factor + 1.0)
	t_value = a2Myr(a_value)
	
	return t_value

def a2Myr(a_factor):
	# Do some interpolation
	a_fname = 'data/output_list_5Myr.txt' 
	a_file = open(a_fname, 'r')
	t_step = 5.0	# Mega-Years

	all_a = a_file.readlines()
	n_a = len(all_a)
	a_interp = np.zeros((n_a))
	t_interp = np.zeros((n_a))

	for iz in range(0, n_a):
		this_a = float(all_a[iz].rstrip()) 
		this_t = (1.0 + iz)* t_step
		a_interp[iz] = this_a
		t_interp[iz] = this_t

	f_a = interpolate.interp1d(a_interp, t_interp)
	t_value = f_a(a_factor) 

	return t_value

def G():
    G = 4.302e-6            # kpc / Msun  * (km 2 / s^2)
    G_const_Mpc_Msun_s = 4.51737014558e-48
    G = G_const_Mpc_Msun_s * (s2Myr()**-2) * 1.0e+6 # 1.0e+6 kpc to Mpc conversion

    return G

def e_unit():
    e_unit=(3.2407e-2) ** 2
    return e_unit

def km2mpc():
    km2mpc = 3.24078e-20
    return km2mpc

def particle_density(n_grid, box):
    
    std_dens = np.power((n_grid / box), 3.0)

    return std_dens
