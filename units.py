'''
    Python Routines for COsmology and Data I/O (PyRCODIO) v0.2
    Edoardo Carlesi (2020) 
    ecarlesi83@gmail.com

    units.py: this file contains useful information about units and conversion
'''

import math
import numpy as np
from scipy import interpolate


def kpc2km():
    ''' Simple conversion '''

    kpc2km = 3.086e+16
    return kpc2km


def km2kpc():
    ''' Simple conversion '''
    
    return 1.0 / kpc2km()


def Myr2s():
    ''' Simple conversion '''
    
    Myr2s = 3.154e+13
    return Myr2s


def s2Myr():
    ''' Simple conversion '''
    
    return 1.0 / Myr2s()


def Myr2a(time):
    '''
    Reads from a pre-computed table with timesteps of 5 Myrs, then interpolates to the desired output
    WARNING this is computed once within a fixed cosmology so it is not reliable when using other cosmological parameters
    '''

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


def Myr2z(time):
    ''' Simple conversion '''
 
    a_value = Myr2a(time)
    z_value = (1.0 / a_value) - 1.0
    
    return z_value


def a2Myr(a_factor):
    '''
    Reads from a pre-computed table with timesteps of 5 Myrs, then interpolates to the desired output
    WARNING this is computed once within a fixed cosmology so it is not reliable when using other cosmological parameters
    '''
	
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


def z2Myr(z_factor):
    ''' Simple conversion '''

    a_value = 1.0 / (z_factor + 1.0)
    t_value = a2Myr(a_value)
	
    return t_value


def G():
    ''' Returns the gravitational constant '''

    G = 4.302e-6            # kpc / Msun  * (km 2 / s^2)
    G_const_Mpc_Msun_s = 4.51737014558e-48
    G = G_const_Mpc_Msun_s * (s2Myr()**-2) * 1.0e+6 # 1.0e+6 kpc to Mpc conversion

    return G


def e_unit():
    ''' Rescaled energy units '''

    e_unit=(3.2407e-2) ** 2

    return e_unit


def km2mpc():
    ''' Simple conversion '''

    km2mpc = 3.24078e-20
    
    return km2mpc


def particle_density(grid=None, box=None):
    ''' Very idiotic density calculation '''

    std_dens = np.power((grid / box), 3.0)

    return std_dens


def rad2deg():
    ''' This returns a conversion factor from radians to degrees '''

    coeff = 180.0 / np.pi

    return coeff


if __name__ == '__main__':
    ''' Use this for testing local functions '''
    pass



