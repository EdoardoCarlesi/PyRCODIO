import math
import numpy as np
from scipy import interpolate

def km2kpc():
	km2kpc = 3.086e+16
	return km2kpc

def s2Myr():
	s2Myr = 3.154e+13
	return s2Myr

def Myr2z(time):
	# Do some interpolation
	a_fname = 'data/output_list_5Myr.txt' 
	a_file = open(a_fname, 'r')
	t_step = 5.0	# Mega-Years

	all_a = a_file.readlines()
	n_a = len(all_a)
	z_interp = np.zeros((n_a))
	t_interp = np.zeros((n_a))

	for iz in range(0, n_a):
		this_z = 1.0 / (float(all_a[iz].rstrip())) - 1.0
		this_t = (1 + iz) * t_step
		z_interp[iz] = this_z
		t_interp[iz] = this_t

	f_t = interpolate.interp1d(t_interp, z_interp)
	z_value = f_t(time)

	return z_value

def z2Myr(redshift):
	# Do some interpolation
	a_fname = 'data/output_list_5Myr.txt' 
	a_file = open(a_fname, 'r')
	t_step = 5.0	# Mega-Years

	all_a = a_file.readlines()
	n_a = len(all_a)
	z_interp = np.zeros((n_a))
	t_interp = np.zeros((n_a))

	#print 'Total points: ', n_a, n_a * 0.005

	for iz in range(0, n_a):
		this_z = 1.0 / (float(all_a[iz].rstrip())) - 1.0
		this_t = (1.0 + iz)* t_step
		z_interp[iz] = this_z
		t_interp[iz] = this_t

	f_z = interpolate.interp1d(z_interp, t_interp)
	t_value = f_z(redshift) 

	return t_value



