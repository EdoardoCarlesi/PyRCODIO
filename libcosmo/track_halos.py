from utils import *
from units import *
from find_halos import *
import numpy as np
import math

#def halo_trajectory():

# This assumes that halos_zp1 already contains all the closest halos to halo_z, the candidate for which we are looking for a progenitor
def find_progenitor(halo_z, halos_all_zp1, timeStep):
	guess_x = backward_x(halo_z.x, halo_z.v, timeStep)
	
	r_prog =  halo_z.r + 2 * distance(halo_z.x, guess_x)
	#print 'Looking for halos aroud a R=%.3f at %s ' % (r_prog, guess_x)
	halos_zp1 = find_halos_point(guess_x, halos_all_zp1, r_prog)
	
	n_zp1 = len(halos_zp1)

	# The progenitor halo will be evaluated using these four quantities as proxy
	fac_a = 1.5
	fac_d = 0.01 
	fac_m = 2.0
	fac_v = 1.0
		
	pos = halo_z.x
	vel = halo_z.v
	mass = halo_z.m
	eval0 = 10.0

	# No likely progenitors 
	if n_zp1 == 0:
		#TODO enable to return placeholder
		return -1
	else:	
		for izp1 in range(0, n_zp1):
			this_x = halos_zp1[izp1].x
			this_v = halos_zp1[izp1].v
			this_m = halos_zp1[izp1].m

			#print mass, vel, pos
			#print halos_zp1[izp1].info()

			this_a = 1.0 - abs_val(angle(this_v, vel))
			this_d = abs_val(distance(guess_x, this_x) / distance(pos, this_x))
			this_m = abs_val(np.log10(mass/this_m))
			this_v = abs_val(distance(this_v, vel)/module(vel))

			eval_this = fac_a * this_a + fac_m * this_m + fac_d * this_d + fac_v * this_v

			if eval_this < eval0:
			#	print 'Eval:%.3f, Angle: %.3f Dist: %.3f Mass: %.3f, Vel: %.3f' % (eval_this, this_a, this_d, this_m, this_v)
			#	print 'NewEval', eval_this, eval0
				eval0 = eval_this
				best_id = izp1

		return halos_zp1[best_id]


# TODO add/subtract hubble expansion somehow!
def backward_x(old_x, vel, dMyrs):
	new_x = [0.0] * 3

	if (old_x[0] > 500.):
		facMpc = 1.0
	else:
		facMpc = 1000.0

	totT = dMyrs * s2Myr()

	for ix in range(0, 3):	
		new_x[ix] = old_x[ix] * facMpc - totT * vel[ix] / km2kpc()

	return new_x

# TODO add/subtract hubble expansion somehow!
def forward_x(old_x, vel, dMyrs):
	new_x = [0.0] * 3

	if (old_x[0] > 500.):
		facMpc = 1.0
	else:
		facMpc = 1000.0

	totT = dMyrs * s2Myr()

	for ix in range(0, 3):	
		new_x[ix] = old_x[ix] * facMpc + totT * vel[ix] / km2kpc()

	return new_x


# Extrapolate the mass of an object at z+1 based on its z value (useful for placeholder-halos)
def mass_z(m0, z):
		
	a = 0.5
	b = 0.5

	# Check a real formula
	m_z = m0 * pow((1 + z), a) * math.exp(b * z)
	return m_z


