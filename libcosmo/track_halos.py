from utils import *
from units import *
from find_halos import *
import numpy as np
import math


# This assumes that halos_zp1 already contains all the closest halos to halo_z, the candidate for which we are looking for a progenitor
def find_progenitor(halo_z, halos_all_zp1, aFactor, timeStep):
	r_prog = 0.0
	guess_x = backward_x(halo_z.x, halo_z.v, aFactor, timeStep)
	
	#r_prog =  1000. #halo_z.distance(guess_x)
	#r_prog =  halo_z.distance(guess_x)

	r_prog = module(halo_z.v) * timeStep * s2Myr() / km2kpc() / aFactor	# In comoving units
	#print 'r_prog ', halo_z.r, distance(halo_z.x, guess_x), guess_x
	#print 'Looking for halos aroud a R=%.3f (Rphys=%.3f) at %s ' % (r_prog, r_prog * aFactor, guess_x)
	halos_zp1 = find_halos_point(guess_x, halos_all_zp1, r_prog)
	
	n_zp1 = len(halos_zp1)

	# The progenitor halo will be evaluated using these four quantities as proxy
	fac_a = 0.5
	fac_d = 1.0 # Distance from the expected position (in v * dT units)
	fac_r = 1.0 # Distance from the parent halo (in v * dT units)
	fac_m = 2.0
	fac_v = 0.05
	fac_x = 3.0 # Angle between the velocity of the halo and the displacement of the new candidate, to check that no "parallel" object has been 
		    # incorrectly identified
	
	max_eval = 3.5
	
	pos = halo_z.x
	vel = halo_z.v
	mass = halo_z.m
	eval0 = 10.0

	# No likely progenitors 
	if n_zp1 == 0:
		place_halo = Halo()
		place_halo.x = guess_x		
		# This is VERY crude... just taking mass and velocity from previous halo...
		place_halo.v = vel		
		place_halo.m = mass

		return place_halo
	else:	
		for izp1 in range(0, n_zp1):
			halo_x = halos_zp1[izp1].x
			halo_v = halos_zp1[izp1].v
			halo_m = halos_zp1[izp1].m

			#print mass, vel, pos
			#print halos_zp1[izp1].info()
			
			direction = vec_subt(pos, halo_x) 
			
			this_x = 1.0 - abs_val(angle(direction, halo_v))
			this_a = 1.0 - abs_val(angle(halo_v, vel))
			this_d = distance(guess_x, halo_x) / r_prog #distance(pos, this_x))
			this_r = distance(pos, halo_x) / r_prog #distance(pos, this_x))
			this_m = abs_val(np.log10(mass/halo_m))
			this_v = abs_val(distance(halo_v, vel)/module(vel))

			eval_this = fac_a * this_a + fac_m * this_m + fac_d * this_d + fac_v * this_v + fac_r * this_r + this_x * fac_x

			if eval_this < eval0:
				eval0 = eval_this
				best_id = izp1
			#	print 'Eval:%.3f, Angle: %.3f D: %.3f, R: %.3f, Mass: %.3f, Vel: %.3f, Direction: %.3f' % \
			#	(eval_this, this_a, this_d, this_r, this_m, this_v, this_x)
			#	print 'Mnew:%.3e, Mold: %.3e **** Vnew: %s, Vold: %s' % (halo_m, mass, halo_v, vel)
			#	print 'NewEval', eval_this, eval0
		
		# If eval0 is too large then we better replace the halo with a placeholder and see what happens at the next step
		if eval0 > max_eval:
			#print 'Returning halo placeholder...'
			place_halo = Halo()
			place_halo.x = guess_x		
			#print 'Pos: ', guess_x
			# This is VERY crude... just taking mass and velocity from previous halo...
			place_halo.v = vel		
			place_halo.m = mass

			return place_halo
		else:
			return halos_zp1[best_id]


# TODO add/subtract hubble expansion somehow!	- Old_a is the expansion factor at the time where the halo is in old_x position
def backward_x(this_x, this_vel, this_a, d_myrs):
	new_x = [0.0] * 3

	if (this_x[0] > 500.):
		facMpc = 1.0
	else:
		facMpc = 1000.0

	new_t = a2Myr(this_a) - d_myrs
	new_a = Myr2a(new_t)
	tot_t = d_myrs * s2Myr()
	delta_a = abs_val(this_a - new_a)

	#print 'BackwardX: ', this_a, d_myrs, new_a, new_t, tot_t, delta_a, d_myrs

	for ix in range(0, 3):	
		new_x[ix] = this_x[ix] * facMpc - (this_vel[ix] * tot_t) / km2kpc() * 0.677
		#new_x[ix] = - (delta_a / this_a) * (old_x[ix] - new_x[ix])
		#new_x[ix] = old_x[ix] * facMpc - totT * vel[ix] / km2kpc()
		#new_x[ix] = old_x[ix] * facMpc - vel[ix] / km2kpc() * totT * (1 + deltaA)
		# Correct for the Hubble expansion
		#new_x[ix] += (old_x[ix] - new_x[ix]) * deltaA / newA
		#new_x[ix] *= deltaA (old_x[ix] - new_x[ix]) * deltaA / newA

	return new_x



def forward_x(old_x, vel, oldA, dMyrs):
	# Just changing the velocity sign
	new_x = backward_x(old_x, -vel, oldA, dMyrs)

	return new_x


# Extrapolate the mass of an object at z+1 based on its z value (useful for placeholder-halos)
def mass_z(m0, z):
		
	a = 0.5
	b = 0.5

	# Check a real formula
	m_z = m0 * pow((1 + z), a) * math.exp(b * z)
	return m_z


