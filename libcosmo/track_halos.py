from utils import *
from units import *

#def halo_trajectory():

# This assumes that halos_zp1 already contains all the closest halos to halo_z, the candidate for which we are looking for a progenitor
def find_progenitor(halo_z, halos_zp1, timeStep):
	n_zp1 = len(halos_zp1)

	# The progenitor halo will be evaluated using these four quantities as proxy
	angl = 0
	dist = 0 
	mass = 0
	velo = 0
	
	# No likely progenitors 
	if n_zp1 == 0:
		#TODO enable to return placeholder
		return -1
	else:	
		guess_x = backward_x(halo_z.x, halo_z.v, timeStep)
		
		for izp1 in range(0, n_zp1):
			this_x = halos_zp1[izp1].x
			this_v = halos_zp1[izp1].v
			this_m = halos_zp1[izp1].m
			this_a = angle()
			eval_this = this_a


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


