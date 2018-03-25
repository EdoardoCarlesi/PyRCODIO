from utils import *
from units import *

#def halo_trajectory():

# This assumes that halos_zp1 already contains all the closest halos to halo_z, the candidate for which we are looking for a progenitor
def find_progenitor(halo_z, halos_zp1):
	n_zp1 = 0
	guess_x = d_v * d_t
	dist = 0 
	mass = 0
	

def new_x(old_x, vel, dMyrs):
	new_x = [0.0] * 3

	if (old_x[0] > 500.):
		facMpc = 1.0
	else:
		facMpc = 1000.0

	totT = dMyrs * s2Myr()

	for ix in range(0, 3):	
		new_x[ix] = old_x[ix] * facMpc - totT * vel[ix] / km2kpc()

	return new_x



