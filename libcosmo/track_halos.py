from utils import *
from units import *

#def halo_trajectory():

# This assumes that halos_zp1 already contains all the closest halos to halo_z, the candidate for which we are looking for a progenitor
def find_progenitor(halo_z, halos_zp1):
	n_zp1 = 0
	guess_x = d_v * d_t
	dist = 0 
	mass = 0
	

def new_x(old_x, vel, d_z):
	new_x = [0.0] * 3
	# Convert d_z to time interval

	if (old_x[0.0] > 500.):
		facMpc = 1.0
	else:
		facMpc = 1000.0

	d_t = 250. # Myrs
	d_t *= units.s2Myr()

	for ix in range(0, 3):	
		new_x[ix] = old_x[ix] * facMpc - d_t * vel[ix] / units.km2kpc()

	return new_x



