from find_halos import *
from utils import *
from units import *
from halo import *

import numpy as np
import math

'''
	Computes the merger tree for a given number of halos. This is generally the two main LG halos at z=0 and their satellites
	Subhaloes are stored separately if the trace_subs boolean variable is set to true. Then ids_subs contains the ids of the main haloes
	whose subhaloes need to be tracked at each step, using a factor r_subs * Rvir to identify them
'''
def merger_tree(end_snap, ini_snap, min_common, main_halos, main_parts, main_ids, settings, track_subs, ids_subs, r_subs):
	zs = ahf_redshifts(base_path)
	ss = ahf_snapshots(base_path)

	n_halos = len(main_halos)

	# Final list containing HaloThroughZ objects
	all_halo_z = []

	# These are simply all the haloes within a given radius at each step
	all_subs_z = []		

	# Temporary lists that allow intra-snapshot comparisons - only for the haloes being tracked
	old_main_halo = []
	old_main_part = []
	old_main_ids = []
	
	# Initialize halo_z list members 
	for i_halo in range(0, n_halos):
		halo_z = HaloThroughZ(end_snap - ini_snap)
		all_halo_z.append(halo_z)
		#subs_z = SubHaloThroughZ(end_snap - ini_snap)
		#all_subs_z.append(subs_z)

	# Now loop over all the particle and halo files
	for i_snap in range(end_snap, ini_snap, -1):
		# TODO this assumes one MPI task only !!!!!
		(this_part_file, this_halo_file) = settings.get_ahf_files(i_snap, 0)

		this_all_halo = read_ahf(this_halo_file)
		(this_all_ids, this_all_parts) = read_particles(this_part_file)
	
		# Time steps, redshift and expansion factors
		this_z = float(zs[i_snap])
		this_t = z2Myr(this_z)
		this_a = 1.0 / (1.0 + this_z)

		# This is the first step
		if i_snap == end_snap:

			# Loop on the main halos to be tracked (usually the two LG main members)
			for i_main in range(0, n_halos):
				this_halo = main_halos[i_main]
				this_part = main_parts[i_main]
				this_ids  = main_ids[i_main]
				#print i_snap, i_main, this_halo.info()

				all_halo_z[i_main].add_step(this_halo, this_z)

				# Save all the halos for the next step
				old_main_halo.append(this_halo)
				old_main_part.append(this_part)
				old_main_ids.append(this_ids)
			
				if track_subs == True:
					(is_there, entries) = is_there(this_ids, ids_subs)

					if is_there == True:
						sub_halos = find
						this_subs = SubHalos(this_halo, sub_halos, )
						all_subs_z.append(this_subs)

			# Do this once the loop on all the halos is finished
			old_t = this_t

		# Normal steps after the first one used for initialization
		else:	
			# Time step is deltaT, it might be a fixed quantity but we don't know in advance so we compute it every time
			time_step = abs(old_t - this_t)	
	
			for i_main in range(0, n_halos):
				#print i_main, old_main_part
				#print i_main, old_main_ids
				(halo_progenitors, index_progenitors, token_halo) = find_progenitors(old_main_halo[i_main], old_main_part[i_main],\
							this_all_halo, this_all_parts, min_common, this_a, time_step)

				# We only take the main progenitor, all the progenitors are sorted by merit
				this_halo = halo_progenitors[0]			
				this_index = index_progenitors[0]			

				# If no likely progenitor has been found but only a token halo has been returned, then save the old particle list
				# also for the next step
				if token_halo == True:
					this_part = old_part[i_main]
					this_ids = old_ids[i_main]
				else:
					# Find the particle ID list of the main progenitor
					this_part = this_all_parts[this_index]
					this_ids = this_all_ids[this_index]
				
				# Now save the particles and halo for the next step
				old_main_halo[i_main] = this_halo
				old_main_part[i_main] = this_part
				old_main_ids[i_main] = this_ids
			
				# Store the halo into the general structure
				all_halo_z[i_main].add_step(this_halo, this_z)
				
			old_t = this_t

	return all_halo_z



def find_halo_id(idnum, ahf_halos):
	n_halo = len(ahf_halos)
	i_halo = 0
	h_found = False	

	while (i_halo < n_halo) and (h_found == False):
		if (ahf_halos[i_halo].ID == idnum):
			h_found = True	
		else:
			i_halo += 1

	if h_found == True:
		return i_halo
	else:
		print 'Halo ID: %ld not found.' % idnum
		return -1



# This assumes that halos_zp1 already contains all the closest halos to halo_z, the candidate for which we are looking for a progenitor
def find_progenitors(halo_z, part_z, halos_all_zp1, part_all_zp1, min_common, aFactor, timeStep):
	r_prog = 0.0
	guess_x = backward_x(halo_z.x, halo_z.v, aFactor, timeStep)
	r_prog = module(halo_z.v) * timeStep * s2Myr() / km2kpc() / aFactor	# In comoving units
	#print 'r_prog ', halo_z.r, distance(halo_z.x, guess_x), guess_x
	#print 'Looking for halos aroud a R=%.3f (Rphys=%.3f) at %s ' % (r_prog, r_prog * aFactor, guess_x)
	halos_zp1 = find_halos_point(guess_x, halos_all_zp1, r_prog)
	
	token_halo = False
	id_progenitors = []
	progenitors = []
	n_progenitors = 0
	n_zp1 = len(halos_zp1)
	
	for i_prog in range(0, n_zp1):
		this_id = halos_zp1[i_prog].ID
		this_index = find_halo_id(this_id, halos_all_zp1)
		this_npart = halos_all_zp1[this_index].npart
		this_part = part_all_zp1[this_index]
		
		this_common = compare_particles(part_z, this_part)
	
		if this_common > min_common:
			this_progenitor = halos_zp1[i_prog]
			progenitors.append(this_progenitor)
			id_progenitors.append(this_index)
			n_progenitors += 1
	
		# TODO add dynamical properties of the halos (relative orientation, masses etc.) to determine the likeliest progenitor
		# compare_dynamics()
	
	# If there is no likely progenitor then we place a token halo instead	FIXME	do a better modeling
	if n_progenitors == 0:
		token_halo = True
		dummy_halo = Halo()
		dummy_halo.x = guess_x
		dummy_halo.m = halo_z.m
		dummy_id = -1
		id_progenitors.append(dummy_id)
		progenitors.append(dummy_halo)

	return (progenitors, id_progenitors, token_halo)


# Old_a is the expansion factor at the time where the halo is in old_x position - FIXME improve expansion implementation 
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


# FIXME Extrapolate the mass of an object at z+1 based on its z value (useful for placeholder-halos)
def mass_z(m0, z):
	a = 0.5
	b = 0.5

	# Check a real formula
	m_z = m0 * pow((1 + z), a) * math.exp(b * z)
	return m_z



# This one compares the dynamical state of two haloes to determine whether it is a progenitor or not
def compare_dynamics(halo_z, halo_zp1):
	# The progenitor halo will be evaluated using these four quantities as proxy
	fac_a = 0.5
	fac_d = 1.0 # Distance from the expected position (in v * dT units)
	fac_r = 1.0 # Distance from the parent halo (in v * dT units)
	fac_m = 2.0
	fac_v = 0.05
	fac_x = 3.0 # Angle between the velocity of the halo and the displacement of the new candidate, to check that no "parallel" object has been 
		    # incorrectly identified
	
	pos = halo_z.x
	vel = halo_z.v
	mass = halo_z.m

	# No likely progenitors 
	if n_zp1 == 0:
		this_eval = 100.0
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

	return eval_this


# Particle IDs are already sorted ! This is important for this algorithm to work
def compare_particles(particle_ids1, particle_ids2):

	n_part1 = len(particle_ids1)
	n_part2 = len(particle_ids2)

	common = 0

	i_part1 = 0
	i_part2 = 0
	
	while (i_part1 < n_part1 and i_part2 < n_part2):
		this_id1 = particle_ids1[i_part1]
		this_id2 = particle_ids2[i_part2]

		if this_id1 == this_id2:
			common += 1
			i_part1 += 1 
			i_part2 += 1
		elif (this_id1 < this_id2):
			i_part1 += 1
		else:
			i_part2 += 1

	return common
