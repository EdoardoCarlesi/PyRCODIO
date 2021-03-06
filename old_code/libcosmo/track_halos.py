from .find_halo import *
from .utils import *
from .units import *
from .halos import *

import pickle
import numpy as np
import math



def merger_tree_chunks(end_snap, ini_snap, n_chunks, n_use, min_common, main_halos, main_parts, main_ids, settings):
	zs = settings.redshifts

	n_halos = len(main_halos)

	# Final list containing HaloThroughZ objects
	all_halo_z = [] 

	# Append a list of subhalos to each main id halo being tracked
	#for i_sub in range(0, n_ids):
	#	this_subs = []
	#	all_subs_z[i_sub].append(this_subs)

	# Temporary lists that allow intra-snapshot comparisons - only for the haloes being tracked
	old_main_halo = []
	old_main_part = []
	old_main_ids = []
	
	# Initialize halo_z list members 
	for i_halo in range(0, n_halos):
		halo_z = HaloThroughZ(end_snap - ini_snap)
		all_halo_z.append(halo_z)

	# Now loop over all the particle and halo files
	for i_snap in range(end_snap, ini_snap, -1):
		this_task = 0				# FIXME this assumes one mpi task only
		(this_root, suff_halo, suff_part) = settings.get_chunk_files(i_snap, n_chunks)

		# Time steps, redshift and expansion factors
		this_z = float(settings.redshifts[i_snap * n_chunks])
		this_t = z2Myr(this_z)
		this_a = 1.0 / (1.0 + this_z)

		print(i_snap, n_chunks, this_root, suff_halo, suff_part, this_t, this_z)

		# If this is the last snap just use all the inf
		if i_snap == end_snap:

			# Loop on the main halos to be tracked (usually the two LG main members)
			for i_main in range(0, n_halos):
				this_halo = main_halos[i_main]
				this_part = main_parts[i_main]
				this_ids  = main_ids[i_main]

				all_halo_z[i_main].add_step(this_halo, this_z)

				# Save all the halos for the next step
				old_main_halo.append(this_halo)
				old_main_part.append(this_part)
				old_main_ids.append(this_ids)

			# Do this once the loop on all the halos is finished
			old_t = this_t
	
		# This is not the first step - 
		else:
			#print 'Reading files: ', this_part_file, this_halo_
			this_all_halo = read_ahf_chunks(this_root, suff_halo, n_use)
			n_this_halos = len(this_all_halo)
			(this_all_ids, this_all_parts) = read_particles_chunks(this_root, suff_part, n_use, n_this_halos)

			# Time step is deltaT, it might be a fixed quantity but we don't know in advance so we compute it every time
			time_step = abs(old_t - this_t)	

			for i_main in range(0, n_halos):
				#print 'oldmain: ', i_main, old_main_halo[i_main].info()
				(this_halo, token_halo) = find_progenitors(old_main_halo[i_main], old_main_part[i_main],\
						this_all_halo, this_all_parts, min_common, this_a, time_step)

				this_id = this_halo.ID

				try:
					this_index = this_all_halo[0].id_index[str(this_id)]
				except:
					print('Halo ID %s not found.' % this_id)

				# If no likely progenitor has been found but only a token halo has been returned, then save the old particle list
				# also for the next step
			if token_halo == True:
				this_part = old_main_part[i_main]
				this_ids = old_main_ids[i_main]
			else:	
				# Find the particle ID list of the main progenitor
				this_part = this_all_parts[this_index] 
				this_ids = this_id
					
			# Now save the particles and halo for the next step
			old_main_halo[i_main] = this_halo
			old_main_part[i_main] = this_part
			old_main_ids[i_main] = this_ids
			
			# Store the halo into the general structure
			all_halo_z[i_main].add_step(this_halo, this_z)

		old_t = this_t

	return all_halo_z



'''
	Computes the merger tree for a given number of halos. This is generally the two main LG halos at z=0 and their satellites
	Subhaloes are stored separately if the trace_subs boolean variable is set to true. Then ids_subs contains the ids of the main haloes
	whose subhaloes need to be tracked at each step, using a factor r_subs * Rvir to identify them
'''
def merger_tree(end_snap, ini_snap, min_common, main_halos, main_parts, main_ids, all_halos, settings, track_subs, ids_subs, r_subs):
	zs = settings.redshifts

	n_halos = len(main_halos)

	# Final list containing HaloThroughZ objects
	all_halo_z = [] 

	# How many haloes among the _main ones do we want to track with subhalo distribution?
	n_ids = len(ids_subs)

	# These are simply all the haloes within a given radius at each step
	all_subs_z = [[] for i_sub in range(0, n_ids)]
	track_sub_index = []

	# Append a list of subhalos to each main id halo being tracked
	#for i_sub in range(0, n_ids):
	#	this_subs = []
	#	all_subs_z[i_sub].append(this_subs)

	# Temporary lists that allow intra-snapshot comparisons - only for the haloes being tracked
	old_main_halo = []
	old_main_part = []
	old_main_ids = []
	
	# Initialize halo_z list members 
	for i_halo in range(0, n_halos):
		halo_z = HaloThroughZ(end_snap - ini_snap)
		all_halo_z.append(halo_z)

	# Now loop over all the particle and halo files
	#for i_snap in range(ini_snap, end_snap):
	for i_snap in range(end_snap, ini_snap, -1):
		this_task = 0				# FIXME this assumes one mpi task only
		(this_part_file, this_halo_file) = settings.get_ahf_files(i_snap, this_task)
		
		# Time steps, redshift and expansion factors
		this_z = float(settings.redshifts[i_snap])
		this_t = z2Myr(this_z)
		this_a = 1.0 / (1.0 + this_z)
	
		#print this_t, this_z

		# This index is reinitialized to 0 at every loop
		i_subs = 0

		# If this is the last snap just use all the inf
		if i_snap == end_snap:
		#if i_snap == ini_snap:

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

				#print i_main, this_halo.info()
			
				if track_subs == True:
					(there_is, entries) = is_there(this_ids, ids_subs)

					if there_is == True:
						#print 'Tracking subs, ', there_is, entries, ids_subs, this_ids
						this_r = this_halo.r * r_subs
						sub_halos = find_halos(this_halo, all_halos, this_r)
						this_subs = SubHalos(this_halo, sub_halos)
						all_subs_z[i_subs].append(this_subs)
						track_sub_index.append(i_main)
						i_subs += 1

			# Do this once the loop on all the halos is finished
			old_t = this_t

		# This is not the first step - 
		else:
			print('Reading files: ', this_part_file, this_halo_file)
			this_all_halo = read_ahf(this_halo_file)
			(this_all_ids, this_all_parts) = read_particles(this_part_file)

			#this_t = z2Myr(this_z)

			# Time step is deltaT, it might be a fixed quantity but we don't know in advance so we compute it every time
			time_step = abs(old_t - this_t)	

			for i_main in range(0, n_halos):
				#print 'oldmain: ', i_main, old_main_halo[i_main].info()

				(this_halo, token_halo) = find_progenitors(old_main_halo[i_main], old_main_part[i_main],\
							this_all_halo, this_all_parts, min_common, this_a, time_step)

				this_id = this_halo.ID

				try:
					this_index = this_all_halo[0].id_index[str(this_id)]
				except:	
					this_ids = this_id
					#print 'Halo ID %s not found.' % this_id

				# If no likely progenitor has been found but only a token halo has been returned, then save the old particle list
				# also for the next step
				if token_halo == True:
					this_part = old_main_part[i_main]
					this_ids = old_main_ids[i_main]
				else:	
					# Find the particle ID list of the main progenitor
					this_part = this_all_parts[this_index] 
					
				# Now save the particles and halo for the next step
				old_main_halo[i_main] = this_halo
				old_main_part[i_main] = this_part
				old_main_ids[i_main] = this_ids
			
				# Store the halo into the general structure
				all_halo_z[i_main].add_step(this_halo, this_z)
					
			if track_subs == True:
				for i_sub in range(0, n_ids):
					#print i_sub, track_sub_index
					this_halo = old_main_halo[track_sub_index[i_sub]]
					this_r = this_halo.r * r_subs
					sub_halos = find_halos(this_halo, this_all_halo, this_r)
					this_subs = SubHalos(this_halo, sub_halos)
					#print i_main, len(sub_halos), this_halo.info()
					all_subs_z[i_subs].append(this_subs)
					i_subs += 1

			old_t = this_t

	return (all_halo_z, all_subs_z)



def find_halo_id(idnum, ahf_halos):
	n_halo = len(ahf_halos)
	i_halo = 0
	h_found = False	
	
	id_str = str(idnum)

	try:
		this_id = ahf_halos[0].id_index[id_str]
		#print 'ID: ', this_id
		return this_id
	except:
		#print 'No ID: %ld.' % idnum
		return -1



def find_halo_id_old(idnum, ahf_halos):
	n_halo = len(ahf_halos)
	i_halo = 0
	h_found = False	
	
	id_str = str(idnum)

	while (i_halo < n_halo) and (h_found == False):
		if (ahf_halos[i_halo].ID == idnum):
			h_found = True	
		else:
			i_halo += 1

	if h_found == True:
		return i_halo
	else:
		print('Halo ID: %ld not found.' % idnum)
		print(ahf_halos[i_halo].info())
		return -1



# This assumes that halos_zp1 already contains all the closest halos to halo_z, the candidate for which we are looking for a progenitor
def find_progenitors(halo_z, part_z, halos_all_zp1, part_all_zp1, min_common, aFactor, timeStep):
	r_prog = 0.0
	guess_x = backward_x(halo_z.x, halo_z.v, aFactor, timeStep)
	r_prog = 3.0 * module(halo_z.v) * timeStep * s2Myr() / km2kpc() / aFactor	# In comoving units

	halos_zp1 = find_halos_point(guess_x, halos_all_zp1, r_prog)
	
	token_halo = False
	progenitors = []
	scores = []
	n_progenitors = 0	
	score_zero = 25.
	npart = halo_z.npart
	n_zp1 = len(halos_zp1)
	id_z = halo_z.ID
	m_z = halo_z.m

	#print '----> Find progenitors: ', guess_x, halo_z.x, r_prog, ' possible progenitors: ', n_zp1

	for i_prog in range(0, n_zp1):
		this_m = halos_zp1[i_prog].m
		this_id = halos_zp1[i_prog].ID
		#this_index = find_halo_id_old(this_id, halos_all_zp1)
		this_index = find_halo_id(this_id, halos_all_zp1)
		this_npart = halos_all_zp1[this_index].npart
		this_part = part_all_zp1[this_index]

		this_common = compare_particles(part_z, this_part)
		#print i_prog, this_common

		if this_common > min_common:
			this_progenitor = halos_zp1[i_prog]
			this_npart = this_progenitor.npart
			progenitors.append(this_progenitor)

			share_son = float(this_common) / float(this_npart )
			share_dad = float(this_common) / float(npart)

			this_score = compare_dynamics(halo_z, this_progenitor, aFactor, timeStep)
			this_score /= (share_son * share_dad)
			scores.append(this_score)	
			
			if this_score < score_zero:
			#	print 'Progenitor found'
			#	print 'OriginalID:', this_id, ' step: ', i_prog, ' common: ', this_common, this_score, share_son, share_dad
				return_halo = this_progenitor
				#index_prog = i_prog
				score_zero = this_score
				this_share_son = share_son
				this_share_dad = share_dad
				return_common = this_common
				n_progenitors += 1
			#else:
				#print '(wrong)IDs: %ld, %ld | Ms: %3.3e, %3.3e | common: %d | %7.3f, %7.3f, %7.3f' % \
				#	(id_z, this_id, m_z, this_m, this_common, this_score, share_son, share_dad)
			#	print id_z, this_id, m_z, this_m, ' common: ', this_common, this_score, share_son, share_dad
				#print halo_z.info()
				#print this_progenitor.info()
				#print '...'
	#if n_progenitors > 0:
		#print 'IDs: %ld, %ld | Ms: %3.3e, %3.3e | common: %d | %7.3f, %7.3f, %7.3f' % \
		#	(id_z, return_halo.ID, m_z, return_halo.m, return_common, score_zero, this_share_son, this_share_dad)
	
	# If there is no likely progenitor then we place a token halo instead	FIXME	do a better modeling
	if n_progenitors == 0:
		#print 'Progenitor for %ld (npart=%d) not found, replacing with token halo.' % (id_z, halo_z.npart)
		token_halo = True
		return_halo = compute_token_halo(aFactor, timeStep, guess_x, halo_z)

	elif n_zp1 == 0:
		#print 'No likely progenitors for %ld.' 
		token_halo = True
		return_halo = compute_token_halo(aFactor, timeStep, guess_x, halo_z)

	return (return_halo, token_halo)


def compute_token_halo(a0, tStep, guess_x, halo):
	dummy_halo = Halo()
	dummy_halo.ID = -1	

	if tStep < 100.:
		tStep *= 1000.

	#print 'Compute token halo ', a0, tStep, guess_x
	t1 = a2Myr(a0) - tStep
	a1 = Myr2a(t1)
	z1 = 1.0/ a0 - 1.0
	z2 = 1.0/ a1 - 1.0

#	print a0, a1

	m0 = halo.m / halo_mass_z(z1)
	m_z2 = m0 * halo_mass_z(z2)

	dummy_halo.m = m_z2	
	dummy_halo.x = halo.x	#FIXME write a meaningful guess_x function
	
	frac_m = m_z2 / halo.m
	npart_new = int(halo.npart * frac_m)
	dummy_halo.npart = npart_new
	dummy_halo.v[0] = 0.0
	dummy_halo.v[1] = 0.0
	dummy_halo.v[2] = 0.0

	#print 'Zold: ', z1, '         Znew: ', z2
	#print 'Mold: ', halo.m, '     Mnew: ', m_z2 
	#print 'Nold: ', halo.npart, ' Nnew: ', npart_new

	return dummy_halo

def halo_mass_z(z):
	alpha = 6.1
	beta = 2.7
	
	fac0 = np.power((1 + z), beta)
	fac1 = np.exp(-alpha * (np.sqrt(1 + z) - 1))

	return fac0 * fac1


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
		#new_x[ix] = this_x[ix] * facMpc - (this_vel[ix] * tot_t) / km2kpc() * 0.677
		# TODO correct for the hubble expansion
		# TODO correct for Newton acceleration

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
def compare_dynamics(halo_z, halo_zp1, aFactor, timeStep):
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

	#TODO add some check if ID = -1 OR velocity is 0 then do not care about the dynamics just look at the particles

	# FIXME
	guess_x = backward_x(pos, vel, aFactor, timeStep)

	halo_x = halo_zp1.x
	halo_v = halo_zp1.v
	halo_m = halo_zp1.m
	r_prog = halo_zp1.r
			
	direction = vec_subt(pos, halo_x) 
			
	this_x = 1.0 - abs_val(angle(direction, halo_v))
	this_a = 1.0 - abs_val(angle(halo_v, vel))
	this_d = distance(guess_x, halo_x) / r_prog #distance(pos, this_x))
	this_r = distance(pos, halo_x) / r_prog #distance(pos, this_x))
	this_m = abs_val(np.log10(mass/halo_m))
	this_v = abs_val(distance(halo_v, vel)/module(vel))

	#eval_this = fac_a * this_a + fac_m * this_m + fac_d * this_d + fac_v * this_v + fac_r * this_r + this_x * fac_x
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
