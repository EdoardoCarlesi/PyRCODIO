from .halos import * 
from .utils import *
from .particles import *
import numpy as np 
import sys

#sys.path.append('../libio/')
#from libio.find_files import *
#from libio.read_ascii import *

# Returns an array of all the haloes whose com is located within a given radius from a given halo
def find_halos(halo_c, halo_all, radius):
	n_halos = len(halo_all)
	halo_s = []
	x_c = halo_c.x

	for h in range(0, n_halos):
		halo_0 = halo_all[h]
		d_c = halo_0.distance(x_c)	

		if d_c < radius:
			halo_s.append(halo_0)

	return halo_s

def find_ids(halos, all_ids, all_parts):
	n_halos = len(halos)
	these_ids = []
	these_parts = []

	for i_halo in range(0, n_halos):
		this_id = str(halos[i_halo].ID)
		this_index = all_ids[this_id] 
		these_ids.append(this_id)
		these_parts.append(all_parts[this_index[0]])

	return (these_ids, these_parts)

# Returns an array of all the haloes whose com is located within a given radius from a given point in space, and above a given mass
def find_halos_mass_radius(x_c, halo_all, radius, mass):
	n_halos = len(halo_all)
	halo_s = []

	for h in range(0, n_halos):
		halo_0 = halo_all[h]
		m_c = halo_0.m
		d_c = halo_0.distance(x_c)	
	
		if d_c < radius and m_c > mass:
			halo_s.append(halo_0)

	return halo_s


# Returns an array of all the haloes whose com is located within a given radius from a given point in space
def find_halos_point(x_c, halo_all, radius):
	n_halos = len(halo_all)
	halo_s = []

	for h in range(0, n_halos):
		halo_0 = halo_all[h]
		d_c = halo_0.distance(x_c)	
	
		if d_c < radius:
			halo_s.append(halo_0)

	return halo_s


# Identify a suitable halo pair from a given halo list, a search center, a search radius and a minimum mass/separation for the LG candidates
def find_lg(halos, lgmod):
	center = lgmod.center
	iso_radius = lgmod.d_iso
	m_min = lgmod.m_min
	m_max = lgmod.m_max
	r_min = lgmod.r_min
	r_max = lgmod.r_max
	model_code = lgmod.model_code	
	vrad_max = lgmod.vrad_max

	#print lgmod.info()

	# Center is a three-D variable
	# These are initialized empty
	halos_lg = []		# list of individual halos
	lgs = []		# Haloes paired as LG structures
	halos_center = []	# List of haloes with the right mass range and distance from centrum

	n_halos = len(halos)

	# First identify a set of candidates within the radius and mass range
	for h in range(0, n_halos):
		halo_this = halos[h]

		if halo_this.m > m_min:
			halos_center.append(halo_this)
			
	n_candidates = len(halos_center)
	
	# Total number of LG-halos
	n_all_lgs = 0

	# Now loop on halo center candidates
	for h in range(0, n_candidates):
		halo_lg0 = halos_center[h]
		count_lg = 0
		count_wrong = 0
		
		if halo_lg0.m > m_max:
			count_wrong = 1

		# We need to run the loop on ALL halos, in case there is a third halo lying close
		for i in range(h+1, n_candidates):
			halo_lg1 = halos_center[i]
			dis_this = halo_lg1.distance(halo_lg0.x)

			# This is a subhalo! Keep it				
			if dis_this < halo_lg0.r or dis_this < r_min:
				count_wrong += 0
				#count_wrong += 1
			elif dis_this > r_min and dis_this < r_max:
				halo_lg2 = halo_lg1	# This is a possible candidate
				#vrad = vel_radial(lg1.x, lg2.x, lg1.v, lg2.v)
				vrad = vel_radial(halo_lg2.x, halo_lg0.x, halo_lg2.v, halo_lg0.v) 
				count_lg += 1
				
				# There are too many close-by haloes
				if count_lg > 1:
					count_wrong += 1

				if vrad + (0.67 * dis_this * 0.1) > vrad_max:
					count_wrong += 1

		# A new first & second LG halos have been found:
		if count_wrong == 0 and count_lg == 1:
			already_there = 0

			for j in range(0, n_all_lgs):
				lg1 = halos_lg[2 * j]
				lg2 = halos_lg[2 * j + 1]
			
				if halo_lg0.ID == lg1.ID or halo_lg0.ID == lg2.ID:
					already_there += 1

			# We have a candidate! 
			if already_there == 0:
				n_iso_radius = 0

				# Check also for third haloes within the isolation radius before adding the pair
				com = center_of_mass([halo_lg0.m, halo_lg2.m], [halo_lg0.x, halo_lg2.x])		
				halos_iso = find_halos_point(com, halos, iso_radius)
				nh_iso = len(halos_iso)

				for k in range(0, nh_iso):
					#print halos_iso[k].info()
					if halos_iso[k].m > m_min:
						#print halos_iso[k].distance(com), halos_iso[k].m
						n_iso_radius += 1
	
				if n_iso_radius == 2: 
					# A new first & second LG halos have been found:
					halos_lg.append(halo_lg0)
					halos_lg.append(halo_lg2)
					this_lg = LocalGroup(model_code)
					this_lg.init_halos(halo_lg0, halo_lg2)
					lgs.append(this_lg)
					n_all_lgs += 1

	return lgs


def locate_virgo(ahf_all):
	
	n_ahf = len(ahf_all)
	coord_unit = ahf_all[n_ahf - 1].x[0] + ahf_all[n_ahf - 1].x[1] + ahf_all[n_ahf - 1].x[2] 

	if coord_unit > 10000.:
		virgo_x = [47000.0, 61000.0, 49500.0]
	else:
		virgo_x = [47.0, 61.0, 49.5]

	virgo_r = 5000.
	virgo_m_min = 5.e+13
	virgos = find_halos_point(virgo_x, ahf_all, virgo_r)

	mtotVirgo = 0.0
	m0 = virgo_m_min
	x0 = virgo_x

	for iv in range(0, len(virgos)):
		mtotVirgo += virgos[iv].m
		if virgos[iv].m > m0:
			m0 = virgos[iv].m
			x0 = virgos[iv].x

	print('At position %.3f, %.3f, %.3f found Virgo of mass %.3e. Total mass in a sphere of %.3f kpc/h around it = %.3e ' % \
		(x0[0], x0[1], x0[2], m0, virgo_r, mtotVirgo))

	return (x0, m0, mtotVirgo)


def locate_clusters(ahf_all, box_center, runcode):
	
	n_ahf = len(ahf_all)
	coord_unit = box_center[0] + box_center[1] + box_center[2]

	hubble = 0.677
	#hubble = 1.0
	facMpc = 1000.

	cluster_name = []
	cluster_pos = []
	cluster_dist = []

	cluster_name.append('Virgo             ')
	#cluster_pos.append([-4.67, 16.83, -0.87])
	cluster_pos.append([-4.67 * hubble, 16.83 * hubble, -0.87 * hubble])

	#cluster_name.append('Coma (a)')
	cluster_name.append('Coma               ')
	cluster_pos.append([0.47, 72.55, 10.38])
	#cluster_pos.append([0.47/hubble, 72.55/hubble, 10.38/hubble])

	cluster_name.append('Coma (b)           ')
	cluster_pos.append([-2.43, 68.58, -12.71])

	cluster_name.append('Coma (c)          ')
	cluster_pos.append([-4.27, 74.18, -7.67])

	#cluster_name.append('Leo (a)')
	#cluster_pos.append([-2.41, 71.25, -10.75])

	cluster_name.append('Hercules          ')
	cluster_pos.append([19.28, 57.68, 71.93])

	cluster_name.append('Hercules (a)       ')
	cluster_pos.append([19.28, 57.68, 71.93])

	cluster_name.append('Hercules (b)       ')
	cluster_pos.append([17.49, 63.64, -59.70])

	cluster_name.append('Hercules (c)')
	cluster_pos.append([15.49, 60.94, 74.25])

	cluster_name.append('Perseus             ')
	cluster_pos.append([43.05, -16.89, -21.82])

	cluster_name.append('Perseus-Pisces (a)  ')
	#cluster_name.append('Perseus(a)')
	cluster_pos.append([50.05, -10.89, -12.82])

	cluster_name.append('Perseus-Pisces (b)')
	cluster_pos.append([90.05, -18.39, -15.37])

	cluster_name.append('Perseus-Pisces (c)')
	cluster_pos.append([88.74, -19.15, -19.68])

	#cluster_name.append('Perseus-Pisces (d)')
	#cluster_name.append('Perseus(d)')
	#cluster_pos.append([53.39, -16.06, -5.15])

	#cluster_name.append('Centaurus (a)')
	cluster_name.append('Centaurus           ')
	cluster_pos.append([-34.69, 15.27, -7.77])

	#cluster_name.append('Centaurus (b)')
	#cluster_pos.append([-42.34, 25.11, 2.59])

	cluster_name.append('Hydra')
	cluster_pos.append([-24.67, 21.12, -25.01])

	#cluster_name.append('Pavo-Indus')
	#cluster_pos.append([-22.18, -77.19, 55.09])

	cluster_name.append('Shapley             ')
	cluster_pos.append([-91.88, 34.42, -17.64])

	if coord_unit > 10000.:
		#print box_center	
		dio='porco'
	else:
		print('Rescaling box center units... ')
		for ix in range(0, 3):
			box_center[ix] = box_center[ix] * facMpc
	

	n_cl = len(cluster_name)

	for ic in range(0, n_cl):
		for ix in range(0, 3):
			cluster_pos[ic][ix] = cluster_pos[ic][ix] * hubble * facMpc + box_center[ix]

		#print cluster_pos[ic][:] 

	cluster_r0 = 13000.0 #/ hubble
	cluster_m = 0.5e+14
	
	ahf_x = []
	ahf_m = []

	for ic in range(0, n_cl):
		m0 = cluster_m
		x0 = [-1.0, -1.0, -1.0]
		this_x = cluster_pos[ic]

		d_center  = distance(box_center, this_x)
		cluster_r = cluster_r0 #* (1. + 0.1 * (d_center / (2 * cluster_r0)))
		clusters  = find_halos_point(this_x, ahf_all, cluster_r)

		#print 'Found: ', len(clusters), ' around: ', this_x, cluster_name[ic]

		for iv in range(0, len(clusters)):

#			print clusters[iv].info()

			if clusters[iv].m > m0:
				#m0 = clusters[iv].m
				x0 = clusters[iv].x
				d0 = distance(this_x, x0)
				#print cluster_name[ic], ' ',  clusters[iv].m, clusters[iv].x, d_center/1000., d0/1000.
				#print cluster_name[ic], ' ',  clusters[iv].m/1.e+14, \
				print(str(runcode), '\t', cluster_name[ic], '\t\t', \
                                clusters[iv].m/1.e+14/hubble, '\t', d0/facMpc, '\t\t', \
				(clusters[iv].x[0]-box_center[0])/facMpc,'\t',\
                                (clusters[iv].x[1]-box_center[1])/facMpc,'\t',\
                                (clusters[iv].x[2]-box_center[2])/facMpc) #, \
				#clusters[iv].v[0], clusters[iv].v[1], clusters[iv].v[2],  \
				#clusters[iv].l[0], clusters[iv].l[1], clusters[iv].l[2])#, '\n'
	
		cluster_dist.append(distance(cluster_pos[ic], x0))

		#print cluster_name[ic], cluster_dist[ic], m0#, x0[0] - box_center[0], x0[1] - box_center[1], x0[2] - box_center[2]

		ahf_m.append(m0)
		ahf_x.append(x0)

	return (ahf_x, ahf_m, cluster_name)



def print_subhalos(com, mcut, lg_halos, fname_lg):
    #f_out = open(fname_lg, 'wb')
    f_out = open(fname_lg, 'w')

    for halo in lg_halos:
        if halo.m > mcut:
            this_d = halo.distance(com)

            if this_d > 10.0:
                f_out.write('%s\n' % (halo.line.rstrip('\n')))
            elif this_d == 0.0:
                f_out.write('# Host ID %s\n' % (halo.ID))

    f_out.close()

def print_subhalo_ids(com, mcut, lg_halos, fname_lg):
    f_out = open(fname_lg, 'w')

    for halo in lg_halos:
        if halo.m > mcut:
            this_d = halo.distance(com)
            f_out.write('%s\n' % (halo.ID))

    f_out.close()


'''
# Find the best match for each halo, ahf_now contains all the haloes we want to trace at z and ahf_back all the possible candidates at z+1
def match_progenitors(ahf_now, ahf_back):
	n_now = len(ahf_now)
	n_back = len(ahf_back)

# FIXME check this function maybe it is not needed to track all the subhalos... maybe it is... restructure this anyway FIXME
def halos_and_subhalos_through_z(end_snap, ini_snap, base_path, root_file, suff_halo, suff_part, main_halos, main_parts, r_subs):
	zs = ahf_redshifts(base_path)
	ss = ahf_snapshots(base_path)

	n_halos = len(halos)

	all_halo_z = []
	old_halo = []
	old_part = []
	all_sub_z = [] #* n_halos
	old_sub = []
	old_sub_part = []

#	print old_sub

	for ihz in range(0, n_halos):
		halo_z = halo.HaloThroughZ(end_snap - ini_snap)
		all_halo_z.append(halo_z)
		dummy_list = []
		all_sub_z.append(dummy_list)
		old_sub.append(dummy_list)


	for i_snap in range(end_snap, ini_snap, -1):
		# TODO this assumes one MPI task only !!!!!
		this_part = base_path + root_file + ss[i_snap] + '.0000.z' + zs[i_snap] + suff_part
		this_halo = base_path + root_file + ss[i_snap] + '.0000.z' + zs[i_snap] + suff_halo

		this_halos = read_ahf(this_halo)
		this_parts = read_particles(this_part)

		this_z = float(zs[i_snap])

		this_t = z2Myr(this_z)
		this_a = 1.0 / (1.0 + this_z)

		# This is the first step
		if i_snap == end_snap:

			# Loop on the main halos to track (usually the two LG main members)
			for i_main in range(0, n_halos):
				this_halo = main_halos[i_main]
				this_part = main_parts[i_main]
				all_halo_z[i_main].add_step(this_halo, this_t, this_z)
				
				# Trace all the halos a factor r_sub away from the main halo
				rsub = this_halo.r * r_subs
				this_halo.sub_halos(this_halos, rsub)
				subs = this_halo.sort_sub_mass()
				n_sub = len(subs)

				# Save all the halos for the next step
				old_halo.append(this_halo)
				old_part.append(this_part)
				tmp_subs = []
				tmp_part_subs = []

				# Now track all the subhalos
				for i_sub in range(0, n_sub):
					this_sub = subs[i_sub]
					this_sub_part = subs[i_sub]
					old_sub[i_main].append(this_sub)
					tmp_sub = halo.SubHaloThroughZ(end_snap - ini_snap)
					tmp_sub.add_step(this_sub, this_t, this_z)
					tmp_sub.host.add_step(this_halo, this_t, this_z)
					tmp_subs.append(tmp_sub)

				# Allocate subhalo list for the i_main-th halo
				all_sub_z[i_main] = [Halo.SubHaloThroughZ(end_snap-ini_snap) for i in range(0, n_sub)]

				for i_sub in range(0, n_sub):
					all_sub_z[i_main][i_sub] = tmp_subs[i_sub]

			old_t = this_t

		# Normal steps after the first one
		else:	

			timeStep = abs(old_t - this_t)		
			this_part = []	# this_parts contains ALL the part files, this_part only those of the progenitors

			# FIXME the find progenitor routine has changed now!!!!!
			for i_main in range(0, n_halos):
				# This function returns more than one possible progenitor, sorted by merit
				this_halo = th.find_progenitors(old_halo[i_main], old_part[i_main], this_halos, this_parts, min_common, this_a, timeStep)	
				for i_halo in range(0, len(this_halos)):
					this_id = find_id(this_halo[i_halo].ID, this_halos)
					this_part.append(this_parts[this_id])
					
				old_halo[i_main] = this_halo[0]
				all_halo_z[i_main].add_step(this_halo, this_t, this_z)
				
				for i_sub in range(0, n_sub):
					try:
						this_sub = th.find_progenitor(old_sub[i_main][i_sub], this_halos, this_a, timeStep)	
						doSubs = True
					except:
			# print 'Subhalo[%d][%d] at step z=%.3f has no likely progenitor.' % (i_main, i_sub, this_z)
						doSubs = False

					if doSubs:
						old_sub[i_main][i_sub] = this_sub
						all_sub_z[i_main][i_sub].add_step(this_sub, this_t, this_z)
						all_sub_z[i_main][i_sub].host.add_step(this_halo, this_t, this_z)

			old_t = this_t

	return (all_halo_z, all_sub_z)
'''


