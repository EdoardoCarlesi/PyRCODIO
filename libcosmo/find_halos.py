from halo import *
from utils import *
from find_halos import *
from particles import *
import numpy as np 

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
	radius = lgmod.d_max
	iso_radius = lgmod.d_iso
	m_min = lgmod.m_min
	r_min = lgmod.r_min
	r_max = lgmod.r_max
	
	print lgmod.info()

	# Center is a three-d variable
	# These are initialized empty
	halos_lg = []
	halos_center = []	# List of haloes with the right mass range and distance from centrum

	n_halos = len(halos)

	# First identify a set of candidates within the radius and mass range
	for h in range(0, n_halos):
		halo_this = halos[h]

		# print halo_this.distance(center)
		if halo_this.m > m_min and halo_this.distance(center) < radius+r_max :
			halos_center.append(halo_this)
			
	n_candidates = len(halos_center)
	
	# Total number of LG-halos
	n_all_lgs = 0

	# Now loop on halo center candidates
	for h in range(0, n_candidates):
		halo_lg0 = halos_center[h]
		count_lg = 0
		count_wrong = 0

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
				count_lg += 1
				
				# There are too many close-by haloes
				if count_lg > 1:
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
				#print halo_lg0.info()
				# Check also for third haloes within the isolation radius before adding the pair
				com = center_of_mass([halo_lg0.m, halo_lg2.m], [halo_lg0.x, halo_lg2.x])		
				halos_iso = find_halos_point(com, halos, iso_radius)
				nh_iso = len(halos_iso)
				#print com
				#print iso_radius

				#print 'NH iso: %d' % nh_iso
				#print halos[0].info()
				#print halos_iso[0].info()

				for k in range(0, nh_iso):
					#print halos_iso[k].info()
					if halos_iso[k].m > m_min:
							#print halos_iso[k].distance(com), halos_iso[k].m
						n_iso_radius += 1
	
				if n_iso_radius == 2: 
					# A new first & second LG halos have been found:
					halos_lg.append(halo_lg0)
					halos_lg.append(halo_lg2)
					n_all_lgs += 1

	return halos_lg


def rate_lg_pair(lg1, lg2, box_center):
	# Benchmark (obs.) quantities
	dcbox0 = 5000.	# kpc/h
	rhalo0 = 500.	# kpc/h
	vrad0 = -100.
	mass0 = 3.0e+12
	ratio0 = 1.1
	hubble0 = 67.0

	box = 100.0
	npart = 512

	com = center_of_mass([lg1.m, lg2.m], [lg1.x, lg2.x])
	c12 = vec_subt(box_center, com)
	m12 = lg1.m + lg2.m

	if lg1.m > lg2.m:
		rm12 = lg1.m / lg2.m
	else:
		rm12 = lg2.m / lg1.m

	dcenter = vec_module(c12)
	rhalos = lg1.distance(lg2.x)
	vrad = vel_radial(lg1.x, lg2.x, lg1.v, lg2.v)
	vrad += hubble0 * rhalos/1000.

	# Relative weights to estimate the relevance of each property relative to the "real" LG 
	fac_rh = 1.0
	fac_c = 0.25
	fac_v = 0.25
	fac_m = 1.5
	fac_ra = 1.5

	# How to compute the LG-likeliness factors
	diff_c = abs(dcenter) / dcbox0
	diff_rh = abs(rhalos - rhalo0) / rhalo0
	diff_m = np.log10(m12 / mass0)
	diff_v = abs(vrad0 - vrad) / abs(vrad0)
	diff_ra = abs(rm12 - ratio0) / abs(ratio0)

#	print "LG1: %s" % lg1.info()
#	print "LG2: %s" % lg2.info()

	lg_rate = diff_rh * fac_rh + fac_c * diff_c + diff_m * fac_m + diff_ra * fac_ra + fac_v * diff_v
	
	# Get a penalty for positive vrad
	if vrad > 0.0:
		lg_rate += 10.

	contamin = abs_val((lg1.m/lg1.npart) - simu_pmass(box, npart))/simu_pmass(box, npart)
	
	print 'LG rating: %.3f, Npart: %d & %d,  Res.Factor: %.3f \n' % (lg_rate, lg1.npart, lg2.npart, contamin)
#	print 'Values.   RH: %f, C:%f, M:%e, V:%f, RA:%f' % (rhalos, dcenter, m12, vrad, rm12)
#	print 'Diff: %f, rh: %f, c:%f, m:%f, v:%f, ra:%f' % (lg_rate, diff_rh, diff_c, diff_m, diff_v, diff_ra)

	return lg_rate


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

	print 'At position %.3f, %.3f, %.3f found Virgo of mass %.3e. Total mass in a sphere of %.3f kpc/h around it = %.3e ' % \
		(x0[0], x0[1], x0[2], m0, virgo_r, mtotVirgo)

	return (x0, m0, mtotVirgo)


def locate_clusters(ahf_all, box_center):
	
	n_ahf = len(ahf_all)
	coord_unit = box_center[0] + box_center[1] + box_center[2]

	hubble = 0.677
	facMpc = 1000.

	cluster_name = []
	cluster_pos = []

	cluster_name.append('Virgo')
	cluster_pos.append([-4.67, 16.83, -0.87])

	cluster_name.append('Coma (a)')
	cluster_pos.append([0.47, 72.55, 10.38])

	cluster_name.append('Coma (b)')
	cluster_pos.append([-2.43, 68.58, -12.71])

	cluster_name.append('Coma (c)')
	cluster_pos.append([-4.27, 74.18, -7.67])

	cluster_name.append('Leo (a)')
	cluster_pos.append([-2.41, 71.25, -10.75])

	cluster_name.append('Hercules (a)')
	cluster_pos.append([19.28, 57.68, 71.93])

	cluster_name.append('Hercules (b)')
	cluster_pos.append([17.49, 63.64, -59.70])

	cluster_name.append('Hercules (c)')
	cluster_pos.append([15.49, 60.94, 74.25])

	cluster_name.append('Perseus-Pisces (a)')
	cluster_pos.append([50.05, -10.89, -12.82])

	cluster_name.append('Perseus-Pisces (b)')
	cluster_pos.append([90.05, -18.39, -15.37])

	cluster_name.append('Perseus-Pisces (c)')
	cluster_pos.append([88.74, -19.15, -19.68])

	cluster_name.append('Perseus-Pisces (d)')
	cluster_pos.append([53.39, -16.06, -5.15])

	cluster_name.append('Centaurus (a)')
	cluster_pos.append([-34.69, 15.27, -7.77])

	cluster_name.append('Centaurus (b)')
	cluster_pos.append([-42.34, 25.11, 2.59])

	cluster_name.append('Hydra')
	cluster_pos.append([-24.67, 21.12, -25.01])

	cluster_name.append('Pavo-Indus')
	cluster_pos.append([-22.18, -77.19, 55.09])

	cluster_name.append('Shapley (a)')
	cluster_pos.append([-39.98, 19.87, 3.63])

	cluster_name.append('Shapley (b)')
	cluster_pos.append([-91.88, 34.42, -17.64])

	if coord_unit > 10000.:
		n_cl = len(cluster_name)

		for ic in range(0, n_cl):
			for ix in range(0, 3):
				cluster_pos[ic][ix] = cluster_pos[ic][ix] * hubble * facMpc + box_center[ix]

	cluster_r = 12000.0
	cluster_m = 0.7e+13
	
	ahf_x = []
	ahf_m = []

	for ic in range(0, n_cl):
		m0 = cluster_m
		x0 = [-1.0, -1.0, -1.0]
		this_x = cluster_pos[ic]
		clusters = find_halos_point(this_x, ahf_all, cluster_r)

		for iv in range(0, len(clusters)):

			if clusters[iv].m > m0:
				m0 = clusters[iv].m
				x0 = clusters[iv].x
				
		#print cluster_name[ic], x0, m0

		ahf_m.append(m0)
		ahf_x.append(x0)

	return (ahf_x, ahf_m, cluster_name)



def find_halo_id(idnum, ahf_halos):

	n_halos = len(ahf_halos)
	i_halos = 0
	h_found = False	

	while (i_halos < n_halos) and (h_found == False):
		if (ahf[i_halo].ID == idnum):
			h_found = True	
		else:
			i_halos += 1

	if h_found == True:
		return i_halos
	else:
		print 'Halo ID: %ld not found.' % idnum
		return -1


# Find the best match for each halo, ahf_now contains all the haloes we want to trace at z and ahf_back all the possible candidates at z+1
def match_progenitors(ahf_now, ahf_back):
	n_now = len(ahf_now)
	n_back = len(ahf_back)

