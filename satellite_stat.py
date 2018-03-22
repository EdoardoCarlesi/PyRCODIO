#!/usr/bin/python

import numpy as np
import os

from libcosmo.utils import *
from libcosmo.halo import Halo
from libcosmo.read_ahf import read_ahf
from libcosmo.find_halos import *

resolution='2048'
#resolution='2048b'
#resolution='1024'
run_init = 4
run_end = 5

subrun_init = 0
subrun_end = 10

base_path = '/home/eduardo/CLUES/DATA/'

all_runs = ['00_06', '01_11', '01_13', '34_04', '55_02']
snapshot = 'snapshot_054.0000.z0.000.AHF_halos'

hubble = (67.1 / 1000.)		# this is going to multiply kpc/h distances - thus it is rescaled by a factor of 1000.
#file_lg_header = '# ID_M31(1) ID_MW(2) M_M31(3)[Msun] M_MW(4)[Msun] R_MwM31(5)[Mpc] Vrad_M31(6)[phys, km/s] Nsub_M31(7) Msub_MaxM31(8) Nsub_MW(9) Msub_MaxMW (10) SimuCode(11)\n'
file_lg_header = '# ID_M31(1) ID_MW(2) M_M31(3)[Msun] M_MW(4)[Msun] R_MwM31(5)[Mpc] Vrad_M31(6)[phys, km/s] Nsub_M31(7) Nsub_MW(8) SimuCode(9)\n'
file_sub_header = '# ID_sub(1) 	M_sub(2)[Msun] 	R_sub(3)[Mpc]	Host(4) 	SimuCode(5)\n' 

file_lg_name = 'lg_candidates_'+resolution+'.txt'
file_sub_name = 'lg_subhalos_'+resolution+'.txt'
base_path = '/home/eduardo/CLUES/DATA/' 

# LG Candidates properties .txt file
file_lg_txt = open(base_path + file_lg_name, 'wb')
file_lg_txt.write(file_lg_header)
file_sub_txt = open(base_path + file_sub_name, 'wb')
file_sub_txt.write(file_sub_header)

# Local group selection parameters
center = [50000., 50000., 50000.]
radius = 7500.
iso_radius = 2000.
r_max = 1100.
r_min = 250.
m_min = 4.5e+11
m_max = 10.e+12

r_sub_min = 10.0   # Minimum distance from the main halo - if it's too close it might be double counting the host
m_sub_min = 7.e+9 # This mass is in REAL Msun units, without the /h

vrad_max = -5.0
n_sub_min = 4 # Minimum number of subhalos required to compute the moment of inertia

# Subhalo identification criterion
fac_r = 1.5
np_sub_min = 30

# LGs will be appendend to this as they are found
all_lgs = []

# when several LG-like pairs are found, get the first pair (0) second pair (2) etc.
ind_lg = 0

for subrun_i in range(subrun_init, subrun_end):
	run_num = '%02d' % subrun_i
	
	for run_j in range(run_init, run_end):
		base_run = all_runs[run_j]
		this_file = base_path + resolution + '/' + base_run + '/' + run_num + '/' + snapshot
		print_run = base_run + '_' + run_num
		
		if os.path.exists(this_file):
				print 'Reading in AHF file: ', this_file
				ahf_all = read_ahf(this_file)
				these_lg = find_lg(ahf_all, center, radius, iso_radius, m_min, r_min, r_max)
		else:
				n_lgs = 0;
		
		n_lgs = int(len(these_lg) / 2)
		print 'Found a total of %d LG pairs.' % (n_lgs)

		if n_lgs > 0:
			print these_lg[ind_lg].info()
			print these_lg[ind_lg+1].info()
			all_lgs.append(these_lg[ind_lg])
			all_lgs.append(these_lg[ind_lg+1])

			this_lg1 = these_lg[ind_lg]
			this_lg2 = these_lg[ind_lg+1]

			these_sub1 = find_halos(this_lg1, ahf_all, fac_r * this_lg1.r)
			these_sub2 = find_halos(this_lg2, ahf_all, fac_r * this_lg2.r)

			this_lg1 = these_lg[ind_lg]
			this_lg2 = these_lg[ind_lg+1]
			count_sub1 = 0
			count_sub2 = 0

			masses = []
		
			# Pair properties
			lg_com = center_of_mass(this_lg1.m, this_lg2.m, this_lg1.x, this_lg2.x)
			d12 = this_lg1.distance(this_lg2.x) 
			m12 = this_lg1.m + this_lg2.m
			vel_r = vel_radial(this_lg1.x, this_lg2.x, this_lg1.v, this_lg2.v)
			print vel_r, d12, hubble
			vel_r += hubble * d12
			print vel_r
			h0 = hubble * 10.0
			file_lg_line = '%ld   %ld   %.2e   %.2e   %.2f   %.2f   %d   %d   %s\n' % \
			(this_lg1.ID, this_lg2.ID, this_lg1.m/h0, this_lg2.m/h0, d12/h0, vel_r, \
			 this_lg1.nsub, this_lg2.nsub, print_run)

			if vel_r < vrad_max and m12 < m_max:
				file_lg_txt.write(file_lg_line)
				print file_lg_line			

			coords = [[0] * 3]
		
			for sub_i in range(0, this_lg1.nsub):
				this_sub1 = these_sub1[sub_i]

				if this_sub1.npart > np_sub_min:
					count_sub1 += 1		
			
					masses.append(this_sub1.m)
					new_coords = [0] * 3
					center = this_lg1.x
					new_coords[0] = this_sub1.x[0] - center[0]
					new_coords[1] = this_sub1.x[1] - center[1]
					new_coords[2] = this_sub1.x[2] - center[2]
					coords.append(new_coords)
				
				this_r = this_sub1.distance(this_lg1.x)
				
				#print 'id: %d, np= %d, m= %e, d = %.2f' % (sub_i, this_sub1.npart, this_sub1.m, this_sub1.distance(this_lg1.x))
				if this_sub1.m > m_sub_min * h0 and this_r > r_sub_min:
					file_sub_line = '%ld	%.2e	%.2f	%s	%s\n' % \
					(this_sub1.ID, 	this_sub1.m / h0, this_r / h0, 'M31', print_run)
					file_sub_txt.write(file_sub_line)

			if this_lg1.nsub > n_sub_min:
				moment_inertia(coords, masses)

			#print masses
			#print coords
			print 'LG1 %d) has %d subs over %d particles\n' % (subrun_i, count_sub1, np_sub_min)
	
			for sub_i in range(0, this_lg2.nsub):
				this_sub2 = these_sub2[sub_i]
				if this_sub2.npart > np_sub_min:
					count_sub2 += 1		

					masses.append(this_sub2.m)
					#masses.append(1.0)
					new_coords = [0] * 3
					center = this_lg2.x
					new_coords[0] = this_sub2.x[0] - center[0]
					new_coords[1] = this_sub2.x[1] - center[1]
					new_coords[2] = this_sub2.x[2] - center[2]
					coords.append(new_coords)
	
				this_r = this_sub2.distance(this_lg2.x)

				if this_sub2.m > m_sub_min * h0 and this_r > r_sub_min:
					file_sub_line = '%ld	%.2e	%.2f	%s	%s\n' % \
					(this_sub2.ID, 	this_sub2.m/h0, this_r/h0, 'MW', print_run)
					file_sub_txt.write(file_sub_line)
			
			print 'LG2 %d) has %d subs over %d particles\n' % (subrun_i, count_sub2, np_sub_min)
		
			if this_lg2.nsub > n_sub_min:
				moment_inertia(coords, masses)

		#		print 'id: %d, np= %d, m= %e, d = %.2f' % (sub_i, this_sub2.npart, this_sub2.m, this_sub2.distance(this_lg2.x))

	print ' ----------------------------------- '


	



