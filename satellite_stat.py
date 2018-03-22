#!/usr/bin/python

import numpy as np
import os

from libcosmo.utils import *
from libcosmo.halo import *
from libcosmo.read_ascii import *
from libcosmo.find_halos import *
from libcosmo.std_lg_plot import *

resolution='2048'
#resolution='2048b'
#resolution='1024'
run_init = 4
run_end = 5

subrun_init = 0
subrun_end = 10

base_path = '/home/eduardo/CLUES/DATA/'
outp_path = 'output/'

all_runs = ['00_06', '55_02', '34_13', '45_17', '17_10', '17_01', '64_14', '01_12']

ahf_snap = 'snapshot_054.0000.z0.000.AHF_halos'
#ahf_snap = 'snapshot_t1_054.0000.z0.000.AHF_halos'
snapshot = 'snapshot_054'

#do_plots = "true"
do_plots = "false"

reduce_fac = 8
plot_pos = "false"

hubble = 0.67		
base_path = '/home/eduardo/CLUES/DATA/' 

# Local group selection parameters
center = [50000., 50000., 50000.]

# 17_01
#d_max = 8000.; r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 7.0e+11; m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0

# 00_06
#d_max = 5000.; r_iso = 1900.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -1.0

# 55_02, 34_13
d_max = 7000.; r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 9.0e+11; m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0


lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
lg_dummy = LocalGroup(all_runs[0], hubble)

file_lg_header = lg_dummy.header()
file_sub_header = '# ID_sub(1) 	M_sub(2)[Msun] 	R_sub(3)[Mpc]	Host(4) 	SimuCode(5)\n' 
part_min = 25000

# Subhalo identification criterion
r_sub_min = 10.0   # Minimum distance from the main halo - if it's too close it might be double counting the host
m_sub_min = 7.e+9 # This mass is in REAL Msun units, without the /h
n_sub_min = 4 # Minimum number of subhalos required to compute the moment of inertia
fac_r = 1.0
np_sub_min = 50

# All LG structures will be appended here
all_lgs = []

# when several LG-like pairs are found, get the first pair (0) second pair (2) etc.
ind_lg = 0

for run_j in range(run_init, run_end):
	base_run = all_runs[run_j]

	file_lg_name = outp_path + 'lg_candidates_'+resolution+'_'+base_run+'.txt'
	file_sub_name = outp_path + 'lg_subhalos_'+resolution+'_'+base_run+'.txt'
	file_lg_txt = open(base_path + file_lg_name, 'wb')
	file_lg_txt.write(file_lg_header)
	file_sub_txt = open(base_path + file_sub_name, 'wb')
	file_sub_txt.write(file_sub_header)

	for subrun_i in range(subrun_init, subrun_end):
		run_num = '%02d' % subrun_i
		file_png_name = outp_path + 'lg_' + resolution+'_' + base_run + '_' + run_num + '.png'
		this_file_ahf = base_path + resolution + '/' + base_run + '/' + run_num + '/' + ahf_snap
		this_file_gad = base_path + resolution + '/' + base_run + '/' + run_num + '/' + snapshot
		print_run = base_run + '_' + run_num
		good_lgs = 0
		
		if os.path.exists(this_file_ahf):
				print 'Reading in AHF file: ', this_file_ahf
				ahf_all = read_ahf(this_file_ahf)
				these_lg = find_lg(ahf_all, lg_model)
				n_lgs = int(len(these_lg) / 2)
		else:
				print 'No AHF file: ', this_file_ahf
				n_lgs = 0;
		
		print 'Found a total of %d LG pairs.' % (n_lgs)

		if n_lgs > 0:
			# If there are more candidates we need to find the right one
			rating = 1000
			for ind in range(0, n_lgs):
				mw = these_lg[2 * ind]			
				m31 = these_lg[2 * ind + 1]		

				lg = LocalGroup(print_run, hubble)
				lg.init_halos(m31, mw)				
				 
				if lg.rating() < rating and (lg.LG1.npart > part_min) and (lg.LG2.npart > part_min):
					print 'Old rating: %f new rating %f this index %d' % (rating, lg.rating, ind)
					good_lgs += 1
					rating = lg.rating
					best_lg = lg

			if good_lgs > 0:
				print 'Best LG: ', best_lg.info()

				# Pair properties
				file_lg_line = best_lg.info()

				#if vel_r < vrad_max and m12 < m_max:
				file_lg_txt.write(file_lg_line)

			# Now take care of the substructure
			if good_lgs > 0:
				these_sub1 = find_halos(best_lg.LG1, ahf_all, fac_r * best_lg.LG1.r)
				subs1 = SubHalos(best_lg.LG1, these_sub1, print_run, 'M31')
				subs1.anisotropy("part", np_sub_min)
				subs1.basis_eigenvectors("inertia")
				nc = subs1.basis_eigenvectors("inertia")

				these_sub2 = find_halos(best_lg.LG2, ahf_all, fac_r * best_lg.LG2.r)
				subs2 = SubHalos(best_lg.LG2, these_sub2, print_run, 'MW')
				subs2.anisotropy("part", np_sub_min)
				subs2.basis_eigenvectors("inertia")
				nc = subs2.basis_eigenvectors("inertia")

			# Do some plots at the end
			if do_plots == "true":
				plot_lg(this_file_gad, file_png_name, best_lg.LG1, best_lg.LG2, reduce_fac, 1, plot_pos)



