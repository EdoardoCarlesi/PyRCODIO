#!/usr/bin/python

import numpy as np
import os

from config import *
from libio.read_ascii import *
from libcosmo.utils import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.lg_plot import *

resolution='2048'
#resolution='1024'
run_init = 3
run_end = 4

subrun_init = 0
subrun_end = 10

base_path = '/home/eduardo/CLUES/'
outp_path = 'output/'
env_type = 'zoom'

ahf_snap = 'snapshot_054.0000.z0.000.AHF_halos'; ahf_snap_alt = 'snapshot_054.0000.z0.001.AHF_halos'
snapshot = 'snapshot_054'

#do_plots = "true"
do_plots = "false"

reduce_fac = 8
plot_pos = "false"

hubble = 0.67		
part_min = 25000

# Local group selection parameters
center = [50000., 50000., 50000.]

# Allocate LG Models
all_runs = simu_runs()
(all_lg_models, hash_run) = lg_models()

lg_dummy = LocalGroup(all_runs[0])
sub_dummy = SubHalos(0, '00', '00', [0, 0])

file_sub_header = sub_dummy.header()
file_lg_header = lg_dummy.header()

settings = Settings(base_path, outp_path, env_type, resolution, snapshot)
settings.box_center = center

# Subhalo identification criterion
r_sub_min = 10.0   # Minimum distance from the main halo - if it's too close it might be double counting the host
m_sub_min = 7.e+9 # This mass is in REAL Msun units, without the /h
n_sub_min = 4 # Minimum number of subhalos required to compute the moment of inertia
fac_r = 1.0
np_sub_min = 50


# when several LG-like pairs are found, get the first pair (0) second pair (2) etc.
ind_lg = 0

for run_j in range(run_init, run_end):
	base_run = all_runs[run_j]
	this_run = hash_run[base_run]
	lg_model = all_lg_models[this_run]

	settings.base_run = base_run
	settings.re_init()

	#print base_run, settings.base_run
	(file_lg_name, file_sub_name) = settings.get_zoom_output()

	file_lg_txt = open(file_lg_name, 'wb')
	file_lg_txt.write(file_lg_header)
	file_sub_txt = open(file_sub_name, 'wb')
	file_sub_txt.write(file_sub_header)

	# Subhalo mass function variables
	file_png_mfs1 = outp_path + 'mf_01_' + resolution+'_' + base_run + '.png'
	file_png_mfs2 = outp_path + 'mf_02_' + resolution+'_' + base_run + '.png'
	x_mf1 = []; 	y_mf1 = []
	x_mf2 = []; 	y_mf2 = []
	n_mf = 0

	for subrun_i in range(subrun_init, subrun_end):
		run_num = '%02d' % subrun_i
		settings.init_files(base_run, run_num)

		file_png_name = outp_path + 'lg_' + resolution+'_' + base_run + '_' + run_num + '.png'
		this_file_ahf = settings.ahf_path + ahf_snap; 		this_file_ahf_alt = settings.ahf_path + ahf_snap_alt
		this_file_gad = settings.file_z0
		print_run = base_run + '_' + run_num
		good_lgs = 0
		
		if os.path.exists(this_file_ahf):
				print 'Reading in AHF file: ', this_file_ahf
				ahf_all = read_ahf(this_file_ahf)
				these_lg = find_lg(ahf_all, lg_model)
				n_lgs = int(len(these_lg))

		elif os.path.exists(this_file_ahf_alt):
				print 'Reading in alternative AHF file: ', this_file_ahf_alt
				ahf_all = read_ahf(this_file_ahf_alt)
				these_lg = find_lg(ahf_all, lg_model)
				n_lgs = int(len(these_lg))
		else:
				print 'No AHF file: ', this_file_ahf
				n_lgs = 0;
		
		print 'Found a total of %d LG pairs.' % (n_lgs)

		if n_lgs > 0:
			# If there are more candidates we need to find the right one
			rating = 1000
			for ind in range(0, n_lgs):
				lg = these_lg[ind]
				lg.c_box = settings.box_center

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

				file_sub_line  = subs1.all_info("part", 2 * np_sub_min)
				file_sub_line += subs2.all_info("part", 2 * np_sub_min)
				file_sub_line += '\n ------------------------ \n'
				file_sub_txt.write(file_sub_line)

				(x_m, y_n) = subs1.mass_function()
				x_mf1.append(x_m); 	y_mf1.append(y_n)

				(x_m, y_n) = subs2.mass_function()
				x_mf2.append(x_m); 	y_mf2.append(y_n)
				n_mf += 1
			
			# Do some plots at the end
			if do_plots == "true":
				plot_lg(this_file_gad, file_png_name, best_lg.LG1, best_lg.LG2, reduce_fac, 1, plot_pos)
	
				if (subrun_i == subrun_end-1):
					plot_massfunctions(x_mf1, y_mf1, n_mf, file_png_mfs1)
					plot_massfunctions(x_mf2, y_mf2, n_mf, file_png_mfs2)
					del x_m	;	del y_n
					del x_mf1; 	del x_mf2
					del y_mf1; 	del y_mf2
					n_mf = 0

			# Compute halo evolution

			

