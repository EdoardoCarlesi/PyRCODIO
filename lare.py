#!/usr/bin/python

import numpy as np
import os

from libcosmo.halo import Halo
from libcosmo.read_ahf import * 
from libcosmo.find_halos import *
from libcosmo.utils import *
from libcosmo.std_lg_plot import *
from libics.lare_tools import *

#home_dir = '/z/carlesi/'
home_dir = '/home/eduardo/'
hubble = 0.671		
h0 = hubble 

generate_lare='false'
#generate_lare='true'
#generate_base='false'
generate_base='true'

#do_plots='false'
do_plots='true'
# Factor rescaling the particle distribution in plots
rescale = 8	

#find_virgo='false'
find_virgo='true'

env_type="std"
#env_type="HESTIA"
resolution='512'
#resolution='1024'

# This one specifies the number and the types of LSS modes to be analyzed
lss_init = 1
lss_end = 14

# How many sub-runs per LS realisation
run_init = 0
run_end = 1

# If running on all seeds
ice_init = 0
ice_end = 80
gin_init = 0
gin_end = 30

# Local group selection parameters
center = [50000., 50000., 50000.]
radius = 8000. 
iso_radius = 2300. * h0
r_max = 2000. * h0
r_min = 400. * h0
m_min = 1.e+12 * h0 
m_max = 8.e+12 * h0

ratio_max = 3.
vrad_max = -1.0

# Subhalo identification criterion
fac_r = 1.5
npart_min = 150

# Define more units and stuff
extra_r = 1000.0	# in kpc/h units - this is the size of the additional shell to be placed around the 

file_lg_header = '#SimuCode(0) ID_M31(1) ID_MW(2) M_M31(3)[Msun] M_MW(4)[Msun] R_MwM31(5)[kpc] Vrad_M31(6)[phys, km/s] PosCoM(8-9-10)[Mpc/h]'
	
if find_virgo == "true":
	file_lg_header += 'R_Virgo(11)[Mpc] M_Virgo(12) M_Virgo(<5Mpc)(13)'

file_lg_header += '\n'
file_lare_header = '# SimuCode(1) R_z0(2) X_z0(3) Y_z0(4) Z_z0(5) [Mpc/h] R_ini(6) X_ini(7) Y_ini(8) Z_ini(9) [Mpc/h]\n'
	
if env_type == "HESTIA" :
	base_path = '/z/carlesi/CLUES/'
	sub_path = 'DM_ONLY/'
	sub_path_ahf = 'AHF_output/'
	sub_path_ics = 'ICs/'
	sub_path_snap = 'output/snapdir_127/'
	file_lg_name = 'lg_candidates_HESTIA_'+resolution+'.txt'
	file_lare_name = 'lg_candidates_LaReHESTIA_'+resolution+'.txt'
	file_lgall_name = 'lg_candidates_HESTIAall_'+resolution+'.txt'
	file_codes_name = 'lg_codes_HESTIA_'+resolution+'.txt'
else:
	base_path = '/home/eduardo/CLUES/DATA/' 
	file_lg_name = 'lg_candidates_LGF_'+resolution+'.txt'
	file_lgall_name = 'lg_candidates_LGFall_'+resolution+'.txt'
	file_lare_name = 'lg_candidates_LaReLGF_'+resolution+'.txt'
	file_codes_name = 'lg_codes_'+resolution+'.txt'


# LG Candidates properties .txt file
file_lg_txt = open(base_path + file_lg_name, 'wb')
file_lgall_txt = open(base_path + file_lgall_name, 'wb')
file_codes_txt = open(base_path + file_codes_name, 'wb')

file_lgall_txt.write(file_lg_header)
file_lg_txt.write(file_lg_header)

if generate_lare == "true":
	# Save data on the lagrangian region here
	file_lare_txt = open(base_path + file_lare_name, 'wb')
	file_lare_txt.write(file_lare_header)

if generate_base == 'true':
	base_lss = []
	for iice in range(ice_init, ice_end): 
		for igin in range(gin_init, gin_end):

			if iice < 10:
				nice='0'+str(iice)
			else:
				nice=str(iice)

			if igin < 10:
				ngin='0'+str(igin)
			else:
				ngin=str(igin)

			lss_init = 0
			lss_end = (ice_end - ice_init) * (gin_end - gin_init)
			this_lss = nice + '_' + ngin
			#print this_lss
			base_lss.append(this_lss)
else:
	base_lss = ['00_06', '01_11', '01_13', '17_00', '17_01', '17_02', '17_03', '17_04', '24_06', '25_02', '34_04', '55_02', '68_08']

for irun in range(lss_init, lss_end):
	base_run = base_lss[irun]

	if env_type == "HESTIA" :
		snapshot = 'HESTIA_100Mpc_512_'+base_run+'.127.z0.000.AHF_halos'
		file_z0 = 'snapshot_127'
		base_out = 'lare_z0_' + base_run + '.dat'
		n_ic = 2
		n_z0 = 8

		# LaRe finder - script
		lare_find_sh = home_dir + 'CLUES/GridProperties/scripts/lare_find.sh'

		# LaRe generator - ginnungagap script
		lare_gen_sh = home_dir + 'CLUES/ginnungagap/ginnungagap/zoom/gen_mask.sh'
		lare_dir_sh = home_dir + 'CLUES/ginnungagap/ginnungagap/zoom/'
	else:
		#snapshot = 'snapshot_054.0000.z0.000.AHF_halos'
		#snapshot = 'snapshot_054.AHF_halos'
		snapshot = 'snapshot_054.z0.000.AHF_halos'
		file_z0 = 'snapshot_054'
		base_out = 'lare_z0_' + base_run + '.dat'
		n_ic = 2
		n_z0 = 1

		# LaRe finder - script
		lare_find_sh = home_dir + 'CLUES/GridProperties/scripts/lare_find.sh'

		# LaRe generator - ginnungagap script
		lare_gen_sh = home_dir + 'CLUES/ginnungagap/zoom/gen_mask.sh'
		lare_dir_sh = home_dir + 'CLUES/ginnungagap/zoom/'

	# File names common to all configurations
	file_ic = 'ic_zoom_512_100.000_' + base_run
	plot_out = base_run + '_particles_LG_LV.png'

	# LGs will be appendend to these as they are found
	all_lgs = []
	best_lgs = []

	# when several LG-like pairs are found, get the first pair (0) second pair (2) etc.
	ind_lg = 0

	for run_i in range(run_init, run_end):
		run_num = '%02d' % run_i
		
		if resolution == "1024" or resolution == "2048":
			base_run = base_run + '/' + run_num

		if env_type == "HESTIA" :
			file_ahf_in = base_path + sub_path + '/' + base_run + '/' + sub_path_ahf + snapshot
			file_z0_in = base_path + sub_path + '/' + base_run + '/' + sub_path_snap + file_z0
			file_ic_in = base_path + sub_path + '/' + base_run + '/' + sub_path_ics + file_ic
			#file_ic_in = base_path + '/HESTIA/ICs/' + file_ic
		else:
			file_ahf_in = base_path + resolution + '/' + base_run + '/' + snapshot
			file_z0_in = base_path + resolution + '/' + base_run + '/' + file_z0
			file_ic_in = base_path + resolution + '/ICs/' + file_ic

		if os.path.exists(file_ahf_in):
			
			print ''
			print 'Reading in AHF file: ', file_ahf_in
		
			ahf_all = read_ahf(file_ahf_in)

			if find_virgo == "true":
				(x0, m0, mtotVirgo) = locate_virgo(ahf_all)

			these_lg = find_lg(ahf_all, center, radius, iso_radius, m_min, r_min, r_max)
			n_lgs = int(len(these_lg) / 2)

			print 'Found a total of %d LG pairs for the %s run.' % (n_lgs, base_run)

			if n_lgs > 0:
				# Set to a very high value the variable "goodness of LG". The smallest the value, the closer to the real LG.
				# This is needed in case there are more than one LG candidates in the simulation.
				fac_lg0 = 1000.

				for ind_lg in range(0, n_lgs):
					#file_out = base_path + resolution + '/' + base_run + '/' + base_out + '.' + str(ind_lg)
					# This will be overwritten for the best LG candidate only
					if env_type == "HESTIA" :
						file_out = base_path + '/HESTIA/LaRe/' + base_out
					else:
						file_out = base_path + resolution + '/' + base_run + '/' + base_out
				
					ind_lg = 2 * ind_lg
					lg1 = these_lg[ind_lg]
					lg2 = these_lg[ind_lg+1]
					m12 = lg1.m + lg2.m

					fac_lg12 = rate_lg_pair(lg1, lg2, center)	

					# Pair properties
					lg_com = center_of_mass([lg1.m, lg2.m], [lg1.x, lg2.x])
					d12 = lg1.distance(lg2.x) 
					vel_r = vel_radial(lg1.x, lg2.x, lg1.v, lg2.v)
					#print vel_r, d12, hubble
					vel_r += hubble * d12 / 10.0
					#print vel_r
	
					if lg_com[0] > 500.:
						facMpc = 1000.
					else:
						facMpc = 1.

					file_lg_line = '%s  %ld   %ld   %.2e   %.2e   %7.2f   %7.2f  %5.2f  %5.2f  %5.2f' % \
			(base_run, lg1.ID, lg2.ID, lg1.m/h0, lg2.m/h0, d12/h0, vel_r, lg_com[0]/facMpc, lg_com[1]/facMpc, lg_com[2]/facMpc)

					if find_virgo == "true":
						dVirgo = vec_distance(lg_com, x0) / h0
						file_lg_line += "  %5.2f  %5.2e  %5.2e" % (dVirgo/facMpc, m0/h0, mtotVirgo/h0)
				
					print file_lg_line
					file_lg_line += "\n"

					if lg1.m > lg2.m:
						m_ratio = lg1.m / lg2.m
					else:
						m_ratio = lg2.m / lg1.m

					# Do an additional check on the LG candidates before appending it to the global list of LG candidates
					if vel_r < vrad_max and m12 < m_max and m_ratio < ratio_max:
						file_lgall_txt.write(file_lg_line)
						all_lgs.append(lg1)
						all_lgs.append(lg2)
						print file_lg_line			
		
						# Now check for the likeliest LG-like candidate
						if fac_lg12 < fac_lg0:
							best_lg1 = lg1
							best_lg2 = lg2
							best_lg_line = file_lg_line
							best_lg_com = lg_com
							fac_lg0 = fac_lg12 		

				# The loop on the local groups has finished, we kept only the best one
				# if n_lg > 0 end
				print 'Selected LG candidate: \n'
				print best_lg_line				

				# Now plot / analyse and find LaRe of the best candidate ONLY
				# First write the relevant quantities
				file_lg_txt.write(best_lg_line)
				file_codes_txt.write(base_run+"   ")
			
				# Reset the indicator
		
				# Plot the particle distributions in 3D slices around the LG and in the LV
				if do_plots == "true" and lg1.npart > npart_min and lg2.npart > npart_min:
					print 'Plotting particle distributions to file ', plot_out
					plot_lglv(file_z0_in, ahf_all, plot_out, best_lg1, best_lg2, x0, rescale)
			
					# Only add the LaRe if for the likeliest pair 
					if generate_lare == "true":
						lare_r = extra_r + best_lg1.distance(best_lg2.x) + best_lg1.r + best_lg2.r

						if lare_r > 100.:
							facMpc = 1000.
						else:
							facMpc = 1.

						str_lare_z0 = '%s   %.3f   %.3f   %.3f   %.3f   ' % \
						(base_run, lare_r/facMpc, best_lg_com[0]/facMpc, best_lg_com[1]/facMpc, best_lg_com[2]/facMpc)
	
						if os.path.exists(file_ic_in) or os.path.exists(file_ic_in + '.0'):
							print 'Found IC file: %s to determine LaRe\n' % file_ic_in

								# This one finds the lagrangian region 
							all_lare_x = find_lare(best_lg_com, lare_r, file_ic_in, file_z0_in, \
										n_ic, n_z0, file_out, lare_find_sh)

							#print 'LARE: %s, %s' % (all_lare_x, file_ic_in)
							n_lare_x = len(all_lare_x)

							lare_x = str(all_lare_x[n_lare_x-4]) + ' ' + str(all_lare_x[n_lare_x-3]) + ' ' \
							 + str(all_lare_x[n_lare_x-2]) + ' ' + str(all_lare_x[n_lare_x-1]) + ' '

							if all_lare_x[n_lare_x-4] < 10.0:
								print ' *************** WARNING ************* '
								print 'Problem with coordinate: %f ' % all_lare_x[n_lare_x-4]
								print ' ************************************* '

							str_lare_init = '%.3f   %.3f   %.3f   %.3f\n' % \
					(all_lare_x[n_lare_x-1], all_lare_x[n_lare_x-4], all_lare_x[n_lare_x-3], all_lare_x[n_lare_x-2])

							file_lare_txt.write(str_lare_z0 + str_lare_init)
							print 'Lare X: %s' % lare_x

							# Now generate a mask
							gen_mask(base_run, lare_x, lare_gen_sh, lare_dir_sh)
						else:
							print 'Could not find IC file: %s to generate LaRe\n' % file_ic_in



