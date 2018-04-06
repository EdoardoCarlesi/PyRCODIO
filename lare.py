#!/usr/bin/python

import numpy as np
import os

from libio.read_ascii import * 
from libcosmo.halo import Halo
from libcosmo.find_halos import *
from libcosmo.utils import *
#from libcosmo.std_lg_plot import *
from libio.lare import *
from config import *

home_dir = '/z/carlesi/CLUES/'
#ginn_dir = '/ginnunagap/'
ginn_dir = 'ginnungagap/ginnungagap/'
lare_dir = 'GridProperties/'

#snapshot = 'snapshot_054.0000.z0.000.AHF_halos'
#snapshot = 'snapshot_054.AHF_halos'
snapshot = 'merged_054.AHF_halos'
#snapshot = 'snapshot_054.z0.000.AHF_halos'

#home_dir = '/home/eduardo/'
hubble = 0.671		
h0 = hubble 
facMpc = 1000.

generate_lare='false'
#generate_lare='true'
do_plots='false'
#do_plots='true'

#env_type="std"
env_type="512box100"
#env_type="HESTIA"

resolution='512'
#resolution='1024'

if do_plots == "true":
	outp_dir = 'output_analysis/'
else:
	outp_dir = 'output_txt/'

settings = Settings(home_dir, outp_dir, env_type, resolution, snapshot)

# Factor rescaling the particle distribution in plots
rescale = 10

#find_virgo='false'
find_virgo='true'

# Max/min distance from Virgo in Mpc/h
d_virgo_max = 20000. / h0
d_virgo_min = 7000. / h0

# This one specifies the number and the types of LSS modes to be analyzed
lss_init = 1
lss_end = 14

# How many sub-runs per LS realisation
run_init = 0
run_end = 1

# If running on all seeds
ice_init = 8
ice_end = 12
gin_init = 8
gin_end = 20

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

lg_model = LocalGroupModel(radius, iso_radius, r_max, r_min, m_max, m_min, ratio_max, vrad_max)

# Subhalo identification criterion
fac_r = 1.5
part_min = 30

# Define more units and stuff
extra_r = 500.0	# in kpc/h units - this is the size of the additional shell to be placed around the 

lg_dummy = LocalGroup('00')

file_lg_header = lg_dummy.header()

if find_virgo == "true":
	file_lg_header.rstrip()	
	file_lg_header += 'R_Virgo(11)[Mpc] M_Virgo(12) M_Virgo(<5Mpc)(13)'

file_lg_header += '\n'
file_lare_header = '# SimuCode(1) R_z0(2) X_z0(3) Y_z0(4) Z_z0(5) [Mpc/h] R_ini(6) X_ini(7) Y_ini(8) Z_ini(9) [Mpc/h]\n'
	
#return_ascii_files(env_type, home_dir, data_dir, resolution, outp_dir)

# LG Candidates properties .txt file
file_lg_txt = open(settings.file_lg_name, 'wb')
file_lgall_txt = open(settings.file_lgall_name, 'wb')
file_codes_txt = open(settings.file_codes_name, 'wb')

# Write headers
file_lgall_txt.write(file_lg_header)
file_lg_txt.write(file_lg_header)

if generate_lare == "true":
	# Save data on the lagrangian region here
	file_lare_txt = open(settings.file_lare_name, 'wb')
	file_lare_txt.write(file_lare_header)

base_lss = []

# This is generating all the possible folders and subfolders with the main seed and small scale seed
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

for irun in range(lss_init, lss_end):
	base_run = base_lss[irun]

	# LaRe finder - script
	lare_find_sh = home_dir + lare_dir + '/scripts/lare_find.sh'

	# LaRe generator - ginnungagap script
	lare_gen_sh = home_dir + ginn_dir + 'zoom/gen_mask.sh'
	lare_dir_sh = home_dir + ginn_dir + 'zoom/'

	# File names common to all configurations
	file_ic = 'ic_zoom_512_100.000_' + base_run

	# when several LG-like pairs are found, get the first pair (0) second pair (2) etc.
	ind_lg = 0

	for run_i in range(run_init, run_end):
		run_num = '%02d' % run_i
                print_run = base_run + '_' + run_num
		
		settings.init_files(base_run, run_num)		

		print 'Reading in AHF file: ', settings.file_ahf_in

		if os.path.exists(settings.file_ahf_in):
			print ''
			print 'Reading in AHF file: ', settings.file_ahf_in
		
			ahf_all = read_ahf(settings.file_ahf_in)

			if find_virgo == "true":
				(x0, m0, mtotVirgo) = locate_virgo(ahf_all)

			these_lg = find_lg(ahf_all, lg_model)
			n_lgs = int(len(these_lg) / 2)

			print 'Found a total of %d LG pairs for the %s run.' % (n_lgs, base_run)

			if n_lgs > 0:
				# Set to a very high value the variable "goodness of LG". The smallest the value, the closer to the real LG.
				# This is needed in case there are more than one LG candidates in the simulation.
				# Resetting all the indexes
				rating = 1000
				file_lg_line = ''
				best_lg_line = 'None.'
				save_lg = "false"

				for ilg in range(0, n_lgs):
					mw = these_lg[2 * ilg]			
					m31 = these_lg[2 * ilg + 1]		

					lg = LocalGroup(print_run)
					lg.init_halos(m31, mw)				
					file_lg_line = lg.info()


					if find_virgo == "true":
						file_lg_line.rstrip()
						dVirgo = vec_distance(lg.com, x0) / h0
						file_lg_line += "  %5.2f  %5.2e  %5.2e" % (dVirgo/facMpc, m0/h0, mtotVirgo/h0)
						
					if dVirgo > d_virgo_min and dVirgo < d_virgo_max:
						virgo_condition = True
					else:
							print 'Distance from Virgo is: %f, too large or too small' % dVirgo
							virgo_condition = False
		
					print file_lg_line
					file_lg_line += "\n"
		
					# Do an additional check on the LG candidates before appending it to the global list of LG candidates
					#if vel_r < vrad_max and m12 < m_max and m_ratio < ratio_max and virgo_condition:
					file_lgall_txt.write(file_lg_line)
					lg.rating()

					# Now check for the likeliest LG-like candidate
					if lg.rating < rating and (lg.LG1.npart > part_min) and (lg.LG2.npart > part_min):
						print 'Old rating: %f new rating %f this index %d' % (rating, lg.rating, ilg)
						rating = lg.rating
						best_lg_line = file_lg_line
						best_lg = lg
						save_lg = True

				# The loop on the local groups has finished, we kept only the best one
				if save_lg == True:
					print 'Selected LG candidate: \n'
					print best_lg_line				
					# Now plot / analyse and find LaRe of the best candidate ONLY
					# First write the relevant quantities
					file_lg_txt.write(best_lg_line)
					file_codes_txt.write(settings.base_run+"   ")
					#print best_lg1.npart, best_lg2.npart, best_lg2.m, best_lg1.m
		
					# Plot the particle distributions in 3D slices around the LG and in the LV
					if do_plots == "true" and lg.LG1.npart > npart_min and lg.LG2.npart > npart_min:
						print 'Plotting particle distributions to file ', settings.get_plot_out
						plot_lglv(settings.file_z0_in, ahf_all, settings.plot_out, best_lg1, best_lg2, x0, rescale)
			
					# Only add the LaRe if for the likeliest pair 
					if generate_lare == "true":
						lare_r = extra_r + best_lg1.distance(best_lg.LG2.x) + best_lg.LG1.r + best_lg.LG2.r

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

				else:
					print 'No valid LG candidates in this realisation.\n'


