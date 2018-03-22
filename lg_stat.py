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

i_init = 0
i_end = 100
g_init = 0
g_end = 40
s_init = 0
s_end = 10

base_path = '/home/eduardo/CLUES/DATA/'
outp_path = 'output/'

ahf_snap = 'snapshot_054.0000.z0.000.AHF_halos'
#ahf_snap = 'snapshot_t1_054.0000.z0.000.AHF_halos'
snapshot = 'snapshot_054'

do_plots = "true"
#do_plots = "false"

reduce_fac = 8
plot_pos = "false"

# this is going to multiply kpc/h distances - thus it is rescaled by a factor of 1000.
hubble = (67.1 / 1000.)		
base_path = '/home/eduardo/CLUES/DATA/' 

# Subhalo identification criterion
fac_r = 1.2
np_sub_min = 30

# LGs will be appendend to this as they are found
m_bins = []
n_bins = []

# when several LG-like pairs are found, get the first pair (0) second pair (2) etc.
ind_bins = 0

this_file_png = base_path + outp_path + 'lg_stat_' + resolution + '.png'

for i_i in range(i_init, i_end):
	i_num = '%02d' % i_i

	for i_g in range(g_init, g_end):
		g_num = '%02d' % i_g
		base_num = i_num + '_' + g_num
		this_path = base_path + resolution + '/' + base_num
	
		if os.path.exists(this_path):

			this_file_txt = base_path + outp_path + 'lg_candidates_' + resolution + '_' + base_num + '.txt' 
 
			if os.path.exists(this_file_txt):
				print ind_bins, ') Reading in TXT LG file: ', this_file_txt
				lgs = read_lgs(this_file_txt)
				(this_m_bin, this_n_bin) = bin_lg_sub(lgs)
				
				print this_m_bin
	
				m_bins.append(this_m_bin)
				n_bins.append(this_n_bin)

				ind_bins += 1

			else:
				print 'TXT LG file not found: ', this_file_txt
			
plot_lg_bins(m_bins, n_bins, this_file_png)


#print m_bins
#print n_bins


'''
			for i_s in range(s_init, s_end):
				s_num = '%02d' % i_s
				print_num = base_num + '_' + s_num
				#file_png_name = outp_path + 'stat_' + resolution+'_' + print_num + '.png'
				this_file_ahf =  this_path + '/' + s_num + '/' + ahf_snap
			
				if os.path.exists(this_file_ahf):
					print 'Reading in AHF file: ', this_file_ahf
					#ahf_all = read_ahf(this_file_ahf)
	
'''

