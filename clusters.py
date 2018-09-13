#!/usr/bin/python

from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.particles import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *
from libcosmo.grid import *
from config import *
import time
import pickle
import os.path

base_path='/home/eduardo/CLUES/DATA/'
#sub_dir='TEST'
sub_dir='SIMU_CF3'
extra_tag='CF3k_YH_v1'

plot_slice = True
#plot_slice = False
find_clusters = False
#find_clusters = True

ini_num=70
end_num=75

nbins = 450
f_rescale = 1


box = '400'; num = '256'; box_size=400.0e+3
#box = '500'; num = '256'; box_size=500.0e+3
#box = '600'; num = '256'; box_size=600.0e+3

plot_side = 100.0e+3; thickn = 5000.0; units = 'kpc'
#box_size = 500.0; plot_side = 100.0; thickn = 9.5; units = 'Mpc'

ahf_name='snapshot_054.AHF_halos'
snap_name='snapshot_054'

for run_num in range(ini_num, end_num):

	f_out='plot_'+snap_name+'_'+extra_tag+'_'+str(run_num)+'_bins'+str(nbins)+'.png'

	#print 'N = ', run_num, ' - ', f_out
	print ' ---- ', run_num, ' ----- '

	run_num=str(run_num)
	base_dir=base_path+sub_dir+'/'+box+'/'+num+'/'+extra_tag+'/'+run_num+'/'

	box_center = [0.5 * box_size, 0.5 * box_size, 0.5 * box_size]

	if os.path.isfile(base_dir + ahf_name) and find_clusters == True:
		#print base_dir + ahf_name
		all_halos = read_ahf(base_dir + ahf_name)
		locate_clusters(all_halos, box_center)

	if os.path.isfile(base_dir + snap_name + '.0') and plot_slice == True:
		print base_dir + snap_name
		plot_rho(base_dir + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units)



