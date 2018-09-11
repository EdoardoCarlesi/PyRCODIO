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

ini_num=2
end_num=10

nbins = 2500
f_rescale = 1

#lg_00_06_00.pkl

#lg_base='00_06'
#lg_base='17_10'
#lg_base='01_12'
#lg_base='09_18'
#lg_base='34_13'
#lg_base='45_17'
lg_base='55_02'
#lg_base='37_11'

box = '100'; num = '2048'
box_size = 100.0; plot_side = 1.0; thickn = 2.0; units = 'Mpc'
#box_size = 100.0e+3; plot_side = 10.0e+3; thickn = 5000.0; units = 'kpc'
snap_name='snapshot_054'

for run_num in range(ini_num, end_num):

	subrun = '%02d' % run_num
	lg_tag = lg_base + '_' + subrun
	pkl_name = 'saved/lg_' + lg_tag + '.pkl'	
	f_out='plot_'+snap_name+'_'+lg_tag+'_LG.png'

	if os.path.isfile(pkl_name):
		pkl_file = open(pkl_name, 'r')

#	print 'N = ', lg_tag, ' - ', f_out

		this_lg = pickle.load(pkl_file)
		#print	this_lg.info()

		base_dir=base_path+'/'+num+'/'+lg_base+'/'+subrun+'/'

		box_center = this_lg.get_com()
		
		for i in range(0, 3):
			box_center[i] = box_center[i] * 1.e-3
		#box_center = [0.5 * box_size, 0.5 * box_size, 0.5 * box_size]
		
		print box_center
#		print base_dir + snap_name
		if os.path.isfile(base_dir + snap_name):
			print base_dir + snap_name
			plot_rho(base_dir + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units)



