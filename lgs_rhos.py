#!/usr/bin/python

from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halos import *
from libcosmo.find_halo import *
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

ini_num=0
end_num=1

nbins = 10000
#nbins = 512
f_rescale = 1.0 / 128.0

#all_lg_base = simu_runs()

#all_lg_base=['00_06']
all_lg_base=['37_11']
#lg_base='17_10'
#lg_base='17_13'
#lg_base='01_12'
#lg_base='09_18'
#lg_base='34_13'
#lg_base='45_17'
#lg_base='55_02'
#lg_base='64_14'

#box = '100'; resolution = '4096'; n_files = 1
#box = '100'; resolution = '2048'; n_files = 1
box = '100'; resolution = '1024'; n_files = 4
#box_size = 100.0; plot_side = 1.75; thickn = 3.0; units = 'Mpc'
#box_size = 100.0; plot_side = 1.5; thickn = 3.0; units = 'Mpc'
#box_size = 100.0e+3; plot_side = 10.0e+3; thickn = 1.5e+3; units = 'kpc'
#box_size = 100.0e+3; plot_side = 10.0e+3; thickn = 5000.0; units = 'kpc'
box_size = 100.0; plot_side = 1.25; thickn = 2.5; units = 'Mpc'

#snap_name='snapshot_054'
snap_name='snapshot_054'
#snap_name='ic_arepo_2048_100.000_37_11_00.0'
base_path = '/home/eduardo/CLUES/DATA/1024/00/'

f_out = 'test_lg_gas'
box_center = [46.6, 50.7, 47.8]
#box_center = [46600, 50750, 47800]
#plot_rho(base_path + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units, n_files)
#plot_rho(base_path + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units, n_files)
#plot_rho_slice(base_path + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units, n_files)

'''
for lg_base in all_lg_base:

        # Just choose one LG center to plot all the other sims - 04 is a good average choice
	pkl_name = 'saved/lgs_' + resolution + '_' + lg_base + '_04.pkl'	
	pkl_file = open(pkl_name, 'rb')
	this_lg = pickle.load(pkl_file)
	box_center = this_lg.geo_com()

	for i in range(0, 3):
		box_center[i] = box_center[i] * 1.e-3

	for run_num in range(ini_num, end_num):
	
		subrun = '%02d' % run_num
		lg_tag = lg_base + '_' + subrun
		f_out = 'lg_plot_' + resolution + '_' + lg_tag +'_LG.png'
		base_dir = base_path + '/' + resolution + '/' + lg_base + '/' + subrun + '/'

		print('N = ', lg_tag, ' - ', f_out)

		if os.path.isfile(base_dir + snap_name) or os.path.isfile(base_dir + snap_name + '.0'):
			plot_rho(base_dir + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units, n_files)

		else:
			print('File not found: ', base_dir + snap_name)
'''
	
'''
	This (older) version loops over all the LGs and finds each time all the LGs
	
	subrun = '%02d' % run_num
	lg_tag = lg_base + '_' + subrun
	pkl_name = 'saved/lg_' + lg_tag + '.pkl'	
	f_out='plot_'+snap_name+'_'+lg_tag+'_LG.png'
	print 'N = ', lg_tag, ' - ', f_out

	if os.path.isfile(pkl_name):
		pkl_file = open(pkl_name, 'r')

		this_lg = pickle.load(pkl_file)
		base_dir=base_path+'/'+num+'/'+lg_base+'/'+subrun+'/'
		box_center = this_lg.geo_com()
		
		for i in range(0, 3):
			box_center[i] = box_center[i] * 1.e-3
		
		print box_center

		if os.path.isfile(base_dir + snap_name):
			print base_dir + snap_name
			plot_rho(base_dir + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units)

	else:
		print 'File not found: ', pkl_name

'''
