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

'/home/eduardo/CLUES/DATA/CF2P5/670/512/00'
base_path='/home/eduardo/CLUES/DATA/'
#sub_dir='TEST'
sub_dir='CF3'
#extra_tag='CF3k_YH_v1'
#extra_tag='CF3_YH_v2'
#extra_tag='CF3_YH_v4'
#extra_tag='CF3_YH_h76nocmb'
#extra_tag='CF3_YH_h78'
#extra_tag='CF3_RG_1000'
#extra_tag='CF3_RG_1500'
#extra_tag='CF2P5'
extra_tag='CF3_BGc_256'
#extra_tag='CF3_Brent_v0'

#plot_slice = True
plot_slice = False
#find_clusters = False
find_clusters = True

ini_num=5905
end_num=5906

nbins = 50
#f_rescale = 1
f_rescale = 1./32.
n_files = 8

#box = '400'; num = '256'; box_size=400.0e+3
#box = '500'; num = '512'; box_size=500.0e+3
#box = '500'; num = '512'; box_size=500.0
#box = '670'; num = '512'; box_size=670.0e+3
box = '500'; num = '256'; box_size=500.0e+3
#box = '600'; num = '256'; box_size=600.0e+3

kpc=1.0e+3
#plot_side = 150.0e+3; thickn = 5000.0; units = 'kpc'
#plot_side = 100.0; thickn = 10.0; units = 'Mpc'
#box_size = 670.0; plot_side = 335.0; thickn = 7.5; units = 'kpc'
box_size = 500.0 * kpc; plot_side = 250.0 * kpc; thickn = 6.5 * kpc; units = 'kpc'

#ahf_name='snapshot_054.AHF_halos'; snap_name='snapshot_054'
#ahf_name='snapshot_127.AHF_halos'; snap_name='snapshot_127'
ahf_name='snapshot_019.AHF_halos'; snap_name='snapshot_019'
#ahf_name='snap_019.AHF_halos'; snap_name='snapshot_019'

print('# Run  M(1.e+14)  d(True)  SGx, SGy, SGz	 Vx, Vy, Vz  Lx, Ly, Lz ')

for run_num in range(ini_num, end_num):

	#f_out='plot_'+snap_name+'_'+extra_tag+'_'+str(run_num)+'_bins'+str(nbins)+'.png'
	f_out='plot_'+snap_name+'_'+extra_tag+'_'+str(run_num)+'bins'+str(nbins)+'.png'

	print('N = ', run_num, ' - ', f_out)
	#print ' ---- ', run_num, ' ----- '

	run_num=str(run_num)+'0'
	#base_dir=base_path+sub_dir+'/'+box+'/'+num+'/'+extra_tag+'/'+run_num+'/'
	#base_dir=base_path+sub_dir+'/'+box+'/'+extra_tag+'/'+run_num+'_00/'
	base_dir=base_path+extra_tag+'/'+box+'/'+num+'/'+run_num+'/'

	box_center = [0.5 * box_size, 0.5 * box_size, 0.5 * box_size]
	
	if os.path.isfile(base_dir + ahf_name) and find_clusters == True:
		#print base_dir + ahf_name
		all_halos = read_ahf(base_dir + ahf_name)
		locate_clusters(all_halos, box_center, run_num)
	else:
		print(base_dir + ahf_name, ' not found.')
	
	#if os.path.isfile(base_dir + snap_name + '.0') and plot_slice == True:
	if os.path.isfile(base_dir + snap_name + '.0') and plot_slice == True:
		print(base_dir + snap_name)
		plot_rho(base_dir + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units, n_files)
#        else:
#            print(base_dir + snap_name, ' not found.')

