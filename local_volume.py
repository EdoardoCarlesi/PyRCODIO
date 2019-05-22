from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halos import *
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
sub_dir='2048'
snap_name='snapshot_054'

all_lg_base = simu_runs()

plot_slice = True
#plot_slice = False

ini_num=0
end_num=10

nbins = 1000
f_rescale = 1.0 / 1.0

box = '100'; num = '256'; box_size=100.0
#plot_side = 25.0e+3; thickn = 5000.0; units = 'kpc'
plot_side = 15.0; thickn = 4.0; units = 'Mpc'

for lg_base in all_lg_base:

	for run_num in range(ini_num, end_num):
	
		f_out='plot_'+snap_name+'_'+lg_base+'_0'+str(run_num)+'_bins'+str(nbins)+'.png'
	
		base_dir=base_path+sub_dir+'/'+lg_base+'/0'+str(run_num)+'/'

		box_center = [0.5 * box_size, 0.5 * box_size, 0.5 * box_size]
	
		print base_dir + snap_name

		if os.path.isfile(base_dir + snap_name) and plot_slice == True:
			print base_dir + snap_name
			plot_rho(base_dir + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units)



