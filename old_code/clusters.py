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

#'/home/eduardo/CLUES/DATA/CF2P5/670/512/00'
base_path='/home/eduardo/CLUES/DATA/'
#sub_dir='TEST'
sub_dir='CF3'
#extra_tag='CF3_RG_1500'
extra_tag='CF3_BGc_256'
#extra_tag='CF3_Brent_v0'

ini_num=0
end_num=4

seeds = ['59050', '58151', '57252', '56353', '55454', '54555']

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

#print('# Run  M(1.e+14)  d(True)  SGx, SGy, SGz  Vx, Vy, Vz  Lx, Ly, Lz ')
print('# Run\tM(1.e+14)\td(True)\t\tSGx,\t\tSGy,\t\tSGz')

for i_num in range(ini_num, end_num):

    run_num = seeds[i_num]
    base_dir=base_path+sub_dir+'/'+box+'/'+extra_tag+'/'+run_num+ '/'
    box_center = [0.5 * box_size, 0.5 * box_size, 0.5 * box_size]

    if os.path.isfile(base_dir + ahf_name):
        #print base_dir + ahf_name
        all_halos = read_ahf(base_dir + ahf_name)
        locate_clusters(all_halos, box_center, run_num)
    else:
        print(base_dir + ahf_name, ' not found.')
