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
from config import *
import time
import pickle
from libcosmo.grid import *


#base_dir='/home/eduardo/CLUES/DATA/CF2P5/512/70_00/'
#base_dir='/home/eduardo/CLUES/DATA/128/CF2P5/'
#base_dir='/home/eduardo/CLUES/DATA/TEST/RAND/'; extra_tag='_std'
base_dir='/home/eduardo/CLUES/DATA/TEST/cf2/'; extra_tag='_js'
#base_dir='/home/eduardo/CLUES/DATA/CF2P5/256b/90/'
#base_dir='/home/eduardo/CLUES/DATA/CF2YH/375494/'
#base_dir='/home/eduardo/CLUES/DATA/CF2YH/289090/'
#base_dir='/home/eduardo/CLUES/DATA/cf2JS/289090/'
#base_dir='/home/eduardo/CLUES/DATA/CF2YH/286393/'
#base_dir='/home/eduardo/CLUES/DATA/CF2YH/288191/'
#base_dir='/home/eduardo/CLUES/DATA/QL/'
ahf_name='snapshot_054.AHF_halos'
snap_name='snapshot_054'
#snap_name='snapshot_019'
#f_out='plot_'+snap_name+'_js.png'
f_out='plot_'+snap_name+extra_tag+'.png'

box_size = 100.0e+3; plot_side = 50.0e+3; thickn = 1000.0; units = 'kpc'
#box_size = 100.0; plot_side = 50.0; thickn = 5.0
nbins = 250
f_rescale = 1

box_center = [0.5 * box_size, 0.5 * box_size, 0.5 * box_size]
#all_halos = read_ahf(base_dir + ahf_name)
#locate_clusters(all_halos, box_center)

plot_lv(base_dir + snap_name, box_center, plot_side, f_out, nbins, f_rescale, thickn, units)

