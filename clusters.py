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


base_dir='/home/eduardo/CLUES/DATA/CF2P5/256/00/'
ahf_name='snapshot_054.AHF_halos'
snap_name='snapshot_054'
f_out='plot_256.png'

box_size = 400.0e+3
plot_side = 30.0e+3

box_center = [0.5 * box_size, 0.5 * box_size, 0.5 * box_size]
#all_halos = read_ahf(base_dir + ahf_name)
#locate_clusters(all_halos, box_center)

plot_lv(base_dir + snap_name, box_center, plot_side, f_out)

