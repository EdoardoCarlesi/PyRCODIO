from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.mtree import *
from libcosmo.units import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.particles import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *
from config import *
import random
import time
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



file_single='lcdm_grid64_box100_'
file_suffix='.000064.Vweb-ascii'
box_size = 100000.0
grid_size = 64
base_path = '/home/eduardo/CLUES/DATA/FullBox/vweb/'

sub_ini = 0
sub_end = 1

# Save all the selected merger trees to these files
in_lgs = 'saved/rand_web_lgs_'
out_vweb = 'saved/rand_vweb.pkl'

# Append all trees here - use ALL the DB!	
vweb_ev = []

# Loop on the different box realisations
for i_sub in range(sub_ini, sub_end):
    sub_path = '%02d' % i_sub	

    in_vweb = base_path + file_single + sub_path + file_suffix

    print('Reading file', in_vweb)
    vweb = read_vweb(in_vweb, grid_size, box_size)

    # Load the pre-selected 
    in_lgs_web = in_lgs + sub_path + '.pkl' 
    f_lgs_web = open(in_lgs_web, 'rb')
    all_lgs_web = pickle.load(f_lgs_web)
    n_lgs = len(all_lgs_web)

    indexes = [1, 2, 3]
    #print(vweb.evals[0, indexes], vweb.evals[1, indexes], vweb.evals[2, indexes])   
    #print(vweb.evals[0, 1, 2, 3] )

    #norm = 4096.0
    norm = 512.0
    #norm = 128.0
    ev1 = [];     ev2 = [];     ev3 = []
    print('Merger tree stats for %d pairs, run = %s .' % (n_lgs, sub_path))

    # Read the vweb file
    for lg in all_lgs_web: 
        this_com = lg.get_com()
        indexes = np.zeros((3), dtype=int)

        # Convert to indexes
        for ix in range(0, 3):
            indexes[ix] = int((grid_size-1) * this_com[ix] / box_size) #* grid_size
            
        these_evals = vweb.evals[:, indexes[0], indexes[1], indexes[2]]   
        vweb_ev.append(these_evals)
        ev1.append(these_evals[0]/norm)
        ev2.append(these_evals[1]/norm)
        ev3.append(these_evals[2]/norm)
        #print(vweb.evals[0, indexes], vweb.evals[1, indexes], vweb.evals[2, indexes])   
        #print(indexes, this_com, these_evals)
        #print(these_evals)

#print(ev1)    
print(np.median(ev1))

col_std = 'blue'
n_bins = 15
#plt.hist(ev1, n_bins, facecolor=col_std, normed=1, alpha=0.5)
plt.hist(ev1, n_bins, facecolor=col_std, alpha=0.75, density=True)
#plt.hist(ev1, n_bins, facecolor=col_std, alpha=0.75, density=True)
#plt.hist(ev2, n_bins, facecolor=col_std, alpha=0.5)
#plt.hist(ev3, n_bins, facecolor=col_std, alpha=0.5)
fig_ev1 = 'test_ev1'

#print(len(ev1))

ax = plt.gca()
ax.set_xlabel('$\lambda_1$')
#ax.set_xlabel('t [Gyr]')
ax.set_ylabel('$n$')
#ax.axis([x_min, x_max, y_min, y_max])
#ax.set_ylabel('$n( < \\tau _F)$')
#plt.title('Formation time MW')
#plt.yticks(ticks)
plt.tight_layout()
plt.savefig(fig_ev1)
plt.clf()
plt.cla()
