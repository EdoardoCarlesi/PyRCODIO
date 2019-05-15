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

#do_read_vweb = True
do_read_vweb = False
n_bins = 50

file_single='lcdm_grid64_box100_'
file_suffix='.000064.Vweb-ascii'
box_size = 100000.0
grid_size = 64
base_path = '/home/eduardo/CLUES/DATA/FullBox/vweb/'

sub_ini = 0
sub_end = 5

# Save all the selected merger trees to these files
out_vweb = 'saved/rand_vweb.pkl'
out_ev1 = 'saved/rand_ev1.pkl'
out_ev2 = 'saved/rand_ev2.pkl'
out_ev3 = 'saved/rand_ev3.pkl'
in_lgs = 'saved/rand_web_lgs_'

# Append all trees here - use ALL the DB!	
vweb_ev = []
ev1 = [];     ev2 = [];     ev3 = []

# Loop on the different box realisations
if do_read_vweb == True:
    for i_sub in range(sub_ini, sub_end):
        sub_path = '%02d' % i_sub	

        in_vweb = base_path + file_single + sub_path + file_suffix

        # Load the pre-selected 
        in_lgs_web = in_lgs + sub_path + '.pkl' 
        f_lgs_web = open(in_lgs_web, 'rb')
        all_lgs_web = pickle.load(f_lgs_web)

        print('Reading file', in_vweb)
        vweb = read_vweb(in_vweb, grid_size, box_size)

        n_lgs = len(all_lgs_web)

        indexes = [1, 2, 3]

        #norm = 4096.0
        norm = 512.0
        #norm = 128.0
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

    f_ev1 = open(out_ev1, 'wb'); f_ev2 = open(out_ev2, 'wb'); f_ev3 = open(out_ev3, 'wb')
    pickle.dump(ev1, f_ev1)
    pickle.dump(ev2, f_ev2)
    pickle.dump(ev3, f_ev3)

else:
    f_ev1 = open(out_ev1, 'rb'); f_ev2 = open(out_ev2, 'rb'); f_ev3 = open(out_ev3, 'rb')
    ev1 = pickle.load(f_ev1)
    ev2 = pickle.load(f_ev2)
    ev3 = pickle.load(f_ev3)
    print('Loading', len(ev1), ' eigenvalues.')

col_std = 'blue'
plt.hist(ev1, n_bins, facecolor=col_std, alpha=0.75, density=True)
fig_ev1 = 'dist_ev1'
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


col_std = 'blue'
plt.hist(ev2, n_bins, facecolor=col_std, alpha=0.75, density=True)
fig_ev2 = 'dist_ev2'
ax = plt.gca()
ax.set_xlabel('$\lambda_2$')
#ax.set_xlabel('t [Gyr]')
ax.set_ylabel('$n$')
#ax.axis([x_min, x_max, y_min, y_max])
#ax.set_ylabel('$n( < \\tau _F)$')
#plt.title('Formation time MW')
#plt.yticks(ticks)
plt.tight_layout()
plt.savefig(fig_ev2)
plt.clf()
plt.cla()

col_std = 'blue'
plt.hist(ev3, n_bins, facecolor=col_std, alpha=0.75, density=True)
fig_ev3 = 'dist_ev3'
ax = plt.gca()
ax.set_xlabel('$\lambda_3$')
#ax.set_xlabel('t [Gyr]')
ax.set_ylabel('$n$')
#ax.axis([x_min, x_max, y_min, y_max])
#ax.set_ylabel('$n( < \\tau _F)$')
#plt.title('Formation time MW')
#plt.yticks(ticks)
plt.tight_layout()
plt.savefig(fig_ev3)
plt.clf()
plt.cla()
