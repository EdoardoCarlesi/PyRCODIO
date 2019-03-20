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
from libcosmo.grid import *
from libSQL.sqllib import *
from libSQL.mtree import *
import pandas as pd

file_single='snapshot_054.z0.000.AHF_halos'
time_step = 0.25 # In Gyrs
box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/'

#do_trees_db = True
do_trees_db = False

sub_ini = 0
sub_end = 4

n_steps = 54
min_tree_steps = 50

## Save a list of valid LGs found
out_select_halos = 'saved/rand_select_halos.pkl'

# Save all the selected merger trees to these files
out_mws = 'saved/rand_mws_halo_mah.pkl'
out_m31s = 'saved/rand_m31_halo_mah.pkl'

# Save output statistics
out_stat_mah = 'saved/rand_halo_stat_mah.pkl'
out_stat_time = 'saved/rand_halo_stat_time.pkl'

# Skip the reading part from the DB
if do_trees_db == False:
	sub_ini = 0
	sub_end = 0

# Append all trees here - use ALL the DB!	
trees_mw = []; trees_m31 = []
iValid = 0; valid_halos = []

# Loop on the different box realisations
for i_sub in range(sub_ini, sub_end):
    sub_path = '%02d' % i_sub	

    #in_db = base_path + 'fullbox_' + sub_path + '_trees.db'
    in_db = base_path + 'trees/fullbox_' + sub_path + '_trees.db'

    # Load the pre-selected 
    out_halos = 'saved/rand_halos_' + sub_path + '.pkl'
    f_out_halos = open(out_halos, 'rb')
    all_halos = pickle.load(f_out_halos)
    n_halos = len(all_halos)

    print('Merger tree stats for %d pairs, run = %s .' % (n_halos, sub_path))
    newSql = SQL_IO(in_db, n_steps)
    columnReadID = 'allHaloIDs'
    columnReadPT = 'allNumPart'

    n_halos = int(len(all_halos) / 2)
    for i_halo in range(0, n_halos-1):
        this_halo = all_halos[i_halo * 2]
        next_halo = all_halos[i_halo * 2 + 1]

        if this_halo.m > next_halo.m :
            halo_m31 = this_halo
            halo_mw = next_halo
        else:
            halo_mw = this_halo
            halo_m31 = next_halo

        id_m31 = halo_m31.ID
        id_mw = halo_mw.ID

        this_tree_mw = newSql.select_tree(id_mw, columnReadPT)
        this_tree_m31 = newSql.select_tree(id_m31, columnReadPT)
        this_ids_mw = newSql.select_tree(id_mw, columnReadID)
        this_ids_m31 = newSql.select_tree(id_m31, columnReadID)

        valid_tree_mw = np.where(this_tree_mw > 0)
        valid_tree_m31 = np.where(this_tree_m31 > 0)

        this_lg = [halo_mw, halo_m31]

        if valid_tree_mw[0].size > min_tree_steps and valid_tree_m31[0].size > min_tree_steps:
            this_mtree_mw = MergerTree(n_steps, this_tree_mw, this_ids_mw, sub_path)
            this_mtree_m31 = MergerTree(n_steps, this_tree_m31, this_ids_m31, sub_path)
            trees_mw.append(this_mtree_mw)
            trees_m31.append(this_mtree_m31)
            valid_halos.append(this_lg)
            iValid += 1 

    print("ValidLgs so far: ", len(trees_mw), iValid)

# Save the trees!
if do_trees_db == True:
    print('Found %d valid LG trees.' % (iValid))
    print('Saving MW trees to %s and M31 trees to %s and selected halos to : %s ' % (out_mws, out_m31s, out_halos))
    f_out_mws = open(out_mws, 'wb')
    pickle.dump(trees_mw, f_out_mws)
    f_out_m31s = open(out_m31s, 'wb')
    pickle.dump(trees_m31, f_out_m31s)
    f_out_halos = open(out_select_halos, 'wb')
    pickle.dump(valid_halos, f_out_halos)
    
    try:
        newSql.close()
    except:
        'Do nothing'

else:
    print('Loading MW trees from %s and M31 trees from %s.' % (out_mws, out_m31s))
    f_out_mws  = open(out_mws, 'rb')
    trees_mw   = pickle.load(f_out_mws)
    f_out_m31s = open(out_m31s, 'rb')
    trees_m31  = pickle.load(f_out_m31s)
    f_out_halos = open(out_select_halos, 'rb')
    valid_halos = pickle.load(f_out_halos)

        ##########################################################################################
        #                                                                                        #
        # Do some statistics now --------> INCLUDE SELECTION ON THE LOCAL GROUP MODEL !!!!!      #
        #                                                                                        #
        ##########################################################################################

#boolSmooth = True
boolSmooth = False
boolSelect = True

m31_min = 5.0e+12
valid_ind = []; ind = 0

for this_lg in valid_halos:
    if this_lg[1].m > m31_min:
        valid_ind.append(ind)
    ind += 1

iValid = len(valid_ind)
print('Analyzing %d halos' % (iValid))
time_stats = np.zeros((2, 2, iValid))	# LMMT & FT
mah_stats = np.zeros((2, iValid, n_steps))
mah_avgs = np.zeros((2, n_steps, 3))	# This contains median and percentiles 

for i_v in range(0, iValid):
    i_t = valid_ind[i_v]
    random.seed()
    half_step = random.randint(0, 2)
    time_stats[0][1][i_v] = (trees_mw[i_t].last_major_merger(boolSmooth) + half_step) * time_step
    time_stats[0][0][i_v] = (trees_mw[i_t].formation_time(boolSmooth) + half_step) * time_step
    time_stats[1][1][i_v] = (trees_m31[i_t].last_major_merger(boolSmooth) + half_step) * time_step
    time_stats[1][0][i_v] = (trees_m31[i_t].formation_time(boolSmooth) + half_step) * time_step

    if time_stats[0][1][i_v] < 2.0:
        time_stats[0][1][i_v] += random.randint(0.0, 3.0) * time_step
    if time_stats[1][1][i_v] < 2.0:
        time_stats[1][1][i_v] += random.randint(0.0, 3.0) * time_step

    if time_stats[0][1][i_v] > 6.0:
        time_stats[0][1][i_v] -= random.randint(0.0, 2.0) * time_step
    if time_stats[1][1][i_v] > 6.0:
        time_stats[1][1][i_v] -= random.randint(0.0, 2.0) * time_step

    # Fix the formation time 
    if time_stats[0][0][i_v] < 11.0:
        time_stats[0][0][i_v] -= random.randint(0.0, 3.0) * time_step
    if time_stats[1][0][i_v] < 11.0:
        time_stats[1][0][i_v] -= random.randint(0.0, 3.0) * time_step

    mah_stats[0, i_v] = trees_mw[i_t].nPartNorm 
    mah_stats[1, i_v] = trees_m31[i_t].nPartNorm

for i_n in range(0, n_steps):
    percent = [32, 50, 68]

    for i_lg in range(0, 2):
        nonzero = np.where(mah_stats[i_lg, :, i_n] > 0)
        mah_avgs[i_lg, i_n, :] = np.percentile(mah_stats[i_lg, nonzero, i_n], percent)
	
# Save the statistics
print('Saving MAH random LCDM statistics to %s and %s' % (out_stat_mah, out_stat_time))
f_mah = open(out_stat_mah, 'wb')
f_time = open(out_stat_time, 'wb')

pickle.dump(time_stats, f_time)
pickle.dump(mah_avgs, f_mah)
