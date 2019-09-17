from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.mtree import *
from libcosmo.units import *
from libcosmo.halos import *
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
sub_end = 5

n_steps = 54
min_tree_steps = 50

## Save a list of valid LGs found
out_select_lgs = 'saved/rand_select_lgs.pkl'

# Save all the selected merger trees to these files
out_mws = 'saved/rand_mws_all_mah.pkl'
out_m31s = 'saved/rand_m31_all_mah.pkl'

# Save output statistics
out_stat_mah = 'saved/rand_out_stat_mah.pkl'
out_stat_time = 'saved/rand_out_stat_time.pkl'

# Skip the reading part from the DB
if do_trees_db == False:
	sub_ini = 0
	sub_end = 0

# Append all trees here - use ALL the DB!	
trees_mw = []; trees_m31 = []
iValid = 0; valid_lgs = []

# Loop on the different box realisations
for i_sub in range(sub_ini, sub_end):
    sub_path = '%02d' % i_sub	

    #in_db = base_path + 'fullbox_' + sub_path + '_trees.db'
    in_db = base_path + 'trees/fullbox_' + sub_path + '_trees.db'

    # Load the pre-selected 
    out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
    f_out_lgs = open(out_lgs, 'rb')
    all_lgs = pickle.load(f_out_lgs)
    n_lgs = len(all_lgs)

    print('Merger tree stats for %d pairs, run = %s .' % (n_lgs, sub_path))
    newSql = SQL_IO(in_db, n_steps)
    columnReadID = 'allHaloIDs'
    columnReadPT = 'allNumPart'

    for this_lg in all_lgs:
        this_tree_mw = newSql.select_tree(this_lg.LG1.ID, columnReadPT)
        this_tree_m31 = newSql.select_tree(this_lg.LG2.ID, columnReadPT)
        this_ids_mw = newSql.select_tree(this_lg.LG1.ID, columnReadID)
        this_ids_m31 = newSql.select_tree(this_lg.LG2.ID, columnReadID)

        valid_tree_mw = np.where(this_tree_mw > 0)
        valid_tree_m31 = np.where(this_tree_m31 > 0)

        if valid_tree_mw[0].size > min_tree_steps and valid_tree_m31[0].size > min_tree_steps:
            this_mtree_mw = MergerTree(n_steps, this_tree_mw, this_ids_mw, sub_path)
            this_mtree_m31 = MergerTree(n_steps, this_tree_m31, this_ids_m31, sub_path)
            trees_mw.append(this_mtree_mw)
            trees_m31.append(this_mtree_m31)
            valid_lgs.append(this_lg)

            iValid += 1 

    print("ValidLgs so far: ", len(valid_lgs), iValid)

# Save the trees!
if do_trees_db == True:
    print('Found %d valid LG trees.' % (iValid))
    print('Saving MW trees to %s and M31 trees to %s and selected lgs to : %s ' % (out_mws, out_m31s, out_lgs))
    f_out_mws = open(out_mws, 'wb')
    pickle.dump(trees_mw, f_out_mws)
    f_out_m31s = open(out_m31s, 'wb')
    pickle.dump(trees_m31, f_out_m31s)
    f_out_select_lgs = open(out_select_lgs, 'wb')
    pickle.dump(valid_lgs, f_out_select_lgs)

    try:
        newSql.close()
    except:
        'Do nothing'
    
    all_valid_lgs = valid_lgs

else:
    print('Loading MW trees from %s and M31 trees from %s.' % (out_mws, out_m31s))
    f_out_mws  = open(out_mws, 'rb')
    trees_mw   = pickle.load(f_out_mws)
    f_out_m31s = open(out_m31s, 'rb')
    trees_m31  = pickle.load(f_out_m31s)
    f_out_select_lgs = open(out_select_lgs, 'rb')
    all_valid_lgs = pickle.load(f_out_select_lgs)
    
iAll = len(all_valid_lgs)
print('Found %d trees.' % iAll)

        ##########################################################################################
        #                                                                                        #
        # Do some statistics now --------> INCLUDE SELECTION ON THE LOCAL GROUP MODEL !!!!!      #
        #                                                                                        #
        ##########################################################################################

boolSmooth = False
#boolSmooth = True
boolSelect = True

mass_min = 0.5e+12
mass_max = 10.0e+12
ratio_max = 10.0
ratio_min = 1.0
dist_min = 300.0
dist_max = 1500.0

valid_ind = []; ind = 0; iValid = 0
if boolSelect == True:

    for this_lg in all_valid_lgs:
        
        if this_lg.m_tot() > mass_min and this_lg.m_ratio() > ratio_min and this_lg.r_halos() > dist_min and \
                this_lg.m_tot() < mass_max and this_lg.m_ratio() < ratio_max and this_lg.r_halos() < dist_max:
        
            mmt_mw = trees_mw[ind].last_major_merger(boolSmooth) * time_step
            mmt_m31 = trees_m31[ind].last_major_merger(boolSmooth) * time_step
            mmt_mw = trees_mw[ind].last_major_merger(boolSmooth) * time_step
            mmt_m31 = trees_m31[ind].last_major_merger(boolSmooth) * time_step
            #threshold = random.randrange(0.0, 6.0) + 4.0 + 0.5 * random.randrange(0.0, 6.0) 
            valid_ind.append(ind)
            iValid += 1
            
            '''
            random.seed()
            thr0 = 6.0
            thr1 = random.uniform(0.0,13.0) 
    
            if mmt_mw < thr0 and mmt_m31 < thr0:
                valid_ind.append(ind)
                iValid += 1
            elif mmt_mw < thr1 and mmt_m31 < thr1:
                valid_ind.append(ind)
                iValid += 1
            '''

        ind += 1

else:
    iValid = len(all_valid_lgs)
    for i in range(0, iValid):
        valid_ind.append(i)

print("Using %d pairs out of %d." % (iValid, iAll))

time_stats = np.zeros((2, 2, iValid))	# LMMT & FT
mah_stats = np.zeros((2, iValid, n_steps))
mah_avgs = np.zeros((2, n_steps, 3))	# This contains median and percentiles 

all_ft = []
all_mmt = []
thr0 = 7.0

for i_v in range(0, iValid):
    i_t = valid_ind[i_v]
    time_stats[0][1][i_v] = trees_mw[i_t].last_major_merger(boolSmooth) * time_step
    time_stats[0][0][i_v] = trees_mw[i_t].formation_time(boolSmooth) * time_step
    time_stats[1][1][i_v] = trees_m31[i_t].last_major_merger(boolSmooth) * time_step
    time_stats[1][0][i_v] = trees_m31[i_t].formation_time(boolSmooth) * time_step

    mah_stats[0, i_v] = trees_mw[i_t].nPartNorm 
    mah_stats[1, i_v] = trees_m31[i_t].nPartNorm

    '''
    thr0 = 7.0
    if time_stats[0][0][i_v] > thr0:
        random.seed()
        add1 = random.uniform(0.0,4.0)
        time_stats[0][0][i_v] += add1
        
    if time_stats[1][0][i_v] > thr0:
    '''

for i_n in range(0, n_steps):
    percent = [32, 50, 68]
#    percent = [15, 35, 55]

    for i_lg in range(0, 2):
        nonzero = np.where(mah_stats[i_lg, :, i_n] > 0)
        mah_avgs[i_lg, i_n, :] = np.percentile(mah_stats[i_lg, nonzero, i_n], percent)
	
# Save the statistics
print('Saving MAH random LCDM statistics to %s and %s' % (out_stat_mah, out_stat_time))
f_mah = open(out_stat_mah, 'wb')
f_time = open(out_stat_time, 'wb')

pickle.dump(time_stats, f_time)
pickle.dump(mah_avgs, f_mah)
