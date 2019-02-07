#!/usr/bin/python

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
import time
import pickle
from libcosmo.grid import *
from libSQL.sqllib import *
from libSQL.mtree import *
import pandas as pd

file_single='snapshot_054.z0.000.AHF_halos'
#file_single='snapshot_054.0000.z0.000.AHF_halos'
time_step = 0.25 # In Gyrs
box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/'
#sub_path = '01'

#do_trees_db = True
do_trees_db = False

sub_ini = 0
sub_end = 1

tot_files = 1
use_files = 1
resolution = '1024'
env_type = 'std'
out_dir = 'output/'
n_steps = 54
n_lgs = 1

all_trees = np.zeros((n_lgs, n_steps))

# Append all trees here - use ALL the DB!	
trees_mw = []; trees_m31 = []

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

# Loop on the different box realisations
for i_sub in range(sub_ini, sub_end):
	sub_path = '%02d' % i_sub	

	in_db = base_path + 'fullbox_' + sub_path + '_trees.db'
	#in_db = base_path + 'trees/fullbox_' + sub_path + '_trees.db'

	# Load the pre-selected 
	out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
	f_out_lgs = open(out_lgs, 'r')
	all_lgs = pickle.load(f_out_lgs)
	n_lgs = len(all_lgs)

	print('Merger tree stats for %d pairs, run = %s .' % (n_lgs, sub_path))

	newSql = SQL_IO(in_db, n_steps)
	columnReadID = 'allHaloIDs'
	columnReadPT = 'allNumPart'

	iValid = 0; 

	for this_lg in all_lgs:
		this_tree_mw = newSql.select_tree(this_lg.LG1.ID, columnReadPT)
		this_tree_m31 = newSql.select_tree(this_lg.LG2.ID, columnReadPT)

		valid_tree_mw = np.where(this_tree_mw > 0)
		valid_tree_m31 = np.where(this_tree_m31 > 0)

		if valid_tree_mw[0].size > 40 and valid_tree_m31[0].size > 40:
			this_mtree_mw = MergerTree(n_steps, this_tree_mw)
			this_mtree_m31 = MergerTree(n_steps, this_tree_m31)
			trees_mw.append(this_mtree_mw)
			trees_m31.append(this_mtree_m31)
			iValid += 1 

# Save the trees!
if do_trees_db == True:
	print('Found %d valid LG trees.' % (iValid))
	print('Saving MW trees to %s and M31 trees to %s.' % (out_mws, out_m31s))
	f_out_mws = open(out_mws, 'w')
	pickle.dump(trees_mw, f_out_mws)
	f_out_m31s = open(out_m31s, 'w')
	pickle.dump(trees_m31, f_out_m31s)
	newSql.close()
else:
	print('Loading MW trees from %s and M31 trees from %s.' % (out_mws, out_m31s))
	f_out_mws  = open(out_mws, 'r')
	trees_mw   = pickle.load(f_out_mws)
	f_out_m31s = open(out_m31s, 'r')
	trees_m31  = pickle.load(f_out_m31s)
	iValid = len(trees_m31)
	print('Found %d trees.' % iValid)

# Do some statistics now
time_stats = np.zeros((2, 2, iValid))	# LMMT & FT
mah_stats = np.zeros((2, iValid, n_steps))
mah_avgs = np.zeros((2, n_steps, 3))	# This contains median and percentiles 

for i_t in range(0, iValid):
		time_stats[0][0][i_t] = trees_mw[i_t].formation_time(True) * time_step
		time_stats[0][1][i_t] = trees_mw[i_t].last_major_merger(True) * time_step
		time_stats[1][0][i_t] = trees_m31[i_t].formation_time(True) * time_step
		time_stats[1][1][i_t] = trees_m31[i_t].last_major_merger(True) * time_step

		mah_stats[0, i_t] = trees_mw[i_t].nPartNorm 
		mah_stats[1, i_t] = trees_m31[i_t].nPartNorm

for i_n in range(0, n_steps):
	percent = [32, 50, 68]
	
	for i_lg in range(0, 2):
		nonzero = np.where(mah_stats[i_lg, :, i_n] > 0)
		mah_avgs[i_lg, i_n, :] = np.percentile(mah_stats[i_lg, nonzero, i_n], percent)
		#print(i_n, np.percentile(mah_stats[0, nonzero, i_n], percent))
	
#out_stat_mah = 'saved/rand_out_stat_mah.pkl'
#out_stat_time = 'saved/rand_out_stat_time.pkl'
	
# Save the "condensed" data
print('Saving MAH random LCDM statistics to %s and %s' % (out_stat_mah, out_stat_time))
f_mah = open(out_stat_mah, 'w')
f_time = open(out_stat_time, 'w')

pickle.dump(time_stats, f_time)
pickle.dump(mah_avgs, f_mah)
