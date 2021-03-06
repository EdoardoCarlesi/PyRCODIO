#!/usr/bin/python

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
import time
import pickle
from libcosmo.grid import *
from libSQL.sqllib import *
from libSQL.mtree import *
import pandas as pd
import subprocess as sbp
import os

box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/LGF/1024/'
list_tmp = base_path + 'dirs.tmp'
#this_db='/home/eduardo/CLUES/DATA/LGF/trees/lgf_fb_trees.db'
#this_db='/home/eduardo/CLUES/DATA/LGF/trees/lgf_all_trees.db'
this_db='/home/eduardo/CLUES/DATA/LGF/trees/lgf_n500_trees.db'

sub_ini = 0
sub_end = 1

tot_files = 1
use_files = 1
n_steps = 54
n_lgs = 1

time_step = 0.25	# GYrs
all_trees = []
min_tree_size = 40

list_sub_sh = "cd " + base_path + "; ls -d ??_?? > " + list_tmp
os.system(list_sub_sh)
paths_file = open(list_tmp, 'r')
paths = paths_file.readlines()

ip = 0

columnReadID = 'allHaloIDs'
columnReadPT = 'allNumPart'

out_trees_mw = 'saved/all_trees_lgf_mw.pkl'
out_trees_m31 = 'saved/all_trees_lgf_m31.pkl'
out_stat_mah = 'saved/all_lgf_stat_mah.pkl'
out_stat_time = 'saved/all_lgf_stat_time.pkl'
out_valid_trees = 'saved/all_lgf_valid_trees.pkl'

valid_lgs = []
trees_mw_cs = []
trees_m31_cs = []

#read_db = True
read_db = False

if read_db == False:
	paths = []
else:
	# Load the SQL database containing all the LGF trees
	newSql = SQL_IO(this_db, n_steps)

iValid = 0; iBroken = 0;

for path in paths:
	sub_path = path.replace("\n", "")
	ip += 1
	
	out_lgs = 'saved/lgs_' + sub_path + '.pkl'

	if os.path.isfile(out_lgs):
		out_size = os.path.getsize(out_lgs)

		if out_size > 0:
			f_out_lgs = open(out_lgs, 'rb')
			all_lgs = pickle.load(f_out_lgs)
			n_lgs = len(all_lgs)

			#if n_lgs > 0:
			#	print("LGF subpath %s" % (sub_path))

	#print('Merger tree stats for %d pairs, run = %s' % (n_lgs, sub_path))

	for this_lg in all_lgs:
		this_mw = newSql.select_tree(this_lg.LG1.ID, columnReadPT)
		this_m31 = newSql.select_tree(this_lg.LG2.ID, columnReadPT)
		ids_mw = newSql.select_tree(this_lg.LG2.ID, columnReadID)
		ids_m31 = newSql.select_tree(this_lg.LG2.ID, columnReadID) 

		valid_mw = np.where(this_mw > 0)
		valid_m31 = np.where(this_m31 > 0)

		if valid_mw[0].size > min_tree_size and valid_m31[0].size > min_tree_size:
                        this_mtree_mw = MergerTree(n_steps, this_mw, ids_mw, sub_path)
                        this_mtree_m31 = MergerTree(n_steps, this_m31, ids_m31, sub_path)
                        trees_mw_cs.append(this_mtree_mw)
                        trees_m31_cs.append(this_mtree_m31)
                        valid_lgs.append(this_lg)
                        iValid += 1 
try:
	newSql.close()

except:
	'Do nothing'

if read_db == True:
        f_trees_mw = open(out_trees_mw, 'wb')
        f_trees_m31 = open(out_trees_m31, 'wb')
        f_valid_trees = open(out_valid_trees, 'wb')
        pickle.dump(trees_mw_cs, f_trees_mw)
        pickle.dump(trees_m31_cs, f_trees_m31)
        pickle.dump(valid_lgs, f_valid_trees)
        all_valid_lgs = valid_lgs
else:
        f_trees_mw = open(out_trees_mw, 'rb')
        f_trees_m31 = open(out_trees_m31, 'rb')
        f_valid_trees = open(out_valid_trees, 'rb')
        trees_mw_cs = pickle.load(f_trees_mw)
        trees_m31_cs = pickle.load(f_trees_m31)
        all_valid_lgs = pickle.load(f_valid_trees)

iAll = len(all_valid_lgs)
print('Loaded %d valid trees.' % (iAll))

    #
    #  DO SOME SELECTION ON THE LG PROPERTIES!!!! 
    #   

#boolSmooth = False
boolSmooth = True
boolSelect = True

mass_min = 0.5e+12
mass_max = 4.0e+12
ratio_max = 3.5
ratio_min = 1.0
dist_min = 300.0
dist_max = 1200.0

valid_ind = []; ind = 0
if boolSelect == True:

    for this_lg in all_valid_lgs:
        if this_lg.m_tot() > mass_min and this_lg.m_ratio() > ratio_min and this_lg.r_halos() > dist_min and \
                this_lg.m_tot() < mass_max and this_lg.m_ratio() < ratio_max and this_lg.r_halos() < dist_max:
        
            valid_ind.append(ind)
            iValid += 1

        ind += 1

else:
    iValid = len(all_valid_lgs)
    for i in range(0, iValid):
        valid_ind.append(i)

print("Using %d pairs out of %d." % (iValid, iAll))

# Do some statistics now
time_stats = np.zeros((2, 2, iValid))	# LMMT & FT
mah_stats = np.zeros((2, iValid, n_steps))
mah_avgs = np.zeros((2, n_steps, 3))	# This contains median and percentiles 

all_ft = []
all_mmt = []

for i_v in range(0, iValid):
    i_t = valid_ind[i_v]
    time_stats[0][0][i_v] = (trees_mw_cs[i_t].formation_time(boolSmooth) -1) * time_step
    time_stats[0][1][i_v] = (trees_mw_cs[i_t].last_major_merger(boolSmooth) ) * time_step
    time_stats[1][0][i_v] = (trees_m31_cs[i_t].formation_time(boolSmooth) -1) * time_step
    time_stats[1][1][i_v] = (trees_m31_cs[i_t].last_major_merger(boolSmooth) ) * time_step
    mah_stats[0, i_v] = trees_mw_cs[i_t].nPartNorm 
    mah_stats[1, i_v] = trees_m31_cs[i_t].nPartNorm
#    print(trees_mw_cs[i_t].smooth_tree())

    this_mmt = trees_mw_cs[i_t].formation_time(boolSmooth)
    this_ft = trees_mw_cs[i_t].last_major_merger(boolSmooth)
#    print(trees_mw_cs[i_t].nPartNorm, trees_mw_cs[i_t].nPartNorm[this_mmt], trees_mw_cs[i_t].nPartNorm[this_ft])
    all_ft.append(this_ft * time_step)
    all_mmt.append(this_mmt * time_step)



#print(np.percentile(all_ft, [25, 50, 75]))


for i_n in range(0, n_steps):
	percent = [35, 55, 75]
	#percent = [32, 50, 68]
	
	for i_lg in range(0, 2):
		nonzero = np.where(mah_stats[i_lg, :, i_n] > 0)
		mah_avgs[i_lg, i_n, :] = np.percentile(mah_stats[i_lg, nonzero, i_n], percent)
	
# Save the "condensed" data
print('Saving MAH random LCDM statistics to %s and %s' % (out_stat_mah, out_stat_time))
f_mah = open(out_stat_mah, 'wb')
f_time = open(out_stat_time, 'wb')

pickle.dump(time_stats, f_time)
pickle.dump(mah_avgs, f_mah)

