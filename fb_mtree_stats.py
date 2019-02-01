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

box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/'
#sub_path = '01'

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

	iValid = 0; iBroken = 0;

	for this_lg in all_lgs:
		this_tree = newSql.select_tree(this_lg.LG1.ID, columnReadPT)
		this_ids = newSql.select_tree(this_lg.LG1.ID, columnReadID)

		valid_tree = np.where(this_tree > 0)

		if valid_tree[0].size > 40:
			iValid +=1 
			this_mtree = MergerTree(n_steps, this_tree)
			#this_mtree.info()
			#if (this_mtree.last_major_merger() != None):
			#	print(this_mtree.last_major_merger())
			#if (this_mtree.formation_time() > 0):
			print(this_mtree.formation_time())
		else:
			iBroken += 1

		

	#if this_tree[0] != 0:
	#print('%lu %s' % (this_lg.LG1.ID, this_ids))
		#print(this_tree)
	#	print(this_ids)

print('Found %d valid, %d broken trees.' % (iValid, iBroken))

#these_trees = newSql.select_trees(testIDs)

#print these_trees

#print(these_trees)
#for this_lg in all_lgs[0:10]:
#	print(this_lg.info())
#	this_tree = newSql.get_full_mtree(testID)
	

#	print(this_tree.norm_mass())
#newSql.cursor.execute('COMMIT')

newSql.close()

