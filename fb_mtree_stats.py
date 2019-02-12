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
from libSQL.sqllib import *
from libSQL.mtree import *
import pandas as pd

#file_single='snapshot_054.0000.z0.000.AHF_halos'

tot_files = 1
use_files = 1

n_steps=175

file_ahf='/home/edoardo/CLUES/DATA/SIMULATIONS/LGF/1024/00_06/snapshot_176.0000.z0.000.AHF_halos'
#in_db = '/home/eduardo/CLUES/DATA/FullBox/fullbox_00_trees.db'
in_db = '/home/edoardo/devel/MetroC++/output/lgf_test_trees.db'
#out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
out_lgs = 'saved/test_lgs_00_06.pkl'
f_out_lgs = open(out_lgs, 'r')
all_lgs = pickle.load(f_out_lgs)

testID = all_lgs[0].LG1.ID

newSql = SQL_IO(in_db, n_steps)

#newSql.cursor.execute('BEGIN TRANSACTION')
testIDs = []

columnReadID = 'allHaloIDs'
columnReadPT = 'allNumPart'

#for this_lg in all_lgs[0:1000]:

all_halos = read_ahf(file_ahf)


for halo in all_halos:
	if halo.npart > 5000:
		this_id = newSql.select_tree(halo.ID, columnReadID)
		this_tree = newSql.select_tree(halo.ID, columnReadPT)
		
		print(this_id)

		#for iT in range(0, len(this_id)):
		#	print(iT, this_id[iT], this_tree[iT])


'''
for this_lg in all_lgs:
	testIDs.append(this_lg.LG1.ID)
	this_tree = newSql.select_tree(this_lg.LG1.ID, columnRead)
	print(this_tree.reset_index().values) #/float(this_tree[0]))
these_trees = newSql.select_trees(testIDs)
'''

#print these_trees

#print(these_trees)
#for this_lg in all_lgs[0:10]:
#	print(this_lg.info())
#	this_tree = newSql.get_full_mtree(testID)
	

#	print(this_tree.norm_mass())
#newSql.cursor.execute('COMMIT')

newSql.close()

