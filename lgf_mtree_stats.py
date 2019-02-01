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

#file_single='snapshot_054.z0.000.AHF_halos'
file_single='snapshot_054.0000.z0.000.AHF_halos'
#file_single='snapshot_054.0000.z0.001.AHF_halos'

box_size = 100000.0
#base_path = '/home/eduardo/CLUES/DATA/LGF/1024/'; sub_path='00_10mpi'
base_path = '/home/eduardo/CLUES/DATA/2048/00_06/'; sub_path = '00'
root_file = 'snapshot_'
suff_halo = '.z0.000.AHF_halos'
suff_part = '.z0.000.AHF_particles'
this_ahf_file = base_path + sub_path + '/' + file_single

tot_files = 1
use_files = 1

m_max = 2.e+15
m_min = 1.e+11

# Local group selection parameters
iso_radius = 2200.
radius = 1000. 
r_max = 1500.
r_min = 350. 
m_min = 5.e+11  
m_max = 5.0e+12 
ratio_max = 5.0
vrad_max = 0.0

resolution = '1024'
env_type = 'std'
out_dir = 'output/'
n_steps = 54

#in_db = '/home/eduardo/CLUES/DATA/FullBox/trees/fullbox_01_trees.db'
#in_db = '/home/eduardo/CLUES/DATA/FullBox/fullbox_00_trees.db'
#in_db = '/home/eduardo/CLUES/DATA/LGF/lgf_all_trees.db'
in_db = '/home/eduardo/CLUES/DATA/LGF/lgf_fb_trees.db'

#out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
#f_out_lgs = open(out_lgs, 'r')
#all_lgs = pickle.load(f_out_lgs)
#testID = all_lgs[0].LG1.ID

print(this_ahf_file)

halos = read_ahf(this_ahf_file)
#print(halos)
newSql = SQL_IO(in_db, n_steps)

#newSql.cursor.execute('BEGIN TRANSACTION')

columnReadID = 'allHaloIDs'
columnReadPT = 'allNumPart'

incompl = 0

#for this_lg in all_lgs[0:1000]:
for halo in halos:
#for halo in halos[0:10]:
	thisID = halo.ID
	this_tree = newSql.select_tree(thisID, columnReadPT)
	this_ids = newSql.select_tree(thisID, columnReadID)
	nInc = np.where(this_tree < 10)

	thisSize = nInc[0].size

	if thisSize > 5 and this_tree[0] > 100:
		print thisID, this_tree[0:10]#, this_ids[0:10]
		incompl += 1

	#if (nInc[0].size > 100):
		#print nInc[0].size	
		#print(len(nInc))

print('Total incomplete: ', incompl, ' out of ', len(halos))

#	print halo.m_dm
	#if this_tree[0] > 300:
#	print(thisID)
#	print(this_tree)
#		print(this_ids)
#	halo.info()
'''
	if halo.m_dm > 1.e+11:
	#testIDs.append(this_lg.LG2.ID)
		this_tree = newSql.select_tree(this_ID, columnReadPT)
		this_ids = newSql.select_tree(this_ID, columnReadID)
'''
	#if this_tree[0] != 0:
	#print('%lu %s' % (this_lg.LG1.ID, this_ids))
		#print(this_tree)
	#	print(this_ids)
#these_trees = newSql.select_trees(testIDs)
#print these_trees
#print(these_trees)
#for this_lg in all_lgs[0:10]:
#	print(this_lg.info())
#	this_tree = newSql.get_full_mtree(testID)
#	print(this_tree.norm_mass())
#newSql.cursor.execute('COMMIT')

newSql.close()

