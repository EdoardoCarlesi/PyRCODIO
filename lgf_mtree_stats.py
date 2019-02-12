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
import subprocess as sbp
import os

file_single='snapshot_054.z0.000.AHF_halos'
#file_single='snapshot_054.0000.z0.000.AHF_halos'

box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/LGF/1024/'
list_tmp = base_path + 'dirs.tmp'
this_db='../lgf_fb_trees.db'

sub_ini = 0
sub_end = 1

tot_files = 1
use_files = 1
n_steps = 54
n_lgs = 1

all_trees = np.zeros((n_lgs, n_steps))
list_sub_sh = "cd " + base_path + "; ls -d ??_?? > " + list_tmp
os.system(list_sub_sh)
paths_file = open(list_tmp, 'r')
paths = paths_file.readlines()

ip = 0

# Load the SQL database containing all the LGF trees
in_db = base_path + this_db
newSql = SQL_IO(in_db, n_steps)
columnReadID = 'allHaloIDs'
columnReadPT = 'allNumPart'


for path in paths:
	sub_path = path.replace("\n", "")
	ip += 1
	
	out_lgs = 'saved/lgs_' + sub_path + '.pkl'

	if os.path.isfile(out_lgs):
		out_size = os.path.getsize(out_lgs)

		if out_size > 0:
			f_out_lgs = open(out_lgs, 'r')
			all_lgs = pickle.load(f_out_lgs)
			n_lgs = len(all_lgs)

			if n_lgs > 0:
				print("LGF subpath %s" % (sub_path))


	# Load the pre-selected 

	print('Merger tree stats for %d pairs, run = %s .' % (n_lgs, sub_path))
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
			print(this_mtree.last_major_merger(True), this_mtree.last_major_merger(False))
			#print(this_mtree.formation_time(False), this_mtree.formation_time(False))
		else:
			iBroken += 1

		
print('Found %d valid, %d broken trees.' % (iValid, iBroken))

newSql.close()
'''
