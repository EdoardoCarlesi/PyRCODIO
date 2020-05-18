'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    test.py: this file is used to test new routines and functionalities
'''

import read_files as rf
import tools as t
import halo_utils as hu
import pandas as pd

file_halo = '/media/edoardo/Elements/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
file_tree = '/media/edoardo/Elements/CLUES/DATA/trees/2048/00_06/00/df_all_ids.csv'

halo = rf.read_ahf_halo(file_halo)
tree = rf.read_csv_tree(file_tree)

id0 = halo['ID'].loc[0]
hh = halo[halo['ID'] == id0]
this_halo = hu.Halo(hh)

this_halo.assign_subhalos(halo)

hs = hu.HaloHistory(10)

hs.halos.append(hh)
hs.trajectory_around_host()

print(tree[tree])

#print(this_halo.info())
#print(this_halo.distance(pos))





