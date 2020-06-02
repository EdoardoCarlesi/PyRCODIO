'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_sub_mass_function.py: wrapper for utilites to compute and plot subhalo mass functions
'''

import read_files as rf
import tools as t
import halo_utils as hu
import pandas as pd
import seaborn as sns
import pickle as pkl

base_ahf = '/media/edoardo/data1/DATA/'

file_ahf = '/media/edoardo/data1/DATA/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
file_mah = '/media/edoardo/data1/DATA/00_06/00/halo_'
form_mah = '.allinfo'
data_path = '/home/edoardo/CLUES/PyRCODIO/data/'

time = rf.read_time(data_path)
print(time)

all_halo_mah = []

mahs = rf.read_mah_halo(file_ahf, file_mah, time)



#print(mahs[0].m_t())
#print(mahs[0].m_max())
#print(mahs[0].x_t())
print(mahs[0].formation_time())
#print(mahs[0].m_t_norm())
#print(mahs[0].last_major_merger())

#print(mahs[0].head())
#print(mahs[0]['ID'])
#cols = str(mahs[0].columns)
#row = str(mahs[0].iloc[[1]])
#print(row.split())
#print(cols.split())
#mahs[0].columns = cols.split()
#print(mahs[0].head()) 
#print(mahs[0]['#ID(1)']) 
#print(mahs[0]['Mvir(4)']) 
#print(mahs[0]['HostHalo(2)']) 

"""
    TEST READ IN ROUTINES FOR MAH and AHF

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
"""




