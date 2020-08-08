'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    rs2csv.py: convert (and compress) RS halo catalogs to csv files
'''

import read_files as rf
import halo_utils as hu
import seaborn as sns
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os
import dask.dataframe as dd
import time

catalog = '/srv/cosmdata/multidark/BigMD_3840_Planck1/ROCKSTAR/catalogs/out_0.list'
catalog_compress = '/z/carlesi/STORE/MultiDark/RockStarCSV/BigMD_3840_Planck1/out_0_csv'
rs_df = rf.read_rs_halo(with_dask=True, read_file=catalog)
with_dask = False

if with_dask == True:
    catalog_compress = '/home/edoardo/CLUES/DATA/RS/out_3_csv'
else:
    catalog_compress = '/home/edoardo/CLUES/DATA/RS/out_3.csv'

time0 = time.time()
print('Reading: ', catalog, ' and compressing to csv file: ', catalog_compress)
rs_df = rf.read_rs_halo(with_dask=with_dask, read_file=catalog)

if with_dask == True:
    print('Total number of halos read: ', rs_df.shape[0].compute())
else:
    print('Total number of halos read: ', rs_df['DescID'].count())

print(rs_df.head())
time1 = time.time()
print('Took: ', '%5.3f' % (time1 - time0), ' seconds')

cols_drop = ['DescID', 'Mvir_all', 'M200b', 'M200c', 'M500c', 'M2500c', 'Xoff', 'Voff', 
            'A[x]', 'A[y]', 'A[z]', 'b_to_a(500c)', 'c_to_a(500c)', 'A[x](500c)', 'A[y](500c)', 'A[z](500c)', 'T/|U|', 'M_pe_Behroozi', 'M_pe_Diemer', 'Halfmass_Radius']

rs_df = rs_df.drop(labels=cols_drop, axis=1) #, inplace=True)

print('Col cut: ')
print(rs_df.head())
time2 = time.time()
print('Took: ', '%5.3f' % (time2 - time1), ' seconds')

m_cut = 4.0e+11

rs_df = rs_df[rs_df['Mvir'] > m_cut]

if with_dask == True:
    print('Total halos: ', rs_df.shape[0].compute())
else:
    print('Total number of halos read: ', rs_df['#ID'].count())

print('Mass cut: ')
print(rs_df.head())
time3 = time.time()
print('Took: ', '%5.3f' % (time3 - time2), ' seconds')

rs_df.to_csv(catalog_compress)

