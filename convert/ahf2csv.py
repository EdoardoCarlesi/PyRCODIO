'''
    Python Routines for COsmology and Data I/ (PyRCODIO) v0.2
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    ahf2csv.py: convert (and compress) AHF halo catalogs to csv files
'''

import pandas as pd
import sys
sys.path.insert(1, '/home/edoardo/CLUES/PyRCODIO/')
import read_files as rf

this_ahf = sys.argv[1]
mpi = sys.argv[2]

out_file = this_ahf + '.csv'

halo_df = rf.read_ahf_halo(this_ahf, file_mpi=mpi)
halo_df.to_csv(out_file)




