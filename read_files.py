'''
    Python Routines for Cosmology and Data I/O (PyRCODIO) 
    Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com
'''

import pandas as pd
import numpy as np

"""
    This function assumes that the halo catalog format is AHF and is correcting accordingly.
    This can be of course very easily modified
"""
def read_ahf_halo(file_name):
    halo = pd.read_csv(file_name, sep='\t')
    halo.shift(1, axis=1)

    # Rearrange the columns, the first one is being read incorrectly so we need to split
    halo['ID'] = halo['#ID(1)'].apply(lambda x: x.split()[0])
    halo['HostHalo'] = halo['#ID(1)'].apply(lambda x: x.split()[1])

    # There's a NaN unknown column being read at the end of the file
    new_columns = halo.columns[2:].insert(len(halo.columns[2:-2]), '0')

    # Now drop some useless stuff
    halo.drop('#ID(1)', axis=1, inplace=True)
    halo.columns = new_columns
    halo.drop('0', axis=1, inplace=True)

    return halo


"""
    This is sort of useless but allows to keep some simmetry in the program reading routines
"""
def read_csv_tree(file_name):
    tree = pd.read_csv(file_name)

    return tree

