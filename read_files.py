'''
    Python Routines for Cosmology and Data I/O (PyRCODIO) 
    Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com
    
    read_files.py: read different file formats and types 
'''

import scipy as sp
import halo_utils as hu
import pandas as pd
import numpy as np

import sys
sys.path.append('pygadgetreader/')
from pygadgetreader import *

"""
    This function assumes that the halo catalog format is AHF and is correcting accordingly.
    This can be of course very easily modified
"""
def read_ahf_halo(file_name, file_mpi=True):
    halo = pd.read_csv(file_name, sep='\t')
    halo.shift(1, axis=1)

    # MPI produced files have a slightly different formatting
    if file_mpi == True:

        # Rearrange the columns, the first one is being read incorrectly so we need to split
        halo['ID'] = halo['#ID(1)'].apply(lambda x: x.split()[0])
        halo['HostHalo'] = halo['#ID(1)'].apply(lambda x: x.split()[1])
        
        # There's a NaN unknown column being read at the end of the file
        new_columns = halo.columns[2:].insert(len(halo.columns[2:-2]), '0')

        # Now drop some useless stuff
        halo.drop('#ID(1)', axis=1, inplace=True)
        halo.columns = new_columns
        halo.drop('0', axis=1, inplace=True)

    else:
        halo.rename(columns={"#ID(1)":"ID", "hostHalo(2)":"HostHalo"}, inplace=True)

    halos = []

    n_rows = halo.shape[0]

    for i in range(0, n_rows):
        #this_h = pd.DataFrame(data=halo.loc[i], columns=new_columns)
        this_h = data=halo.loc[i]
        h = hu.Halo(this_h)
        halos.append(h)

    # Halos is a list of Halo objects, halo is a DataFrame type
    return [halos, halo]


"""
    Read directly the full MAH of every single halo once the AHF file is given.
    We also need to pass the time data about a/z/t that will be read by each object.
"""
def read_mah_halo(id_list, mah_path, time):
    mah_format = '.allinfo'
    all_mah = []
    cols = None
    n_cols = 0
    head_count = 0

    for ID in id_list:
        file_mah = mah_path + str(ID) + mah_format

        try:
            # Read header. There is some issue with delimiters so we need to read this separately. Only one time
            head = pd.read_csv(file_mah, sep='\s+', nrows=0)
            cols = head.columns
            n_cols = len(head.columns)
            head_count = 1
            break
        except:
            'This file does not exist'
            #print('File: ', file_mah, ' does not exist.')

    if head_count == 0:
        print('There are no suitable MAH files in folder ', mah_path)
        print('Please make sure there are enough files with the right format. Exiting...')
        id_list = []

    for ID in id_list:
        file_mah = mah_path + str(ID) + mah_format

        try:
            # Read the rest of the file
            this_mah = pd.read_csv(file_mah, sep='\t', skiprows=1, header=None)

            def split_id(x):
                x = x.split(' ')
                new_x = list(filter(None, x)) 
                return new_x[0]

            def split_host(x):
                x = x.split(' ')
                new_x = list(filter(None, x)) 
                return new_x[1]

            new_mah = pd.DataFrame() 

            # For sure there is a smarter way to do this but I could not figure it out so far
            new_mah['ID'] = this_mah.loc[:, 0].apply(lambda x: split_id(x))   
            new_mah['HostHalo'] = this_mah.loc[:, 0].apply(lambda x: split_host(x))   
            
            # Also this could be done more efficiently but let's not waste too much time on this part
            for i in range(2, n_cols):
                new_mah[cols[i]] = this_mah.loc[:, i-1]

            halo_mah = hu.HaloHistory(time)
            halo_mah.load_full_mah(new_mah)

            all_mah.append(halo_mah)

        except:
            'The output file for this halo was not produced (too few particles to be traced, most likely). Do nothing'

    return all_mah



"""
    This is sort of useless but allows to keep some simmetry in the program reading routines
"""
def read_csv_tree(file_name):
    tree = pd.read_csv(file_name)

    return tree



"""
    Read a redshift file (containing the z info of each snapshot) and a tabulated file containing a(t) each 5 Myrs
"""
def read_time(data_path):
    f_a2myr = data_path + 'output_list_5Myr.txt'
    f_z = data_path + 'redshifts_sims.txt'

    df_z = pd.read_csv(f_z, header = None)
    df_a2myr = pd.read_csv(f_a2myr, header = None)

    def z2a(z):
        return (1.0 / (1.0 + z))

    df_a = df_z.apply(lambda z: z2a(z))

    n_t = df_a2myr.shape[0]
    t2a = np.zeros((2, n_t))

    T0 = 13.75   # Time after the big bang
    step = 0.005 # 5 Megayears
    t2a[1, :] = df_a2myr.values.T
    d0 = abs(T0 - n_t * step)
    
    for i in range(1, n_t+1):
        t2a[0, i-1] = d0 + i * step

    n_z = df_z.shape[0]
    time = np.zeros((3, n_z))

    for i in range(0, n_z):
        time[0, i] = df_z.values[n_z -i -1]
        time[1, i] = df_a.values[n_z -i -1]

    time[2, :] = np.interp(time[1, :], t2a[1, :], t2a[0, :])

    return time


'''
    Return the full particle content of a snapshot as a dataframe
'''
def read_snap(file_name=None, velocity=False, part_types=[1], n_files=1):

    # Make this a list if it is not already a list
    if isinstance(part_types, list) == False:
        part_types = [part_types]

    full_data = None

    # Loop over snapshots and particle types
    for f in range(0, n_files):
        if n_files > 1:
            this_file = file_name + '.' + str(f)
        else:
            this_file = file_name

        for part_type in part_types:

            try:
                # Read positions only
                particles = readsnap(this_file, 'pos', part_type)
                print('Reading positions from file: ', this_file, ' ptype = ', part_type)
                
                # Add particle type info
                ptype = np.zeros((len(particles), 1))
                ptype.fill(part_type)

                # Do we want to read in velocity data as well? Or positions only?
                if velocity == True:
                    cols = ['X', 'Y', 'Z', 'Vx', 'Vy', 'Vz', 'Type']

                    velocities = readsnap(this_file, 'vel', part_type)
                    print('Reading velocities from file: ', this_file, ' ptype = ', part_type)

                    # Convert to dataframe
                    if full_data == None:
                        part_vel = np.concatenate((particles, velocities), axis=1)
                        full_data = np.concatenate((part_vel, ptype), axis=1)
                    else:
                        new_data = np.concatenate((particles, velocities), axis=1)
                        new_data = np.concatenate((new_data, ptype), axis=1)
                        full_data = np.concatenate((full_data, new_data), axis=0)

                # Reading only positions
                else:
                    cols = ['X', 'Y', 'Z', 'Type']

                    # Concatenate
                    if full_data == None:
                        full_data = np.concatenate((particles, ptype), axis=1)
                    else:
                        new_data = np.concatenate((particles, ptype), axis=1)
                        full_data = np.concatenate((full_data, new_data), axis=0)

            # Skip file
            except:
                print('There are no particles of type: ', part_type, ' in file: ', this_file)

    try:
        n_part = len(full_data)
        print('Found ', n_part, ' particles in total, for type(s): ', part_types)
        part_df = pd.DataFrame(data=full_data, columns=cols)

    except:
        print('Warning. No particles found in file: ', this_file)
        part_df = pd.DataFrame()

    # Return the selected particles' properties in a dataframe
    return part_df
