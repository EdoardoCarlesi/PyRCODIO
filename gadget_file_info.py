'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_particles_extract.py: extract particle data within a high-res region in the ICs or in a snapshot
'''

import read_files as rf
import sys

file_in = sys.argv[1]
part_type = sys.argv[2]
 
     part_df = rf.read_snap(file_name=this_file, velocity=velocity, part_types=part_type, n_files=n_files)

