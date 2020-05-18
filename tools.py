'''
    Python Routines for COsmology and Data I/O (PyRCODIO)
    Pandas Upgrade
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    tools.py: various utilities and simple computational routines used throughout the code
'''

import pandas as pd
import numpy as np

'''
    Compute the Euclidean distance between two points in space
'''
def distance(x, center):
    dist = 0;

    for i in range(0, len(x)):
        dist += (x[i] - center[i])**2.0

    dist = np.sqrt(dist)

    return dist


'''
    Given a (halo) center at pos_host, find all the (sub) halos within a given radius
'''
def find_halos(pos_host, radius, catalog):

    new_key = 'Distance'

    def dist(x, c):
        return distance(x, c)

    catalog[new_key] = catalog[['Xc(6)', 'Yc(7)', 'Zc(8)']].T.apply(dist, c=pos_host).T
    
    return catalog[catalog[new_key] < radius]    

