'''
    Python Routines for COsmology and Data I/O (PyRCODIO)
    Pandas Upgrade
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    particles.py: functions and more to manipulate ensambles of particles
'''

import tools as t
import pandas as pd
import numpy as np
import scipy as sp


'''
    Given a point in space return the particles around it.
    Input: parts is a data structure
'''
def find_particles(parts, center, radius):

    new_key = 'Distance'

    def dist(x, c):
        return t.distance(x, c)

    parts[new_key] = parts[['Xc(6)', 'Yc(7)', 'Zc(8)']].T.apply(dist, c=center).T

    return catalog[catalog[new_key] < radius]


'''
    Given a dataframe of particles, compute the overdensity around a given center, at a given radius
'''
def overdensity(part_df=None, R=None, center=None, rho0=None):

    if center == None:
        center = 'center of mass'

    volume = (4.0 / 3.0) * np.pi * R * R * R
    parts = find_particles(part_df, center, R)
    nparts = len(parts)

    delta = (np.power(nparts, 3.0) / volume) / rho0

    return delta


'''
    Data must be a dataframe and axes a list of columns
'''
def inertia_tensor(data=None, axes=None):

    X = axes[0]
    Y = axes[0]
    Z = axes[0]

    I = np.zeros((3, 3))

    for i_ax in range(0, 3):
        j_ax = (i_ax + 1) % 3
        k_ax = (i_ax + 2) % 3

        I[i_ax, i_ax] = np.sum(data[j_ax] ** 2.0 + data[k_ax] **2.0)
 

    return I



def triaxiality():

    return t


def simu_rho0(box=None, npart=None):

    rho0 = np.power(npart/box, 3)

    return rho0
