'''
    Python Routines for COsmology and Data I/O (PyRCODIO)
    Pandas Upgrade
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    halo_utils.py: various utilities and simple computational routines used to compute halo trajectories and so on
'''

import tools as t
import pandas as pd
import numpy as np

'''
    This is the main class that stores all the halo properties
'''
class Halo:
    halo_df = pd.DataFrame()
    host_df = pd.DataFrame()
    subhalos = []

    # We must assign the main halo, if this is a subhalo then the host can be also assigned
    def __init__(self, halo, host=None):
        self.halo_df = halo

        if host != None:
            self.host_df = host

    # Print the head of this DF - it hols only one element (row) so it will print everything
    def info(self):
        print(self.halo_df.head())

    # Return the position vector
    def pos(self):
        return self.halo_df[['Xc(6)', 'Yc(7)', 'Zc(8)']].T.values

    # Return the velocity vector
    def vel(self):
        return self.halo_df[['VXc(9)', 'VYc(10)', 'VZc(11)']].T.values

    # Simply return r
    def r(self):
        return self.halo_df[['Rvir(12)']].values[0, 0]

    # Simply return the mass
    def m(self):
        return self.halo_df[['Mvir(4)']].values[0, 0]

    # Compute the distance from a given point
    def distance(self, c):
        return t.distance(self.pos(), c)
    
    # Shortcut to the dataframe structure
    def find(self, key):
        return self.halo_df[key]

    # This function looks for subhalos of the host given a halos dataframe halo catalog
    def assign_subhalos(self, halos, r_sub = None):
        if r_sub == None:
            r_sub = self.r()

        self.subhalos = t.find_halos(self.pos(), r_sub, halos)


'''
    This class is a wrapper that contains (for a halo or a subhalo) the MAH, plus positions, velocities and more
'''
class HaloHistory:
    n_steps = 0

    host = None
    halos = []
    subhalos = []
    redshift = []
    time = []
    n_mergers = []

    # Constructor, we need only the number of steps, the rest is allocated by default
    def __init__(self, n_steps, host = None):
        self.n_steps = n_steps

        # Initialize some empty lists
        self.halos = []
        self.subhalos = []
        self.redshift = []

        # Time in Gigayears
        self.time = []

        # Number of mergers per time step
        self.n_mergers = []

        # If the halo is a subhalo then we can add an host
        if host != None:
            self.host = host

    # Simply return the mass evolution as a function of time
    def m_t(self):
        m_t = np.zeros((self.n_steps))

        for i in range(0, self.n_steps):
            m_t[i] = halos[i].m()

        return m_t
    
    # Return the max mass value
    def m_max(self):
        m_max = np.max(self.m_t())
        return m_max

    # Halo velocity as a function of time
    def v_t(self):
        vel_t = np.zeros((self.n_steps, 3))

        for i in range(0, self.n_steps):
            vel_t[i] = self.halos[i].vel()

        return vel_t

    # Halo position as a function of time
    def x_t(self):
        pos_t = np.zeros((self.n_steps, 3))

        for i in range(0, self.n_steps):
            pos_t[i] = self.halo[i].pos()

        return pos_t

    # Last major merger time, assuming a mass ratio of 0.1
    def last_major_merger(self):
        lmm = 0

        return lmm

    # Formation time, with a mass at half time of 0.5
    def formation_time(self):
        ft = 0

        return ft

    # Compute the path of the halo in the host halo reference system
    def trajectory_around_host(self):

        if self.host != None:
            tah = self.x_t() - self.host.x_t()
        else:
            print('Halo ID ', self.halos[0]['ID'].values[0], ' is not a subhalo.')
            return -1

        return tah





