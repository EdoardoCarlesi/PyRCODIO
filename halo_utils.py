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
    mainhalo = None
    time = None
    subhalos = []
    n_mergers = []

    # Constructor, we need only the number of steps, the rest is allocated by default
    def __init__(self, time, host = None):
        self.n_steps = 0

        # Initialize some empty lists
        self.mainhalo = None
        self.subhalos = []

        # a, z and Time in Gigayears 
        self.time = time

        # Number of mergers per time step
        self.n_mergers = []

        # If the halo is a subhalo then we can add an host
        if host != None:
            self.host = host

    # When reading the complete MAH from a dataframe
    def load_full_mah(self, mah_df):
        self.mainhalo = mah_df
        self.n_steps = mah_df.shape[0]

    def set_time(self, time):
        self.time = time

    # Simply return the mass evolution as a function of time
    def m_t(self):
        mvir_name = 'Mvir(4)'

        return self.mainhalo[mvir_name]
 
    # Return the mass normalized to m(z=0)
    def m_t_norm(self):
        m_t = self.m_t() / self.m_t()[0]

        return m_t
    
    # Return the max mass value
    def m_max(self):
        m_max = np.max(self.m_t())

        return m_max
 
    # Return the max mass value
    def t_m_max(self):
        m_max = np.max(self.m_t())
        m_ind = np.where(self.m_t() == m_max)
        t_m_max = self.time[2, m_ind]

        return t_m_max

    # Halo velocity as a function of time
    def v_t(self):
        vel_t = np.zeros((self.n_steps, 3))

        vel_t[:, 0] = self.mainhalo['VXc(9)']
        vel_t[:, 1] = self.mainhalo['VYc(10)']
        vel_t[:, 2] = self.mainhalo['VZc(11)']

        return vel_t

    # Halo position as a function of time
    def x_t(self):
        pos_t = np.zeros((self.n_steps, 3))

        pos_t[:, 0] = self.mainhalo['Xc(6)']
        pos_t[:, 1] = self.mainhalo['Yc(7)']
        pos_t[:, 2] = self.mainhalo['Zc(8)']

        return pos_t

    # Last major merger time, assuming a mass ratio of 0.1
    def last_major_merger(self):
        m_t = self.m_t_norm().values
        thr = 0.1

        for i in range(1, self.n_steps):
            m0 = m_t[i-1]
            m1 = m_t[i]
            mdiff = abs((m0 - m1) / m0)

            if mdiff > thr:
                lmm = 0.5 * (self.time[2, i-1] + self.time[2, i])
                break

        return lmm

    # Formation time, with a mass at half time of 0.5
    def formation_time(self):
        m_t = self.m_t_norm().values
        f_t = m_t[np.where(m_t > 0.5)]
        n_ft = len(f_t)

        # Take half step in between
        f_t = 0.5 * (self.time[2, n_ft-1] + self.time[2, n_ft])

        return f_t

    # Compute the path of the halo in the host halo reference system
    def trajectory_around_host(self):
        if self.host != None:
            tah = self.x_t() - self.host.x_t()
        else:
            print('Halo ID ', self.halos[0]['ID'].values[0], ' is not a subhalo.')
            return -1

        return tah





