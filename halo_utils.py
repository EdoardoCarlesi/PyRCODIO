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

    # This one holds information about z/a/t
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

'''
    A Local Group model is a way of selecting LG like pairs based on a series of criteria.
    They are all listed here with some generic numerical values
'''
class LocalGroupModel:
    d_iso = 2000.
    r_max = 1250.
    r_min = 250.
    m_max = 5.e+12
    m_min = 5.e+11
    mtot_max = 8.e+12
    mratio_max = 4.
    vrad_max = -1.0
    center = [50000.0] * 3
    model_code = '00'

    def __init__(self, d_iso, r_max, r_min, m_max, m_min, mratio_max, vrad_max):
        self.d_iso = d_iso
        self.r_max = r_max
        self.r_min = r_min
        self.m_max = m_max
        self.m_min = m_min
        self.mratio_max = mratio_max
        self.vrad_max = vrad_max

    def info(self, dump=False):
        lg_mod_info = "D_iso: %.3f, Rh_max: %.3f, Rh_min: %.3f, M_min: %.3f\n" % \
                        (self.d_iso, self.r_max, self.r_min, self.m_min)
        if dump == True:
            print(lg_mod_info)

        return lg_mod_info


'''
    A local group object is a pair of halos, we use this class to be able to quickly compute some of their properties.
'''
class LocalGroup:
    # Initialize to some numerical value
    vrad = -100.
    r = 770.
    d_cbox = 0.0
    d_virgo = 0.0
    rating = 0.0
    hubble = 0.67

    # Init to None
    LG1 = None
    LG2 = None
    com = None
    code = None
    ahf_file = None

    # The constructor needs two halos (lg1, lg2)
    def __init__(self, lg1, lg2):

        if lg1.m > lg2.m:
            self.LG1 = lg1
            self.LG2 = lg2
        else:
            self.LG1 = lg2
            self.LG2 = lg1

        if lg1.x[0] != 0 and lg2.x != 0:
            self.r = self.r_halos()
            self.vrad = self.v_radial()
            self.com = self.get_com()

    def rating(self):
        self.rating = rate_lg_pair(self.LG1, self.LG2)
        return self.rating

    def geo_com(self):
        return 0.5 * (self.LG1.pos() + self.LG2.pos())

    def com(self):
        self.com = center_of_mass([self.LG1.m, self.LG2.m], [self.LG1.x, self.LG2.x])
        return self.com

    def lg_member(self, n_member):
        if n_member == 0:
            return self.LG1
        if n_member == 1:
            return self.LG2

    def r_halos(self):
        self.r = self.LG1.distance(self.LG2.pos())
        return self.r

    def v_radial(self):
        self.vrad = vel_radial(self.LG1.pos(), self.LG2.pos(), self.LG1.vel(), self.LG2.vel())
        self.vrad += self.hubble * self.r * 0.1
        return self.vrad

    def m_ratio(self):
        if self.LG1.m > self.LG2.m:
            return (self.LG1.m / self.LG2.m)
        else:
            return (self.LG1.m / self.LG2.m)

    def header(self):
        n_head = 0
        header = '#'
        header += ' SimuCode('+ str(n_head) +')' ; n_head = +1
        header += ' ID_M31('+ str(n_head) +')' ; n_head = +1
        header += ' ID_MW('+ str(n_head) +')' ; n_head = +1
        header += ' R_MWM31('+ str(n_head) +')' ; n_head = +1
        header += ' Vrad('+ str(n_head) +')' ; n_head = +1
        header += ' Nsub_M31('+ str(n_head) +')' ; n_head = +1
        header += ' Nsub_MW('+ str(n_head) +')' ; n_head = +1
        header += ' X_com('+ str(n_head) +')' ; n_head = +1
        header += ' Y_com('+ str(n_head) +')' ; n_head = +1
        header += ' Z_com('+ str(n_head) +')' ; n_head = +1

        return header

    def info(self):
        h0 = self.hubble
        kpc = 1000.
        file_lg_line = '%s  %ld   %ld   %7.2e   %7.2e   %7.2f   %7.2f   %5d   %5d  %5.2f  %5.2f  %5.2f  %d  %d' % \
                (self.code, self.LG1.ID, self.LG2.ID, self.LG1.m/h0, self.LG2.m/h0, self.r/h0, self.vrad, \
                        self.LG1.nsub, self.LG2.nsub, self.com[0]/kpc, self.com[1]/kpc, self.com[2]/kpc, \
                        self.LG1.npart, self.LG2.npart)

        if self.LG1.m_star > 0.0:
            extra_mstar = '   %7.2e   %7.2e' % (self.LG1.m_star/h0, self.LG2.m_star/h0)
            extra_mgas = '   %7.2e   %7.2e' % (self.LG1.m_gas/h0, self.LG2.m_gas/h0)
            file_lg_line = file_lg_line + extra_mstar + extra_mgas

        return file_lg_line

    def m_tot(self):
        return (self.LG1.m + self.LG2.m)


'''
    This functions allows to select among "best" and "worst" LG candidates using fiducial benchmark values.
    These numbers and parameters are - to a certain extent - arbitrary.
'''
def rate_lg_pair(lg1, lg2):
    # Benchmark (obs.) quantities
    rhalo0 = 500.   # kpc/h
    vrad0 = -100.
    mass0 = 3.0e+12
    ratio0 = 1.1
    hubble0 = 67.0

    npart = 512

    com = center_of_mass([lg1.m, lg2.m], [lg1.x, lg2.x])
    m12 = lg1.m + lg2.m

    if lg1.m > lg2.m:
        rm12 = lg1.m / lg2.m
    else:
        rm12 = lg2.m / lg1.m

    rhalos = lg1.distance(lg2.x)
    vrad = vel_radial(lg1.x, lg2.x, lg1.v, lg2.v)
    vrad += hubble0 * rhalos/1000.

    # Relative weights to estimate the relevance of each property relative to the "real" LG
    fac_rh = 1.0
    fac_c = 0.25
    fac_v = 0.25
    fac_m = 1.5
    fac_ra = 1.5

    # How to compute the LG-likeliness factors
    diff_rh = abs(rhalos - rhalo0) / rhalo0
    diff_m = np.log10(m12 / mass0)
    diff_v = abs(vrad0 - vrad) / abs(vrad0)
    diff_ra = abs(rm12 - ratio0) / abs(ratio0)

    lg_rate = diff_rh * fac_rh + diff_m * fac_m + diff_ra * fac_ra + fac_v * diff_v

    # Get a penalty for positive vrad
    if vrad > 5.0:
        lg_rate += 10.

    return lg_rate


'''
    Given a point in space return the halo IDS around it, very simple.
'''
def halo_ids_around_center(halos, center, radius):
    ids = []

    for h in halos:
        d = h.distance(center)
        if d < radius:
            ids.append(h.halo_df['ID'])

    return ids

