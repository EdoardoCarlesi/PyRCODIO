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

    # We must assign the main halo, if this is a subhalo then the host can be also assigned. Halo is a dataframe!
    def __init__(self, halo, host=None):
        self.halo_df = halo

        if host != None:
            self.host_df = host

    # Print the head of this DF - it hols only one element (row) so it will print everything
    def info(self):
        print(self.halo_df.head())

    # Return the ID of the halo
    def ID(self):
        return self.halo_df[['ID']].values[0]

    # Number of subhalos
    def nsub(self):
        return self.halo_df[['numSubStruct(3)']].values[0]

    # Number of particles
    def npart(self):
        return self.halo_df[['npart(5)']].values[0]

    # Vmax
    def vmax(self):
        return self.halo_df[['Vmax(17)']].values[0]

    # r2 (scale radius)
    def r2(self):
        return self.halo_df[['r2(14)']].values[0]

    # lambda (scale radius)
    def lambdap(self):
        return self.halo_df[['lambda(20)']].values[0]

    # b (triaxiality b)
    def b(self):
        return self.halo_df[['b(25)']].values[0]

    # c (triaxiality c)
    def c(self):
        return self.halo_df[['c(26)']].values[0]

    # c_NFW (concentration)
    def c_NFW(self):
        return self.halo_df[['cNFW(43)']].values[0]

    # Return the position vector
    def pos(self):
        return self.halo_df[['Xc(6)', 'Yc(7)', 'Zc(8)']].T.values

    # Return the velocity vector
    def vel(self):
        return self.halo_df[['VXc(9)', 'VYc(10)', 'VZc(11)']].T.values

    # Simply return r
    def r(self):
        return float(self.halo_df[['Rvir(12)']].values[0])

    # Simply return the mass
    def m(self):
        return self.halo_df[['Mvir(4)']].values[0]

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

        self.subhalos = find_halos(halos, self.pos(), r_sub)


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
    vtan = 0.
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

        try:
            if lg1.m() > lg2.m():
                self.LG1 = lg1
                self.LG2 = lg2
            else:
                self.LG1 = lg2
                self.LG2 = lg1

            if lg1.pos()[0] != 0 and lg2.pos()[0] != 0:
                self.r = self.r_halos()
                self.vrad = self.vel_radial()
                self.vtan = self.vel_tangential()
                self.com = self.mass_com()
        except:
            'We might need some void instance of LocalGroup'

    def rating(self):
        self.rating = rate_lg_pair(self.LG1, self.LG2)
        return self.rating

    def geo_com(self):
        return 0.5 * (self.LG1.pos() + self.LG2.pos())

    def mass_com(self):
        self.com = t.center_of_mass([self.LG1.m(), self.LG2.m()], [self.LG1.pos(), self.LG2.pos()])
        return self.com

    def lg_member(self, n_member):
        if n_member == 0:
            return self.LG1
        if n_member == 1:
            return self.LG2

    def r_halos(self):
        self.r = self.LG1.distance(self.LG2.pos())
        return self.r

    def vel_radial(self, hubble=True):
        self.vrad = t.vel_radial(self.LG1.pos(), self.LG2.pos(), self.LG1.vel(), self.LG2.vel())
            
        # Add the effects of the Hubble expansion by default
        if hubble == True:
            self.vrad += self.hubble * self.r * 0.1
        return self.vrad

    def vel_tangential(self):
        vrad0 = self.vel_radial(hubble=False)
        vtot = t.vec_module(t.vec_subt(self.LG1.vel(), self.LG2.vel()))
        self.vtan = np.sqrt(vtot**2.0 - vrad0**2.0)

        return self.vtan

    def m_ratio(self):
        if self.LG1.m > self.LG2.m:
            return (self.LG1.m / self.LG2.m)
        else:
            return (self.LG1.m / self.LG2.m)

    def header(self, csv=True):
        
        if csv == True:
            header =  'M_M31,' 
            header += 'M_MW,'
            header += 'R,'
            header += 'Vrad,'
            header += 'Vtan,'
            header += 'Nsub_M31,'
            header += 'Nsub_MW,'
            header += 'Npart_M31,'
            header += 'Npart_MW,'
            header += 'Vmax_MW,'
            header += 'Vmax_M31,'
            header += 'lambda_MW,'
            header += 'lambda_M31,'
            header += 'cNFW_MW,'
            header += 'Xc_LG,'
            header += 'Yc_LG,'
            header += 'Zc_LG,'
            header += 'simu_code,'
            header += 'sub_code'

        else:
            # TODO : implement some readable header
            #header += ' X_com, '
            #header += ' Y_com, '
            #header += ' Z_com, '
            n_head = 0
            

        return header

    def info(self, dump=True, csv=True):
        h0 = self.hubble
        kpc = 1000.

        if csv == True:
            file_lg_line = '%7.2e, %7.2e, %7.2f, %7.2f, %7.2f, %d, %d, %d, %d, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f, %5.2f,%s' % \
                       (self.LG1.m(), self.LG2.m(), self.r_halos(), \
                        self.vel_radial(), self.vel_tangential(), self.LG1.nsub(), self.LG2.nsub(), \
                        self.LG1.npart(), self.LG2.npart(), self.LG1.vmax(), self.LG2.vmax(),\
                        self.LG1.lambdap(), self.LG2.lambdap(), self.LG1.c_NFW(), self.LG2.c_NFW(),\
                        self.geo_com()[0], self.geo_com()[1], self.geo_com()[2], self.code)
        else:
            # TODO : implement some readable header
            n_head = 0

        # TODO implement this
        ''' 
        if self.LG1.m_star > 0.0:
            extra_mstar = '   %7.2e   %7.2e' % (self.LG1.m_star/h0, self.LG2.m_star/h0)
            extra_mgas = '   %7.2e   %7.2e' % (self.LG1.m_gas/h0, self.LG2.m_gas/h0)
            file_lg_line = file_lg_line + extra_mstar + extra_mgas
        '''

        if dump == True:
            print(file_lg_line)

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
    Given a point in space return the halo IDS around it.
    Input: halos is a list of halo objects
'''
def halo_ids_around_center(halos, center, radius):
    ids = []

    for h in halos:
        d = h.distance(center)
        if d < radius:
            ids.append(h.halo_df['ID'])

    return ids

'''
    Given a center, find all the halo within a given radius
    Input: catalog is a DataFrame
'''
def find_halos(catalog, center, radius):

    new_key = 'Distance'
    cols = ['Xc(6)', 'Yc(7)', 'Zc(8)']

    def dist(x, c):
        return t.distance(x, c)

    # Check units
    mpc2kpc = 1.0e+3
    n_col = int(0.5 * len(catalog))

    x_sum = np.sum(catalog[cols].iloc[n_col])
    c_sum = np.sum(center)

    # If units are Mpc convert to kpc and return the "right" units
    if x_sum < 1.0e+4:
        catalog[cols] = catalog[cols].apply(lambda x: x * mpc2kpc)

        # Also the viral radius needs to be ported to the correct units
        catalog['Rvir(12)'] = catalog['Rvir(12)'].apply(lambda x: x * mpc2kpc)
        n_sum = np.sum(catalog[cols].iloc[n_col])

    catalog[new_key] = catalog[cols].T.apply(dist, c=center).T
    #print(radius)
    #print(catalog[new_key])

    return catalog[catalog[new_key] < radius]    

'''
    Given a box and a model, return a list of possible local groups
    Taken from the old PyRCODIO routines
'''
def find_lg(halos, lgmod, center, radius):
    center = lgmod.center
    iso_radius = lgmod.d_iso
    m_min = lgmod.m_min
    m_max = lgmod.m_max
    r_min = lgmod.r_min
    r_max = lgmod.r_max
    model_code = lgmod.model_code	
    vrad_max = lgmod.vrad_max

    # First select only halos within a given radius from the center
    halos_center = find_halos(halos, center, radius)

    # Refine the selection choosing only halos above the minimum mass
    halos_mass = halos_center[halos_center['Mvir(4)'] > m_min]
    
    n_candidates = halos_mass['Mvir(4)'].count()
    print('Looking for candidates among ', n_candidates, ' halos.')

    # Sanity check
    #print(halos_mass.info())
    #print(halos_mass.head())
    #pos = halos_mass.iloc[1, 3:6].values
    #print('Pos: ', pos)

    # These are initialized empty
    lgs = []		# Haloes paired as LG structures

    # Now loop on halo center candidates
    for h in range(0, n_candidates):
        # We want to create a single halo object that needs a DataFrame for the constructor
        halo_lg0 = Halo(halos_mass[halos_mass['ID'] == halos_mass.iloc[h, -3]])
        count_lg = 0    
        count_wrong = 0
		
        if halo_lg0.m() > m_max:
            count_wrong = 1
    
        # We need to run the loop on ALL halos, in case there is a third halo lying close
        for i in range(h+1, n_candidates):
            halo_lg1 = Halo(halos_mass[halos_mass['ID'] == halos_mass.iloc[i, -3]])
            dis_this = halo_lg1.distance(halo_lg0.pos())

            # This is a subhalo! Keep it				
            if dis_this < halo_lg0.r() or dis_this < r_min:
                # Just do nothing
                count_wrong += 0
            
            # This is a possible second candidate
            elif dis_this > r_min and dis_this < r_max:
                halo_lg2 = halo_lg1	
                vrad = t.vel_radial(halo_lg2.pos(), halo_lg0.pos(), halo_lg2.vel(), halo_lg0.vel()) 
                count_lg += 1
				
                # There are too many close-by haloes
                if count_lg > 1:
                    count_wrong += 1

                # Check if they satisfy the v radial condition
                if vrad + (0.67 * dis_this * 0.1) > vrad_max:
                    count_wrong += 1

		# A new first & second LG halos have been found:
                if count_wrong == 0 and count_lg == 1:
                    already_there = 0

                    '''
                    for j in range(0, n_all_lgs):
                        lg1 = halos_lg[2 * j]
                        lg2 = halos_lg[2 * j + 1]
			
                        if halo_lg0.ID == lg1.ID or halo_lg0.ID == lg2.ID:
                            already_there += 1
                    '''

                    # We have a candidate! 
                    if already_there == 0:
                        n_iso_radius = 0

                    # Check also for third haloes within the isolation radius before adding the pair
                    this_lg = LocalGroup(halo_lg0, halo_lg2)
                    com = this_lg.geo_com()
                    halos_iso = find_halos(halos, com, iso_radius)
                    nh_iso = len(halos_iso)

                    n_iso_radius = halos_iso[halos_iso['Mvir(4)'] > m_min]['Mvir(4)'].count()
	
                    # Make sure there are only two halos within this radius
                    if n_iso_radius == 2: 
                        lgs.append(this_lg)

    return lgs
