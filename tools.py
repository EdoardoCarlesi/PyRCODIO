'''
    Python Routines for COsmology and Data I/O (PyRCODIO) v0.2
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    tools.py: various utilities and simple computational routines used throughout the code
'''

import read_files as rf
import pandas as pd
import numpy as np
import random

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def select_box(data=None, radius=None, col=None, x_col=None, center=None):
    """ Given a center and a radius, select halos in a dataframe """

    condition = (data[x_col[0]] > (center[0] - radius)) & (data[x_col[0]] < (center[0] + radius)) &\
            (data[x_col[1]] > (center[1] - radius)) & (data[x_col[1]] < (center[1] + radius)) &\
            (data[x_col[2]] > (center[2] - radius)) & (data[x_col[2]] < (center[2] + radius))

    selected = data[condition]
    return selected


def select_sphere(data=None, radius=None, col=None, x_col=None, center=None):
    """ Given a center and a radius, select halos in a dataframe """

    all_x = data[x_col].T.values
    all_x_d = np.sum((all_x - center) ** 2.0, axis=0)
    data[col] = all_x_d
    data[col] = data[col].apply(lambda x: np.sqrt(x))
    selected = data[data[col] < radius]

    return selected


def apply_distance(data=None, x_col=None, center=None, col=None):

    all_x = data[x_col].T.values
    all_x_d = np.sum((all_x - center) ** 2.0, axis=0)
    data[col] = all_x_d
    data[col] = data[col].apply(lambda x: np.sqrt(x))

    return data[col]

def spatial_pca(data=None, cols=None):
    """ Do a PCA analysis of the coordinates to find out asymmetries in the halo distribution """

    scaler = StandardScaler()
    x = data[cols].values
    x = scaler.fit_transform(x)
    
    n_x = len(x)

    if n_x > 3:

        pca = PCA(n_components=3)
        principal = pca.fit_transform(x)
        axs = pca.explained_variance_ratio_

        axx = axs #/ axs[0]
    
    else:
        axx = [0.32, 0.33, 0.34]

    return axx


def std_pca(x=None) : 
    """ Do a PCA analysis of the coordinates to find out asymmetries in the halo distribution """

    #scaler = StandardScaler()
    #x = scaler.fit_transform(x)
    n_comp = len(x)
    pca = PCA(n_components=n_comp)
    principal = pca.fit_transform(x)
    axs = pca.explained_variance_ratio_

    axx = axs / axs[0]

    return axx


def triaxiality(a, b, c):
    """ Franx et al. 1991 defintion """

    return (a ** 2.0 - b **2.0) / (a ** 2.0 - c ** 2.0)


def inertia_tensor(x=None, w=[]): 
    """ Compute the moment of inertia of a mass distribution of halos and get the eigenvalues """

    I = np.zeros((3, 3))

    # This is the unweighted Inertia Tensor, just set w to one
    if len(w) < 1:
        w = np.ones(len(x[0, :]))

    I[0][0] = np.sum(w * (x[1, :] **2 +  x[2, :] **2))
    I[1][1] = np.sum(w * (x[0, :] **2 +  x[2, :] **2))
    I[2][2] = np.sum(w * (x[1, :] **2 +  x[0, :] **2))

    I[1][0] = -np.sum(w * (x[1, :] * x[0, :]))
    I[1][2] = -np.sum(w * (x[1, :] * x[2, :]))
    I[0][2] = -np.sum(w * (x[0, :] * x[2, :]))

    I[0][1] = I[1][0]
    I[2][1] = I[1][2]
    I[2][0] = I[0][2]

    evs = np.linalg.eigvals(I)

    evs /= max(evs)

    return evs


def bin_df(data=None, x_bins=None, col='Mvir(4)', binmode='log'):
    """ Count the number of entries of a data structure given an array of bins """

    nbins = len(x_bins)
    binned = np.zeros((nbins-1))
    col_bin = 'binned'

    if binmode == 'log':
        data['Mlog'] = np.log10(data[col])
        col = 'Mlog'

    for i in range(0, nbins-1):
        binned[i] = len(data[(data[col] > x_bins[i])]) 

    return binned


def gen_bins(nbins=None, binmax=None, binmin=None, binmode='log'):
    """ Simple tool to generate bin intervals """

    if binmode == 'log':
        bmax = np.log10(binmax)
        bmin = np.log10(binmin)

    else:
        bmax = binmax
        bmin = binmin

    bins = np.zeros((nbins))
    step = (bmax - bmin) / float(nbins)
    bins[0] = bmin

    for i in range(1, nbins):
        bins[i] = bins[i-1] + step

    return bins


def density(data=None, size=None, col='Mvir(4)', shape='cube'):
    """ Given a halo df and a volume determine the matter density """

    mtot = np.sum(data[col])

    if shape == 'cube':
        vol = size ** 3.0

    elif shape == 'sphere':
        vol = 4.0 / 3.0 * (np.pi) * (size ** 3.0)

    dens = mtot / vol

    return dens


def distance(x, c):
    """ Compute the Euclidean distance between two points in space """

    return np.sqrt(np.sum((x-c)**2.0))


def shift(center, r):
    """ Given a set of coordinates, randomly shift them by a maximum of 'r' amount """

    new_c = []

    for c in center:
        eps = random.randrange(-r, r)
        c = c + eps
        new_c.append(c)

    return np.array(new_c)


def module(vec):
    """ Very basic operation, there is for sure some quicker way of implementing this but whatever """

    return np.sqrt(np.sum(vec **2.0))


def find_nearest_node_index(x=None, grid=None, box=None):
    """ Given a point x in space, find the nearest grid point once a grid has been placed on the box """

    cell = box / grid
    ix = np.floor(x[0] / cell)
    iy = np.floor(x[1] / cell)
    iz = np.floor(x[2] / cell)

    index = int(ix + grid * iy + grid * grid * iz)

    return index


def angle(v1, v2):
    """ Returns the angle between two vectors """

    mv1 = [0.0] * 3
    mv2 = [0.0] * 3
    mod1 = module(v1)
    mod2 = module(v2)

    for i in range(0, 3):
        mv1[i] = v1[i] / mod1
        mv2[i] = v2[i] / mod2

    v12 = dot_prod(mv1, mv2)

    return v12


def center_of_mass(m, x):
    """ Yet another simple function """

    n = len(m)
    com = [0.0] * 3

    for j in range(0, 3):
        mtot = 0.0

        for i in range(0, n):
            mtot += m[i]
            com[j] += x[i][j] * m[i]

    com[j] /= mtot

    return com


def vel_radial(x1, x2, v1, v2):
    """ Radial velocity component of a two-object system """

    x12 = x2 - x1
    n12 = np.sqrt(np.sum(x12 * x12))
    r12 = np.sum((v2 - v1) * x12)
    nr12 = r12 / n12

    return nr12


def mass_function(masses):
    """ Given a vector of masses, it returns a tuple of mass vs cumulative number of objects """

    n_m = len(masses)
    y_n = [0 for im in range(0, n_m)]
    x_m = sorted(masses)

    for im in range(0, n_m):
        y_n[im] = n_m - im

    return x_m, y_n


def particles_com(part_df, cols=['X', 'Y', 'Z'], mass_types=1):
    """ Compute the center of mass of a particle distribution """

    com = [0.0] * 3

    if mass_types > 1:
        print('Error. Only available for ONE particle mass for all the distribution')
        return 0

    for i, col in enumerate(cols):
        com[i] = part_df[col].mean()

        # Convert to kpc in case
        if com[i] < 1.e+4:
            com[i] = com[i] * 1.e+3

    return np.array(com)


def find_slab(part_df=None, center=None, side=None, thick=None, rand_seed=69, reduction_factor=1.0, z_axis=2, units='kpc', cols=['X', 'Y', 'Z']):
    """
    Find the particles belonging to a slab around a given point in space.
    Slab size, thickness and so on need to be specified.
    """

    # Set some parameters
    kpcThresh = 1.e+4
    kpc2Mpc = 1.e-3

    minima = np.zeros((3))
    maxima = np.zeros((3))

    # Select the two axes for the 2D projection
    ax0 = (z_axis + 1) % 3
    ax1 = (z_axis + 2) % 3
    ax2 = z_axis

    # Column names
    col0 = cols[ax0]
    col1 = cols[ax1]
    col2 = cols[ax2]

    n_part = len(part_df)

    # Sanity check on the units
    half_n = int(n_part * 0.5)
    sum_coord = part_df[col0].iloc[half_n] + part_df[col1].iloc[half_n] + part_df[col2].iloc[half_n]

    # Make sure the units are consistent
    if sum_coord < kpcThresh:
        side = side * kpc2Mpc
        center = center * ([kpc2Mpc] *3)
        thick = thick * kpc2Mpc

        #print(part_df[part_df['Type'] == 4.0].head())
        #print(sum_coord, center, thick)

    # Set the minima and maxima for the particles to be used in the plot
    minima[ax0] = center[ax0] - side * 0.5
    minima[ax1] = center[ax1] - side * 0.5
    minima[ax2] = center[ax2] - thick * 0.5

    maxima[ax0] = center[ax0] + side * 0.5
    maxima[ax1] = center[ax1] + side * 0.5
    maxima[ax2] = center[ax2] + thick * 0.5

    # Find the particles in the slab
    condition_x = (part_df[col0] > minima[ax0]) & (part_df[col0] < maxima[ax0])
    condition_y = (part_df[col1] > minima[ax1]) & (part_df[col1] < maxima[ax1])
    condition_z = (part_df[col2] > minima[ax2]) & (part_df[col2] < maxima[ax2])
    part_select = part_df[(condition_x) & (condition_y) & (condition_z)]

    print('Found: ', len(part_select), ' particles in the slab')

    # Now select a random subsample of the full particle list
    if reduction_factor < 1.0:
        part_select = part_select.sample(frac=reduction_factor, random_state=rand_seed)

        print('The number of particles to be used has been reduced to: ', len(part_select))

    # Return the selected particles' properties in a dataframe
    return part_select


def smooth_web(vweb, x_point=None, smooth_length=1.5, smooth_type='avg'):
    """ vweb is a DataFrame containing all the web information """

    # TODO this needs to be completed
    x_col = ['x', 'y', 'z']
    new_col = 'Distance'

    # Take the simple average of all points within a smoothing_length distance
    if smooth_type == 'avg':
        '''
        vweb[x_col].apply(lambda x: distance(x_point))
        smooth_length
        '''
        pass

    return smooth


def check_units(data=None, cols=None):
    """ Check if the units used are consistent """

    n_pts = int(len(data) * 0.5)

    vals = data[cols].iloc[n_pts]

    # If this is true, then the units are
    if np.sum(vals) < 1.e+4:
        factor = 1.0e+3
    else:
        factor = 1.0

    data[cols] = data[cols].apply(lambda x: x * factor)
    #print(data.head())

    return factor


def ahf_header():
    """ Just in case you were wondering how the header of an AHF file looks like... """

    ahf_header = ['numSubStruct(3)', 'Mvir(4)', 'npart(5)', 'Xc(6)', 'Yc(7)', 'Zc(8)',
       'VXc(9)', 'VYc(10)', 'VZc(11)', 'Rvir(12)', 'Rmax(13)', 'r2(14)',
       'mbp_offset(15)', 'com_offset(16)', 'Vmax(17)', 'v_esc(18)', 'sigV(19)',
       'lambda(20)', 'lambdaE(21)', 'Lx(22)', 'Ly(23)', 'Lz(24)', 'b(25)',
       'c(26)', 'Eax(27)', 'Eay(28)', 'Eaz(29)', 'Ebx(30)', 'Eby(31)',
       'Ebz(32)', 'Ecx(33)', 'Ecy(34)', 'Ecz(35)', 'ovdens(36)', 'nbins(37)',
       'fMhires(38)', 'Ekin(39)', 'Epot(40)', 'SurfP(41)', 'Phi0(42)',
       'cNFW(43)', 'n_gas(44)', 'M_gas(45)', 'lambda_gas(46)',
       'lambdaE_gas(47)', 'Lx_gas(48)', 'Ly_gas(49)', 'Lz_gas(50)',
       'b_gas(51)', 'c_gas(52)', 'Eax_gas(53)', 'Eay_gas(54)', 'Eaz_gas(55)',
       'Ebx_gas(56)', 'Eby_gas(57)', 'Ebz_gas(58)', 'Ecx_gas(59)',
       'Ecy_gas(60)', 'Ecz_gas(61)', 'Ekin_gas(62)', 'Epot_gas(63)',
       'n_star(64)', 'M_star(65)', 'lambda_star(66)', 'lambdaE_star(67)',
       'Lx_star(68)', 'Ly_star(69)', 'Lz_star(70)', 'b_star(71)', 'c_star(72)',
       'Eax_star(73)', 'Eay_star(74)', 'Eaz_star(75)', 'Ebx_star(76)',
       'Eby_star(77)', 'Ebz_star(78)', 'Ecx_star(79)', 'Ecy_star(80)',
       'Ecz_star(81)', 'Ekin_star(82)', 'Epot_star(83)', 'Unnamed: 83', 'ID',
       'HostHalo']

    return ahf_header


def rs_header():
    """ RockStar file header """

    rs_header = ['#ID', 'DescID', 'Mvir', 'Vmax', 'Vrms', 'Rvir', 'Rs', 'Np',
            'X', 'Y', 'Z', 'VX', 'VY', 'VZ', 'JX', 'JY', 'JZ', 'Spin', 'rs_klypin',
            'Mvir_all', 'M200b', 'M200c', 'M500c', 'M2500c', 'Xoff', 'Voff', 'spin_bullock',
            'b_to_a', 'c_to_a', 'A[x]', 'A[y]', 'A[z]', 'b_to_a(500c)', 'c_to_a(500c)',
            'A[x](500c)', 'A[y](500c)', 'A[z](500c)', 'T/|U|', 'M_pe_Behroozi', 'M_pe_Diemer', 'Halfmass_Radius']

    return rs_header


def header_rs2ahf(rs_head):
    """ Convert header from ahf to rockstar """

    ahf_head = ahf_header()
    #rs_head = rs_header()

    rs2ahf = dict()

    rs2ahf['#ID'] = ahf_head[82]
    rs2ahf['Mvir'] = ahf_head[1]
    rs2ahf['Vmax'] = ahf_head[14]
    rs2ahf['Rvir'] = ahf_head[9]
    rs2ahf['Rs'] = ahf_head[11]
    rs2ahf['Np'] = ahf_head[2]
    rs2ahf['X'] = ahf_head[3]
    rs2ahf['Y'] = ahf_head[4]
    rs2ahf['Z'] = ahf_head[5]
    rs2ahf['VX'] = ahf_head[6]
    rs2ahf['VY'] = ahf_head[7]
    rs2ahf['VZ'] = ahf_head[8]
    rs2ahf['JX'] = ahf_head[19]
    rs2ahf['JY'] = ahf_head[20]
    rs2ahf['JZ'] = ahf_head[21]
    rs2ahf['b_to_a'] = ahf_head[22]
    rs2ahf['c_to_a'] = ahf_head[23]
    rs2ahf['spin_bullock'] = ahf_head[18]
    rs2ahf['Spin'] = ahf_head[17]

    new_header = []
    for hd in rs_head:
        if rs2ahf.get(hd) == None:
            rs2ahf[hd] = hd

        new_header.append(rs2ahf[hd])

    return new_header


def periodic_boundaries(data=None, slab_size=10.0, box=100.0, x_cols=['Xc(6)', 'Yc(7)', 'Zc(8)']):
    '''
    Take into account the halos at the boundaries adding them as a sort of buffer to the existing box
    '''

    df_all_tmp = []
    slab_other = box - slab_size
    half_box = 0.5 * box

    # Find all the halos within the buffer regions
    for col in x_cols:
        df_tmp = data[(data[col] < slab_size) | (data[col] > slab_other)]
        df_all_tmp.append(df_tmp)

    # Merge them into a new DF
    df_tmp_new = pd.concat(df_all_tmp)
    #print(len(df_tmp_new))
    df_tmp_new.drop_duplicates(inplace=True)
    #print(len(df_tmp_new))

    def rescale_pos(x):
        if x < slab_size:
            x = float(x + box)
        elif x > slab_other:
            x = float(x - box)

        return x

    # Now correct for the new positions that mirror the periodic boundaries
    for col in x_cols:
        #print(df_tmp_new[col].apply(lambda x: rescale_pos(x)))
        df_tmp_new[col] = df_tmp_new[col].apply(lambda x: rescale_pos(x))

    print(df_tmp_new[x_cols].head(20))

    # Add them all to the old dataframe
    df_new = pd.concat([data, df_tmp_new])

    return df_new


if __name__ == '__main__':
    ''' Use this for testing local functions '''
    pass
