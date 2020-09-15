import scipy as sp
import numpy as np
import units as us
import pandas as pd
import read_files as rf
import dask.dataframe as dd


def sincos(theta=None, radians=False):

    if radians == False:
        theta = np.deg2rad(theta)

    sc = (np.sin(theta) * (theta - np.sin(theta))) / np.power((1.0 - np.cos(theta)), 2.0)

    return sc

def vt_r(v=None, t=None, r=None):

    '''
        Velocity [km/s]
        Time [Gyr]
        R [kpc]
    '''

    v = v * us.s2Myr() / 1.e+3 / us.km2kpc()

    vtr = (v * t) / r

    return vtr

def Rmax(r=None, theta=None, radians=False):

    if radians == False:
        theta = np.deg2rad(theta)

    return 2 * r / (1.0 - np.cos(theta))

def mass_estimate(r=710, v=-118, t=13.75, n_pts=5000):
    
    GN = us.G()
    n_pts = n_pts
    delta = 360.0 / n_pts
    th = np.zeros((n_pts))
    sc = np.zeros((n_pts))

    vtr = vt_r(v=v, t=t, r=r)

    for i in range(0, n_pts):
        th[i] = i * delta
        sc[i] = sincos(i * delta) - vtr

    # Take the first point after the change of sign to find out the zero of the function (approximate)
    try:
        idx = np.where(sc < 0.0)
        id_cross = idx[0][0]
        theta = th[id_cross]

        # Rescale to Mpc from kpc units
        Rm = Rmax(r=r, theta=theta)/1.0e+3
    
        # Always convert to radians
        theta = np.deg2rad(theta)
        Mtot = ((Rm ** 3.0) * ((theta - np.sin(theta)) ** 2.0)) / (GN * 8.0 * (t ** 2.0))
    
    # If for some reason the equation has no solution then set a negative mass as a flag
    except:
        Mtot = -10.0

    return Mtot

# TODO: correct the wrong masses extrapolated with the Timing Argument
def wrong_ta_correction(data=None):

    data['M_TA']

    fac_max = 1.1
    fac_min = 0.9

    return data



'''
    ======= MAIN PROGRAM
    Read files and add timing argument column
'''

import seaborn as sns
import matplotlib.pyplot as plt

vrad_max = -1.0
#lg_df = rf.read_lg_lgf(); file_out = 'output/lg_pairs_512_TA.csv'
#lg_df = rf.read_lg_fullbox(); file_out = 'output/lg_fullbox_TA.csv'
#lg_df = rf.read_lg_rs_fullbox([0, 10]); file_out = 'output/lg_fullbox_rs_0010_TA.csv'

#lg_df = dd.from_pandas(lg_df_pd, npartitions=3)
#lg_df = lg_df[lg_df['Vrad'] < vrad_max]
#lg_df = pd.read_csv('output/lg_pairs_2048_TA.csv')

#lg_df['M_TA'] = lg_df.apply(lambda x: mass_estimate(r=x['R'], v=x['Vrad']), axis=1)
#file_out = 'output/lg_fullbox_rs_TA_0010.csv'
#lg_df.to_csv(file_out)


#lg_df = rf.read_lg_lgf(TA=True); out_file='output/ta_vs_true_lgf.png'
#lg_df = rf.read_lg_rs_fullbox(TA=True); out_file='output/ta_vs_true_rs.png'
lg_df = rf.read_lg_fullbox(TA=True); out_file='output/ta_vs_true_rs.png'
lg_df = lg_df[lg_df['Vrad'] < vrad_max]
lg_df['Mtot'] = lg_df['M_MW'] + lg_df['M_M31']
ta_wrong = lg_df[lg_df['M_TA'] < 0.0]
mass_norm = 1.0e+12

sns.scatterplot(np.log10(lg_df['Mtot']/mass_norm), np.log10(lg_df['M_TA']/mass_norm))
sns.kdeplot(np.log10(lg_df['Mtot']/mass_norm), np.log10(lg_df['M_TA']/mass_norm), n_levels=4)

slopes = np.polyfit(np.log10(lg_df['Mtot']/mass_norm), np.log10(lg_df['M_TA']/mass_norm), 1)

print(slopes)

title = 'Sim vs. TA Mass. Slope: ' + '%.3f' % slopes[0]
plt.title(title)

plt.savefig(out_file)

#plt.show()
'''
'''







