from libcosmo.utils import *
from libcosmo.units import *
from matplotlib import gridspec
import time
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

# RANDOM LGs statistics
rand_stat_mah = 'saved/rand_halo_stat_mah.pkl'
rand_stat_time = 'saved/rand_halo_stat_time.pkl'

# LGF LGs statistics TODO
cs_stat_mah = 'saved/all_lgf_stat_mah.pkl'
cs_stat_time = 'saved/all_lgf_stat_time.pkl'

# Filenames
print('Loading MAH random LCDM statistics from %s and %s' % (rand_stat_mah, rand_stat_time))
print('Loading MAH CS LCDM statistics from %s and %s' % (cs_stat_mah, cs_stat_time))

# Open files
rand_f_mah = open(rand_stat_mah, 'rb')
rand_f_time = open(rand_stat_time, 'rb')
cs_f_mah = open(cs_stat_mah, 'rb')
cs_f_time = open(cs_stat_time, 'rb')

# Load all files
rand_times = pickle.load(rand_f_time)
rand_mahs = pickle.load(rand_f_mah)
cs_times = pickle.load(cs_f_time)
cs_mahs = pickle.load(cs_f_mah)


'''
        GENERAL PLOT PROPERTIES
'''
plt.rc({'text.usetex': True})

mw_cs_id = 0
m31_cs_id = 1
mw_fb_id = 1
m31_fb_id = 0

col_cs = 'red'
col_fb = 'blue'

# Filenames
file_vmax_mah = 'mah_vmax.png'

# General properties
x_min = 0.0; x_max = 12.0
y_min = 0.25; y_max = 1.05
n_snaps = 54; t_step = 0.25

# Set the x axis
time = np.zeros((n_snaps))
for i in range(0, n_snaps):
    time[i] = i * t_step

rand_vmax = np.zeros((2, n_snaps, 3))
cs_vmax = np.zeros((2, n_snaps, 3))

# Set some parameters
m0 = 1.e+12
slope=0.3361033
inter=-1.713984

this_lg10_vmax0 = inter + slope * np.log10(m0 * rand_mahs[0, 0, 1])

for i in range(0, 2):
    for k in range(0, 3):
        for j in range(0, n_snaps):
            rand_lg10_vmax = inter + slope * np.log10(m0 * rand_mahs[i, j, k])
            cs_lg10_vmax = inter + slope * np.log10(m0 * cs_mahs[i, j, k])
            rand_vmax[i, j, k] = np.power(10, rand_lg10_vmax) / np.power(10, this_lg10_vmax0)
            cs_vmax[i, j, k] = np.power(10, cs_lg10_vmax) / np.power(10, this_lg10_vmax0)
            #print(i, j, k, rand_vmax[i, j, k])

plt.figure(0)
axs0 = plt.subplot2grid((5, 8), (0, 0), rowspan=4, colspan=4)
axs1 = plt.subplot2grid((5, 8), (0, 4), rowspan=4, colspan=4)
axs2 = plt.subplot2grid((5, 8), (4, 0), rowspan=1, colspan=4)
axs3 = plt.subplot2grid((5, 8), (4, 4), rowspan=1, colspan=4)

axs0.axis([x_min, x_max, y_min, y_max])
axs0.plot(time, rand_vmax[mw_fb_id, :, 1], color=col_fb)
axs0.fill_between(time, rand_vmax[mw_fb_id, :, 0], rand_vmax[mw_fb_id, :, 2], facecolor=col_fb, alpha=0.2)
axs0.plot(time, cs_vmax[mw_cs_id, :, 1], color=col_cs)
axs0.fill_between(time, cs_vmax[mw_cs_id, :, 0], cs_vmax[mw_cs_id, :, 2], facecolor=col_cs, alpha=0.2)

# Axes and title
axs0.set_title('MW')
axs0.set_ylabel('$V_{MAX}(t) / V_{MAX}(t_0)$')

axs0.set_xticks([])
axs1.set_xticks([])
axs1.set_yticks([])
axs1.axis([x_min, x_max, y_min, y_max])
axs1.plot(time, rand_vmax[m31_fb_id, :, 1], color=col_fb)
axs1.fill_between(time, rand_vmax[m31_fb_id, :, 0], rand_vmax[m31_fb_id, :, 2], facecolor=col_fb, alpha=0.2)
axs1.plot(time, cs_vmax[m31_cs_id, :, 1], color=col_cs)
axs1.fill_between(time, cs_vmax[m31_cs_id, :, 0], cs_vmax[m31_cs_id, :, 2], facecolor=col_cs, alpha=0.2)

# Axes and title
axs1.set_title('M31')

# Plot ratios
ratio_m31 = []; ratio_mw = []; unity = []
for i in range(0, len(time)):
    v_rand = rand_vmax[m31_fb_id, i, 1]
    v_cs = cs_vmax[m31_fb_id, i, 1]
    ratio_m31.append(v_cs / v_rand)
    unity.append(1.0)
    v_rand = rand_vmax[mw_fb_id, i, 1]
    v_cs = cs_vmax[mw_fb_id, i, 1]
    ratio_mw.append(v_cs / v_rand)

axs2.set_xlabel('t [Gyr]')
axs3.set_xlabel('t [Gyr]')
axs2.axis([x_min, x_max, 0.99, 1.75])
axs2.plot(time, ratio_m31, color='black')
axs2.set_yticks([1.0,1.25,1.5,1.75])
axs3.set_yticks([])
axs3.axis([x_min, x_max, 0.99, 1.5])
axs3.plot(time, ratio_mw, color='black')
plt.tight_layout()
plt.subplots_adjust(hspace=0.0, wspace=0.2) 
plt.savefig(file_vmax_mah)
plt.clf()
plt.cla()

'''
'''
