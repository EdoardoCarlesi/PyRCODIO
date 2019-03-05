#!/usr/bin/python

from libcosmo.utils import *
import time
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# RANDOM LGs statistics
rand_stat_mah = 'saved/rand_out_stat_mah.pkl'
rand_stat_time = 'saved/rand_out_stat_time.pkl'

# LGF LGs statistics TODO
cs_stat_mah = 'saved/all_lgf_stat_mah.pkl'
cs_stat_time = 'saved/all_lgf_stat_time.pkl'
	
# Filenames
print('Loading MAH random LCDM statistics from %s and %s' % (rand_stat_mah, rand_stat_time))
print('Loading MAH CS LCDM statistics from %s and %s' % (cs_stat_mah, cs_stat_time))

# Open files
rand_f_mah = open(rand_stat_mah, 'r')
rand_f_time = open(rand_stat_time, 'r')
cs_f_mah = open(cs_stat_mah, 'r')
cs_f_time = open(cs_stat_time, 'r')

# Load all files
rand_times = pickle.load(rand_f_time)
rand_mahs = pickle.load(rand_f_mah)
cs_times = pickle.load(cs_f_time)
cs_mahs = pickle.load(cs_f_mah)

'''
	GENERAL PLOT PROPERTIES
'''
plt.rc({'text.usetex': True})



'''
	2D PLOTS OF LGs MMT & FT
'''

#	time_stats = np.zeros((2, 2, iValid))	---> # LMMT & FT
n_bins_rand = 18
n_bins_cs = 12

# FT bins 2D ---> Rand
ft_fb = rand_times[:, 0, :]
plt.hist2d(ft_fb[0], ft_fb[1], bins=(n_bins_rand, n_bins_rand), cmap=plt.cm.jet, alpha=0.8)

# Contour lines ---> CS
ft_cs = cs_times[:, 0, :]
H, xed, yed = np.histogram2d(ft_cs[1], ft_cs[0], bins=(n_bins_cs, n_bins_cs))
plt.contour(H, colors='black')

# Axes and title
ax = plt.gca()
ax.set_xlabel('MW')
ax.set_ylabel('M31')
plt.title('Formation Time')

# Formation Time 2D
plt.tight_layout()
plt.savefig('histo_ft.png')
plt.clf()
plt.cla()

# MMT bins 2D
mmt_fb = rand_times[:, 1, :]
plt.hist2d(mmt_fb[0], mmt_fb[1], bins=(n_bins_rand, n_bins_rand), cmap=plt.cm.jet)

# Contour lines ---> CS
mmt_cs = cs_times[:, 1, :]
H, xed, yed = np.histogram2d(mmt_cs[1], mmt_cs[0], bins=(n_bins_cs, n_bins_cs))
plt.contour(H, colors='black')
#contour(counts,extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],linewidths=3)

# Axes and title
ax = plt.gca()
ax.set_xlabel('MW')
ax.set_ylabel('M31')
plt.title('Last Major Merger Time')

plt.tight_layout()
plt.savefig('histo_mmt.png')
plt.clf()
plt.cla()

'''
	1D PLOTS OF MAHs
'''
#	mah_avgs = np.zeros((2, n_steps, 3))	---> This contains median and percentiles 

# General properties
x_min = 0.0; x_max = 13.0
y_min = 0.01; y_max = 1.1
n_snaps = 54; t_step = 0.25

# Set the x axis
time = np.zeros((n_snaps))
for i in range(0, n_snaps):
	time[i] = i * t_step

plt.xlabel('GYr'); plt.ylabel('M')

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(time, rand_mahs[0, :, 1], color='red')
axs.fill_between(time, rand_mahs[0, :, 0], rand_mahs[0, :, 2], facecolor='orange', alpha=0.2)

# Overplot
axs.plot(time, cs_mahs[0, :, 1], color='blue')
axs.fill_between(time, cs_mahs[0, :, 0], cs_mahs[0, :, 2], facecolor='cyan', alpha=0.2)

# Axes and title
ax = plt.gca()
ax.set_xlabel('T [Gyr]')
ax.set_ylabel('M/M(z=0)')
plt.title('MW MAH')
plt.tight_layout()
plt.savefig('mah_mw.png')
plt.clf()
plt.cla()

'''
	1D HISTOGRAMS
'''
n_bins = 10
plt.hist(ft_fb[0], n_bins, facecolor='blue', normed=1, alpha=0.5)
plt.hist(ft_fb[1], n_bins, facecolor='red', normed=1, alpha=0.75)
plt.tight_layout()
plt.savefig('test_mah_ft_histo.png')
plt.clf()
plt.cla()



'''
	K-S TEST AND CUMULATIVE DISTRIBUTION
'''

#ft_fb[0]
#ft_fb[1]

#f0 = np.sort(ft_fb[0])
#f1 = np.sort(ft_fb[1])

cdf0 = cdf(ft_fb[0])
cdf1 = cdf(ft_fb[1])

kst = stats.ks_2samp(cdf0[0], cdf1[0])
#print(kst)
#print(f0, f1)
#kst = stats.kstest(f0, 'norm') #f1)
#print(kst)
#stats.kstest(ft_fb[0], ft_fb[1])

axs.axis([0.0, 13.0, 0.0, 1.0])

plt.xlabel('GYr'); plt.ylabel('M')
(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))
axs.plot(cdf0[0], cdf0[1], color='red')
axs.plot(cdf1[0], cdf1[1], color='blue')
plt.tight_layout()
plt.savefig('test_mah_kstest.png')
plt.clf()
plt.cla()




