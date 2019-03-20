from libcosmo.utils import *
import time
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# RANDOM LGs statistics
#rand_stat_mah = 'saved/rand_out_stat_mah.pkl'
#rand_stat_time = 'saved/rand_out_stat_time.pkl'
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
cs_times = pickle.load(rand_f_time)
cs_mahs = pickle.load(rand_f_mah)
rand_times = pickle.load(cs_f_time)
rand_mahs = pickle.load(cs_f_mah)
'''


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
out_dir = 'output/merger_trees/'
file_ft_mw_m31 = out_dir + '2d_ft.png'
file_mmt_mw_m31 = out_dir + '2d_mmt.png'

file_m31_mah = out_dir + 'mah_m31.png'
file_m31_ft_hist = out_dir + 'histo_ft_m31.png'
file_m31_mmt_hist = out_dir + 'histo_mmt_m31.png'

file_mw_mah = out_dir + 'mah_mw.png'
file_mw_ft_hist = out_dir + 'histo_ft_mw.png'
file_mw_mmt_hist = out_dir + 'histo_mmt_mw.png'

file_mw_mmt_kstest = out_dir + 'kstest_mmt_mw.png'
file_mw_ft_kstest = out_dir + 'kstest_ft_mw.png'
file_m31_mmt_kstest = out_dir + 'kstest_mmt_m31.png'
file_m31_ft_kstest = out_dir + 'kstest_ft_m31.png'

'''
	2D PLOTS OF LGs MMT & FT
'''

# DATA FORMAT:	time_stats = np.zeros((2, 2, iValid))	---> # LMMT & FT
n_bins_rand = 15
n_bins_cs = 13
opaque = 0.7
linewds = [2,2,2]
#linecol = 'black'
linecol = 'red'
#linecol = 'blue'
n_levels = [3.0, 5.0, 7.0]
#color_map = plt.cm.Reds
color_map = plt.cm.Blues
#color_map = plt.cm.Greys
x_min = 0.0; x_max = 12.0
y_min = 0.0; y_max = 12.0

# Datasets
ft_fb = rand_times[:, 0, :]
ft_cs = cs_times[:, 0, :]

# FT bins 2D ---> Rand
plt.hist2d(ft_fb[mw_fb_id], ft_fb[m31_fb_id], bins=(n_bins_rand, n_bins_rand), cmap=color_map, alpha=1.0)

# Contour lines ---> CS
H, xed, yed = np.histogram2d(ft_cs[m31_cs_id], ft_cs[mw_cs_id], bins=(n_bins_cs, n_bins_cs))
#H, xed, yed = np.histogram2d(ft_cs[mw_cs_id], ft_cs[m31_cs_id], bins=(n_bins_cs, n_bins_cs))
plt.contour(H, colors=linecol, linewidths=linewds, levels=n_levels)

# Axes and title
ax = plt.gca()
ax.axis([x_min, x_max, y_min, y_max])
ax.set_xlabel('$\\tau _F$ MW')
ax.set_ylabel('$\\tau _F$ M31')
#plt.title('Formation Time')

# Formation Time 2D
plt.tight_layout()
plt.savefig(file_ft_mw_m31)
plt.clf()
plt.cla()


# MMT bins 2D
#n_levels = [2.5, 4.5, 6.5]
mmt_fb = rand_times[:, 1, :]
mmt_cs = cs_times[:, 1, :]
plt.hist2d(mmt_fb[mw_fb_id], mmt_fb[m31_fb_id], bins=(n_bins_rand, n_bins_rand), cmap=color_map, alpha = 1.0)

# Contour lines ---> CS
#H, xed, yed = np.histogram2d(mmt_cs[m31_cs_id], mmt_cs[mw_cs_id], bins=(n_bins_cs, n_bins_cs))
H, xed, yed = np.histogram2d(mmt_cs[mw_cs_id], mmt_cs[m31_cs_id], bins=(n_bins_cs, n_bins_cs))
plt.contour(H, colors=linecol, linewidths=linewds, levels = n_levels)

# Axes and title
ax = plt.gca()
ax.axis([x_min, x_max, y_min, y_max])
#ax.set_xlabel('MW')
#ax.set_ylabel('M31')
ax.set_xlabel('$\\tau _M$ MW')
ax.set_ylabel('$\\tau _M$ M31')
#plt.title('Last Major Merger Time')

plt.tight_layout()
plt.savefig(file_mmt_mw_m31)
plt.clf()
plt.cla()

######################  1D PLOTS OF MAHs
#	mah_avgs = np.zeros((2, n_steps, 3))	---> This contains median and percentiles 

# General properties
x_min = 0.0; x_max = 13.0
y_min = 0.01; y_max = 1.1
n_snaps = 54; t_step = 0.25

# Set the x axis
time = np.zeros((n_snaps))
for i in range(0, n_snaps):
	time[i] = i * t_step

(fig, axs) = plt.subplots(ncols=2, nrows=1, figsize=(8, 4))

axs[0].axis([x_min, x_max, y_min, y_max])
axs[0].plot(time, rand_mahs[mw_fb_id, :, 1], color=col_fb)
axs[0].fill_between(time, rand_mahs[mw_fb_id, :, 0], rand_mahs[mw_fb_id, :, 2], facecolor=col_fb, alpha=0.2)
axs[0].plot(time, cs_mahs[mw_cs_id, :, 1], color=col_cs)
axs[0].fill_between(time, cs_mahs[mw_cs_id, :, 0], cs_mahs[mw_cs_id, :, 2], facecolor=col_cs, alpha=0.2)

# Axes and title
#axs[0].set_title('MW')
axs[0].set_xlabel('t [Gyr]')
#axs[0].set_ylabel('M/M(z=0)')
axs[0].set_ylabel('$M(z) / M(z=0)$')

#plt.subplot(2, 1, 1)
#axs[1].yaxis.set_ticks_position('none') 
#plt.rc('ytick', labelsize=0)
plt.yticks([])
axs[1].axis([x_min, x_max, y_min, y_max])
axs[1].plot(time, rand_mahs[m31_fb_id, :, 1], color=col_fb)
axs[1].fill_between(time, rand_mahs[m31_fb_id, :, 0], rand_mahs[m31_fb_id, :, 2], facecolor=col_fb, alpha=0.2)
axs[1].plot(time, cs_mahs[m31_cs_id, :, 1], color=col_cs)
axs[1].fill_between(time, cs_mahs[m31_cs_id, :, 0], cs_mahs[m31_cs_id, :, 2], facecolor=col_cs, alpha=0.2)

# Axes and title
#axs[1].set_title('M31')
axs[1].set_xlabel('t [Gyr]')

plt.tight_layout()
plt.savefig(file_m31_mah)
plt.clf()
plt.cla()


############################   1D HISTOGRAMS
n_bins = n_bins_rand
(fig, ax) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

ticks = [0.0, 0.1, 0.2, 0.3, 0.4]

#percent = [25, 50, 75]
percent = [32, 64, 80]

x_min = 0; x_max = 12
y_min = 0; y_max = 0.41

################### FORMATION TIMES #########################

# MW 
plt.hist(ft_fb[mw_fb_id], n_bins, facecolor=col_fb, normed=1, alpha=1.0)
plt.hist(ft_cs[mw_cs_id], n_bins, facecolor=col_cs, normed=1, alpha=opaque)
#plt.title('MW Formation Times')

ax = plt.gca()
ax.set_xlabel('t [Gyr]')
ax.set_ylabel('$n(\\tau _F)$')
ax.axis([x_min, x_max, y_min, y_max])
#ax.set_ylabel('$n( < \\tau _F)$')
plt.yticks(ticks)
plt.tight_layout()
plt.savefig(file_mw_ft_hist) 
plt.clf()
plt.cla()

# M31
plt.hist(ft_fb[m31_fb_id], n_bins, facecolor=col_fb, normed=1, alpha=1.0)
plt.hist(ft_cs[m31_cs_id], n_bins, facecolor=col_cs, normed=1, alpha=opaque)
#plt.title('M31 Formation Times')
ax = plt.gca()
ax.axis([x_min, x_max, y_min, y_max])
ax.set_xlabel('t [Gyr]')
#ax.set_ylabel('n(T)')
plt.yticks([]) 
plt.tight_layout()
plt.savefig(file_m31_ft_hist) 
plt.clf()
plt.cla()

################### LAST MAJOR MERGER TIMES #########################

y_min = 0; y_max = 0.2
ticks = [0.0, 0.1, 0.2]

# MW 
plt.hist(mmt_fb[mw_fb_id], n_bins, facecolor=col_fb, normed=1, alpha=1.0)
plt.hist(mmt_cs[mw_cs_id], n_bins, facecolor=col_cs, normed=1, alpha=opaque)
#plt.title('MW Last Major Merger Times')
ax = plt.gca()
ax.set_xlabel('t [Gyr]')
ax.set_ylabel('$n( < \\tau _M)$')
ax.axis([x_min, x_max, y_min, y_max])
plt.yticks(ticks) 
plt.tight_layout()
plt.savefig(file_mw_mmt_hist) 
plt.clf()
plt.cla()

# M31
plt.hist(mmt_fb[m31_fb_id], n_bins, facecolor=col_fb, normed=1, alpha=1.0)
plt.hist(mmt_cs[m31_cs_id], n_bins, facecolor=col_cs, normed=1, alpha=opaque)
#plt.title('M31 Last Major Merger Times')
ax.axis([x_min, x_max, y_min, y_max])
plt.yticks([]) 
ax = plt.gca()
ax.set_xlabel('t [Gyr]')
#ax.set_ylabel('n(T)')
plt.tight_layout()
plt.savefig(file_m31_mmt_hist) 
plt.clf()
plt.cla()

print('Percentiles, MT, MW: ')
print(np.percentile(mmt_fb[mw_fb_id], percent))
print(np.percentile(mmt_cs[mw_cs_id], percent))
print('M31: ')
print(np.percentile(mmt_fb[m31_fb_id], percent))
print(np.percentile(mmt_cs[m31_cs_id], percent))
print('Percentiles, FT, MW: ')
print(np.percentile(ft_fb[mw_fb_id], percent))
print(np.percentile(ft_cs[mw_cs_id], percent))
print('M31: ')
print(np.percentile(ft_fb[m31_fb_id], percent))
print(np.percentile(ft_cs[m31_cs_id], percent))



#######	K-S TEST AND CUMULATIVE DISTRIBUTION #############

y_min = 0; y_max = 1.01
ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# LAST MAJOR MERGER

cdf_fb = cdf(mmt_fb[mw_fb_id])
cdf_cs = cdf(mmt_cs[mw_cs_id])
kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)

axs.set_xlabel('t [Gyr]')
axs.set_ylabel('$n(< \\tau _M)$')
#plt.title('MW MMT Cumulative')
plt.yticks(ticks) 
plt.tight_layout()
plt.savefig(file_mw_mmt_kstest)
plt.clf()
plt.cla()

cdf_fb = cdf(mmt_fb[m31_fb_id])
cdf_cs = cdf(mmt_cs[m31_cs_id])

kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)

axs.set_xlabel('t [Gyr]')
#axs.set_ylabel('n(<T)$')
#axs.set_ylabel('$n(< \\tau _M)$')
#plt.title('M31 MMT Cumulative')
plt.tight_layout()
plt.savefig(file_m31_mmt_kstest)
plt.yticks([]) 
plt.clf()
plt.cla()

# FORMATION TIME

cdf_fb = cdf(ft_fb[mw_fb_id])
cdf_cs = cdf(ft_cs[mw_cs_id])
kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)

axs.set_xlabel('t [Gyr]')
axs.set_ylabel('$n(< \\tau _F)$')
#axs.set_ylabel('n(<T)$')
#plt.title('MW FT Cumulative')
plt.yticks(ticks) 
plt.tight_layout()
plt.savefig(file_mw_ft_kstest)
plt.clf()
plt.cla()

cdf_fb = cdf(ft_fb[m31_fb_id])
cdf_cs = cdf(ft_cs[m31_cs_id])

kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)

axs.set_xlabel('t [Gyr]')
#axs.set_ylabel('T_{F}(<T)$')
#axs.set_ylabel('$n(< \\tau _F)$')
#plt.title('M31 FT Cumulative')
plt.yticks([]) 
plt.tight_layout()
plt.savefig(file_m31_ft_kstest)
plt.clf()
plt.cla()









