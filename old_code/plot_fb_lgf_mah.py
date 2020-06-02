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

# z values
z_file='/home/eduardo/CLUES/DATA/output/output_z.txt'

# RANDOM LGs statistics
#rand_stat_mah = 'saved/rand_out_stat_mah.pkl'
#rand_stat_time = 'saved/rand_out_stat_time.pkl'
rand_stat_mah = 'saved/rand_halo_stat_mah.pkl'
rand_stat_time = 'saved/rand_halo_stat_time.pkl'

# Full Trees
m31_trees = 'saved/all_trees_lgf_m31.pkl'
mw_trees = 'saved/all_trees_lgf_mw.pkl'

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

#file_m31_mah = 'ProjectsPlots/mah_m31.png'
#file_m31_mah_ratio = 'ProjectsPlots/mah_m31_ratio.png'
file_m31_mah = 'mah_m31.png'
file_m31_mah_ratio = 'mah_ratio.png'
file_m31_ft_hist = out_dir + 'histo_ft_m31.png'
file_m31_mmt_hist = out_dir + 'histo_mmt_m31.png'

file_mw_mah = 'ProjectsPlots/mah_mw.png'
file_mw_mah_ratio = 'ProjectsPlots/mah_mw_ratio.png'
file_mw_ft_hist = out_dir + 'histo_ft_mw.png'
file_mw_mmt_hist = out_dir + 'histo_mmt_mw.png'

file_mw_mmt_kstest = out_dir + 'kstest_mmt_mw.png'
file_mw_ft_kstest = out_dir + 'kstest_ft_mw.png'
file_m31_mmt_kstest = out_dir + 'kstest_mmt_m31.png'
file_m31_ft_kstest = out_dir + 'kstest_ft_m31.png'
'''
        2D PLOTS OF LGs MMT & FT
'''

# DATA FORMAT:  time_stats = np.zeros((2, 2, iValid))   ---> # LMMT & FT
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

'''
'''

# Full trees
#m31_trees = 'all_trees_lgf_m31.pkl'
#mw_trees = 'all_trees_lgf_mw.pkl'

f_mt_mw = open(mw_trees, 'rb')
f_mt_m31 = open(m31_trees, 'rb')

mt_mw = pickle.load(f_mt_mw)
mt_m31 = pickle.load(f_mt_m31)



######################  1D PLOTS OF MAHs
#       mah_avgs = np.zeros((2, n_steps, 3))    ---> This contains median and percentiles

# General properties
x_min = 0.0; x_max = 12.0
y_min = 0.01; y_max = 1.1
n_snaps = 54; t_step = 0.25

size_lab = 13

# Set the x axis
time = np.zeros((n_snaps))
for i in range(0, n_snaps):
    time[i] = i * t_step

plt.figure(0)

axs0 = plt.subplot2grid((5, 8), (0, 0), rowspan=4, colspan=4)
axs1 = plt.subplot2grid((5, 8), (0, 4), rowspan=4, colspan=4)
axs2 = plt.subplot2grid((5, 8), (4, 0), rowspan=1, colspan=4)
axs3 = plt.subplot2grid((5, 8), (4, 4), rowspan=1, colspan=4)

axs0.axis([x_min, x_max, y_min, y_max])
axs0.plot(time, rand_mahs[mw_fb_id, :, 1], color=col_fb)
axs0.fill_between(time, rand_mahs[mw_fb_id, :, 0], rand_mahs[mw_fb_id, :, 2], facecolor=col_fb, alpha=0.2)
axs0.plot(time, cs_mahs[mw_cs_id, :, 1], color=col_cs)
axs0.fill_between(time, cs_mahs[mw_cs_id, :, 0], cs_mahs[mw_cs_id, :, 2], facecolor=col_cs, alpha=0.2)

imm=0; ids = []
for imw in range(0, 31):
    mah_mw = mt_mw[imw].smooth_tree()
    v0 = 954
    v1 = 970
    intv = int(mah_mw[1] * 1000)

'''
    if intv != v0 and intv != v1: 
        print(mah_mw[1])
        ids.append(imw)
        imm = imm +1
        axs0.plot(time, mah_mw)
print(imm)
'''

#print(mt_mw[0].smooth_tree())
# Axes and title
axs0.set_title('MW')
axs0.set_ylabel('$M(t) / M(t_0)$', size=size_lab)

axs0.set_xticks([])
axs1.set_xticks([])
axs1.set_yticks([])
axs1.axis([x_min, x_max, y_min, y_max])
axs1.plot(time, rand_mahs[m31_fb_id, :, 1], color=col_fb)
axs1.fill_between(time, rand_mahs[m31_fb_id, :, 0], rand_mahs[m31_fb_id, :, 2], facecolor=col_fb, alpha=0.2)
axs1.plot(time, cs_mahs[m31_cs_id, :, 1], color=col_cs)
axs1.fill_between(time, cs_mahs[m31_cs_id, :, 0], cs_mahs[m31_cs_id, :, 2], facecolor=col_cs, alpha=0.2)

'''
imm=0
for im31 in ids:
    mah_m31 = mt_m31[im31].smooth_tree()
    axs1.plot(time, mah_m31)
'''

# Axes and title
axs1.set_title('M31')

# Plot ratios
ratio_m31 = []; ratio_mw = []; unity = []
for i in range(0, len(time)):
    v_rand = rand_mahs[m31_fb_id, i, 1]
    v_cs = cs_mahs[m31_fb_id, i, 1]
    ratio_m31.append(1.0 + 0.5 * abs(v_cs / v_rand - 1.0))

    v_rand = rand_mahs[mw_fb_id, i, 1]
    v_cs = cs_mahs[mw_fb_id, i, 1]
    ratio_mw.append(1.0 + 0.25 * abs(1.0 -v_cs / v_rand))
    unity.append(1.0)

axs2.set_xlabel('t [Gyr]', size=size_lab)
axs3.set_xlabel('t [Gyr]', size=size_lab)
axs2.axis([x_min, x_max, 0.99, 1.75])
axs2.plot(time, ratio_m31, color='black')
#axs2.plot(time, unity, color='red')
axs2.set_yticks([1.0,1.25,1.5,1.75])
axs3.set_yticks([])
axs3.axis([x_min, x_max, 0.99, 1.5])
axs3.plot(time, ratio_mw, color='black')
#axs3.plot(time, unity, color='red')
plt.tight_layout()
plt.subplots_adjust(hspace=0.0, wspace=0.2) 
plt.savefig(file_m31_mah_ratio)
plt.clf()
plt.cla()

'''
a0 = 5.7
b0 = 2.7

z_time_inv = []; z_time = []
z_in = open(z_file, 'r')
data = z_in.read()
for dat in data.split():
    z_time_inv.append(float(dat))
nz = len(z_time_inv)
for it in range(0, nz):
    zt = z_time_inv[nz - it - 1]
    z_time.append(zt)


print(z_time)
output1 = curve_fit(mahFit, z_time, cs_mahs[mw_cs_id, :, 1], p0 = [a0, b0]) #, sigma = mw_cs_err)
output2 = curve_fit(mahFit, z_time, cs_mahs[mw_cs_id, :, 0], p0 = [a0, b0]) #, sigma = mw_cs_err)
output3 = curve_fit(mahFit, z_time, cs_mahs[mw_cs_id, :, 2], p0 = [a0, b0]) #, sigma = mw_cs_err)
print('1. Fit MW CS:', output1)
print('2. Fit MW CS:', output2)
print('3. Fit MW CS:', output3)
output1 = curve_fit(mahFit, z_time, rand_mahs[mw_fb_id, :, 1], p0 = [a0, b0]) #, sigma = mw_cs_err)
output2 = curve_fit(mahFit, z_time, rand_mahs[mw_fb_id, :, 0], p0 = [a0, b0]) #, sigma = mw_cs_err)
output3 = curve_fit(mahFit, z_time, rand_mahs[mw_fb_id, :, 2], p0 = [a0, b0]) #, sigma = mw_cs_err)
print('1. Fit MW FB:', output1)
print('2. Fit MW FB:', output2)
print('3. Fit MW FB:', output3)

print(z_time)
output = curve_fit(mahFit, z_time, cs_mahs[m31_cs_id, :, 1], p0 = [a0, b0]) #, sigma = mw_cs_err)
print('Fit M31 CS:', output)
output = curve_fit(mahFit, z_time, rand_mahs[m31_fb_id, :, 1], p0 = [a0, b0]) #, sigma = mw_cs_err)
print('Fit M31 FB:', output)

'''




###########################   1D HISTOGRAMS
n_bins = n_bins_rand
(fig, ax) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

ticks = [0.0, 0.1, 0.2, 0.3, 0.4]

percent = [25, 50, 75]
#percent = [32, 64, 80]

x_min = 0; x_max = 12.5
y_min = 0; y_max = 0.41

size_lab = 18
size_lab1 = 15

################### FORMATION TIMES #########################

# MW
plt.hist(ft_fb[mw_fb_id], n_bins, facecolor=col_fb, normed=1, alpha=1.0)
plt.hist(ft_cs[mw_cs_id], n_bins, facecolor=col_cs, normed=1, alpha=opaque)
#plt.title('MW 3ormation Times')

ax = plt.gca()
ax.set_xlabel(' ')
#ax.set_xlabel('t [Gyr]')
ax.set_ylabel('$n(\\tau _F)$', size = size_lab)
ax.axis([x_min, x_max, y_min, y_max])
#ax.set_ylabel('$n( < \\tau _F)$')
plt.title('Formation time MW', size = size_lab1)
plt.xticks([])
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
ax.set_xlabel(' ')
#ax.set_xlabel('t [Gyr]')
#ax.set_ylabel('n(T)')
plt.title('Formation time M31', size = size_lab1)
plt.xticks([])
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
#ax.set_xlabel('t [Gyr]')
ax.set_xlabel(' ')
ax.set_ylabel('$n( < \\tau _M)$', size=size_lab)
ax.axis([x_min, x_max, y_min, y_max])
plt.title('Last major merger time MW', size = size_lab1)
plt.xticks([])
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
plt.xticks([])
plt.yticks([])
ax = plt.gca()
#ax.set_xlabel('t [Gyr]')
ax.set_xlabel(' ')
#ax.set_ylabel(' ')
plt.title('Last major merger time M31', size=size_lab1)
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



####### K-S TEST AND CUMULATIVE DISTRIBUTION #############

y_min = 0; y_max = 1.01
ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
ksnorm = 314.0

# LAST MAJOR MERGER
cdf_fb = cdf(mmt_fb[mw_fb_id])
cdf_cs = cdf(mmt_cs[mw_cs_id])

KS_randSubSample(mmt_cs[mw_cs_id], mmt_fb[mw_fb_id], 10000)
KS_randBootstrap(mmt_fb[mw_fb_id], 314, 10000)

print('=== MMT MW  ===')
pos05 = np.where(cdf_cs[0] == 6.75)
print('CS, half: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 6.75)
print('FB, half: ', 1.0 - cdf_fb[1][pos05[0][0]])
pos05 = np.where(cdf_cs[0] == 10)
print('CS, quart: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 10)
print('FB, quart: ', 1.0 - cdf_fb[1][pos05[0][0]])


kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)
print('KST confidence level = ', kst[0] * 16.48)
print('KST confidence level = ', kst[1] * ksnorm)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)

axs.set_xlabel('t [Gyr]', size  = size_lab)
axs.set_ylabel('$n(< \\tau _M)$', size = size_lab)
#plt.title('MW MMT Cumulative')
plt.yticks(ticks)
plt.tight_layout()
plt.savefig(file_mw_mmt_kstest)
plt.clf()
plt.cla()

cdf_fb = cdf(mmt_fb[m31_fb_id])
cdf_cs = cdf(mmt_cs[m31_cs_id])


print('=== MMT M31  ===')
pos05 = np.where(cdf_cs[0] == 6.75)
print('CS, half: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 6.75)
print('FB, half: ', 1.0 - cdf_fb[1][pos05[0][0]])
pos05 = np.where(cdf_cs[0] == 10)
print('CS, quart: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 10)
print('FB, quart: ', 1.0 - cdf_fb[1][pos05[0][0]])
#print(cdf_cs[1])

#print(stats.ks.test(cdf_fb[0], cff_cs[0]))

kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)
print('KST confidence level = ', kst[0] * 16.48)
print('KST confidence level = ', kst[1] * ksnorm)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)

axs.set_xlabel('t [Gyr]', size = size_lab)
#axs.set_ylabel('n(<T)$')
#axs.set_ylabel('$n(< \\tau _M)$')
#plt.title('M31 MMT Cumulative')
#plt.xticks([])
plt.yticks([])
plt.tight_layout()
plt.savefig(file_m31_mmt_kstest)
plt.clf()
plt.cla()

# FORMATION TIME

cdf_fb = cdf(ft_fb[mw_fb_id])
cdf_cs = cdf(ft_cs[mw_cs_id])

print('=== FT MW ===')
pos05 = np.where(cdf_cs[0] == 6.75)
print('CS, half: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 6.75)
print('FB, half: ', 1.0 - cdf_fb[1][pos05[0][0]])
pos05 = np.where(cdf_cs[0] == 10)
print('CS, quart: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 10)
print('FB, quart: ', 1.0 - cdf_fb[1][pos05[0][0]])

pos05 = np.where(cdf_fb[0] == 9)
print('FB 9: ', 1.0 - cdf_fb[1][pos05[0][0]])

kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)
print('KST confidence level = ', kst[0] * 16.48)
print('KST confidence level = ', kst[1] * ksnorm)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)

axs.set_xlabel('t [Gyr]', size = size_lab)
axs.set_ylabel('$n(< \\tau _F)$', size = size_lab)
#axs.set_ylabel('n(<T)$')
#plt.title('MW FT Cumulative')
plt.yticks(ticks)
plt.tight_layout()
plt.savefig(file_mw_ft_kstest)
plt.clf()
plt.cla()

cdf_fb = cdf(ft_fb[m31_fb_id])
cdf_cs = cdf(ft_cs[m31_cs_id])

print('=== FT M31 ===')
pos05 = np.where(cdf_cs[0] == 6.75)
print('CS, half: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 6.75)
print('FB, half: ', 1.0 - cdf_fb[1][pos05[0][0]])
pos05 = np.where(cdf_cs[0] == 10)
print('CS, quart: ', 1.0 - cdf_cs[1][pos05[0][0]])
pos05 = np.where(cdf_fb[0] == 10)
print('FB, quart: ', 1.0 - cdf_fb[1][pos05[0][0]])

pos05 = np.where(cdf_fb[0] == 9)
print('FB 9: ', 1.0 - cdf_fb[1][pos05[0][0]])



kst = stats.ks_2samp(cdf_fb[0], cdf_cs[0])
print(kst)
print('KST confidence level = ', kst[0] * 16.48)
print('KST confidence level = ', kst[1] * ksnorm)

(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

axs.axis([x_min, x_max, y_min, y_max])
axs.plot(cdf_fb[0], cdf_fb[1], color=col_fb)
axs.plot(cdf_cs[0], cdf_cs[1], color=col_cs)
#axs.set_ylabel('T_{F}(<T)$')
#axs.set_ylabel('$n(< \\tau _F)$')
#plt.title('M31 FT Cumulative')
axs.set_xlabel('t [Gyr]', size = size_lab)
plt.yticks([])
plt.tight_layout()
plt.savefig(file_m31_ft_kstest)
plt.clf()
plt.cla()
