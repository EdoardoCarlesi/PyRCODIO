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

fBaseLG = 'saved/'
fBaseFB = 'saved/'

iIni = 0; iEnd = 10
gIni = 0; gEnd = 2

subPaths = gen_subpaths(iIni, iEnd, gIni, gEnd)

print(subPaths)


'''
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
