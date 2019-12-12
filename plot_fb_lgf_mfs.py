from libcosmo.utils import *
from libcosmo.units import *
import time
import pickle
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import os

fOutMFs = 'saved/lgs_all_mfs.pkl'
fBaseLG = 'saved/lgs_center_5000.0_'
fBaseFB = 'saved/rand_mfs_'

randStr = ['00', '01', '02', '03', '04']
iIni = 0; iEnd = 100
gIni = 0; gEnd = 30

subPaths = gen_subpaths(iIni, iEnd, gIni, gEnd)

allMFs = []
#doLGs = True
doLGs = False

nBins = 30

mf_bins = np.linspace(9.1, 13.5, nBins)
x_bins = np.zeros((nBins), dtype='float')
binX = np.zeros((nBins-1), dtype='float')
n_bins = np.zeros((nBins), dtype='int')

for ib in range(0, nBins):
    x_bins[ib] = np.power(10.0, mf_bins[ib])

for ib in range(0, nBins-1):
    binX[ib] = 0.5 * (x_bins[ib] + x_bins[ib+1])

if doLGs: 
    for sP in subPaths:
        thisPkl = fBaseLG + sP + '.pkl'

        if os.path.isfile(thisPkl):
            print('Found: ', thisPkl)

            thisFile = open(thisPkl, 'rb')
            thisHalos = pickle.load(thisFile)

            ms = []
            for h in thisHalos:
                ms.append(h.m)

            thisMF = mass_function(ms)
            allMFs.append(thisMF)

    print('Saving to file: ', fOutMFs)
    fMFs = open(fOutMFs, 'wb')
    pickle.dump(allMFs, fMFs)
    fMFs.close()
else:

    print('Loading LGs (LGF) from file: ', fOutMFs)
    fMFs = open(fOutMFs, 'rb')
    allMFs = pickle.load(fMFs)
    fMFs.close()

    binMFs = [None] * nBins
    binMFsFB = [None] * nBins
    errMFs = [None] * nBins
    errMFsFB = [None] * nBins

    for i in range(0, nBins):
       binMFs[i] = []
       binMFsFB[i] = []

    for mf in allMFs:
        bx = mf[0]
        by = mf[1]
        this_bins = bin_this(x_bins, bx, by)
#        print(this_bins)
        for i in range(0, nBins-1):
            value = this_bins[1][i]

            if value > 0.0:
                binMFs[i].append(value)
    
    for run in range(0, 5):
        runStr = '%02d' % run
        thisFileFB = fBaseFB + runStr + '.pkl' 

        
        if os.path.isfile(thisFileFB):
            print('Loading LGs (STD) from file: ', thisFileFB)
            thisFB = open(thisFileFB, 'rb')
            thisMFs = pickle.load(thisFB)
    #        print(len(thisMFs))

            for mf in thisMFs:
                bx = mf[0]
                by = mf[1]
                this_bins = bin_this(x_bins, bx, by)

                for i in range(0, nBins-1):
                    value = this_bins[1][i]
     #               print(value)

                    if value > 0.0:
                        binMFsFB[i].append(value)

    perc0 = 20; perc1 = 50; perc2 = 100 - perc0

    minBinsCS = np.zeros((nBins-1));       medBinsCS = np.zeros((nBins-1));     maxBinsCS = np.zeros((nBins-1))
    minBinsFB = np.zeros((nBins-1));       medBinsFB = np.zeros((nBins-1));     maxBinsFB = np.zeros((nBins-1))
#    print(binMFsFB[0])           

    for i in range(0, nBins-1):
        if len(binMFs[i]) > 0:
            minBinCS = np.percentile(binMFs[i], perc0)
            medBinCS = np.percentile(binMFs[i], perc1)
            maxBinCS = np.percentile(binMFs[i], perc2)

            minBinsCS[i] = minBinCS
            medBinsCS[i] = medBinCS
            maxBinsCS[i] = maxBinCS

        if len(binMFsFB[i]) > 0:
            minBinFB = np.percentile(binMFsFB[i], perc0)
            medBinFB = np.percentile(binMFsFB[i], perc1)
            maxBinFB = np.percentile(binMFsFB[i], perc2)

            minBinsFB[i] = minBinFB
            medBinsFB[i] = medBinFB
            maxBinsFB[i] = maxBinFB

    ratioMed = np.zeros((nBins-1))
    ratioMin = np.zeros((nBins-1))
    ratioMax = np.zeros((nBins-1))
    ratioY = np.zeros((nBins-1))
    
    for i in range(0, nBins-1):
        if medBinsFB[i] > 0 and medBinsCS[i] > 0:
            rMed = medBinsCS[i] / medBinsFB[i]
            rMin = minBinsCS[i] / minBinsFB[i] - rMed * 0.5
            rMax = maxBinsCS[i] / maxBinsFB[i] + rMed * 0.5
        else:
            rMed = 1.0

        ratioMed[i] = rMed
        ratioMax[i] = rMax
        ratioMin[i] = rMin
        ratioY[i] = 1.0

    # Plot stuff
    col = 'black'
    col_cs = 'red'
    col_fb = 'blue'

    plt.figure(0)
    ax0 = plt.subplot2grid((8, 8), (0, 0), rowspan=6, colspan=8)
    ax1 = plt.subplot2grid((8, 8), (6, 0), rowspan=2, colspan=8)
    ax0.set_xscale('log')
    ax1.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xticks([])
    ax0.plot(binX, medBinsCS, color=col_cs)
    ax0.fill_between(binX, minBinsCS, maxBinsCS, facecolor=col_cs, alpha=0.4)
    ax0.plot(binX, medBinsFB, color=col_fb)
    ax0.fill_between(binX, minBinsFB, maxBinsFB, facecolor=col_fb, alpha=0.2)
    ax1.plot(binX, ratioY, color=col)
    ax1.plot(binX, ratioMed, color=col_fb)
    ax1.fill_between(binX, ratioMin, ratioMax, facecolor=col_fb, alpha=0.2)

    file_out='mass_functions_cs_vs_fb.png'

    plt.tight_layout()
    plt.savefig(file_out)
    plt.clf()
    plt.cla()


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
