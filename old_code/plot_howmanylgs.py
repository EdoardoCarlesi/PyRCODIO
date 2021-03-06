from libcosmo.utils import *
from libcosmo.units import *
from matplotlib import gridspec
import time
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Open files
fname_lgs = 'LGs_models_density.txt'; file_lg_dens = 'LGs_model_dens.png'; doHalo = False;
#fname_lgs = 'HALOs_models_density.txt'; file_lg_dens = 'halos_model_dens.png'; doHalo = True;
f_lgs = open(fname_lgs, 'r')
lg_dens = f_lgs.readlines()

#print(lg_dens)

x = []
ratio1 = []
ratio2 = []
dens_fb1 = []
dens_fb2 = []
dens_cs = []

err1 = []; err2 = []; err3 = []

for rep in lg_dens:
    rep = rep.replace("\t", " ")
    rep = rep.replace("\n", "")
    rep = rep.split(" ")

    line = []
    for ele in rep:
        if ele != '':
            line.append(ele)

    print(line)
    dens_fb1.append(float(line[0]))
    dens_cs.append(float(line[1]))
    dens_fb2.append(float(line[2]))
    ratio1.append(float(line[1])/float(line[0]))
    ratio2.append(float(line[1])/float(line[2]))

    if doHalo == False:
        error1 = float(line[0]) / np.sqrt(float(line[3]));         err1.append(error1)   
        error2 = float(line[1]) / np.sqrt(float(line[5]));         err2.append(error2)   
        error3 = float(line[2]) / np.sqrt(float(line[4]));         err3.append(error3)   

print(ratio1)
print(dens_fb1)

'''
        GENERAL PLOT PROPERTIES
'''
plt.rc({'text.usetex': True})

col_cs = 'black'
col_fb1 = 'red'
col_fb2 = 'blue'
col_ratio1 = 'red'
col_ratio2 = 'blue'

# Filenames
file_vmax_mah = 'mah_vmax.png'

xticks = [1,2,3,4,5,6, 7]
labels = []
labels.append('M1')
labels.append('M2')
labels.append('M3')
labels.append('M4')
labels.append('M5')
labels.append('M6')
labels.append('')

plt.figure(0)
axs0 = plt.subplot2grid((8, 6), (0, 0), rowspan=6, colspan=6)
axs1 = plt.subplot2grid((8, 6), (6, 0), rowspan=2, colspan=6)

#plt.gca().legend('RAND', 'CS')

#x = [1, 2, 3, 4, 5, 6]
x = [1, 2, 3, 4, 5, 6]

if doHalo == True:
    axs0.plot(x, dens_fb1, color=col_fb1, label='RAND') 
    axs0.plot(x, dens_fb2, color=col_fb2, label='RAND$_ \delta$') 
    axs0.plot(x, dens_cs, color=col_cs, label='CS') 
else:
    axs0.errorbar(x, dens_fb1, yerr=err1, color=col_fb1, label='RAND') 
    axs0.errorbar(x, dens_fb2, yerr=err3, color=col_fb2, label='RAND$_ \delta$') 
    axs0.errorbar(x, dens_cs, yerr=err2, color=col_cs, label='CS') 

axs0.legend(loc = 'upper right', bbox_to_anchor=(0.95, 0.95)) #, shadow=True, ncol=2)

# Axes and title
axs0.set_ylabel('$ N_{h} Mpc ^{-3} h ^ 3$')
#axs0.set_ylabel('$ N_{LG} Mpc ^{-3} h ^ 3$')
axs0.set_yscale('log')

axs1.set_xlabel('Model')
axs1.set_ylabel('ratio')

axs0.set_xticks([])
#axs1.set_xticks(xticks)
#plt.xticks(xticks, labels)

if doHalo == True:
    axs1.set_yticks([0.5, 1.00, 1.5, 2.0])
    axs0.axis([0, 7, 1.5e-3, 1.8e-2])
    axs1.axis([0, 7, 0.4, 2.1])
else:
    axs1.set_yticks([2, 3, 4, 5])
    axs0.axis([0, 7, 9e-7, 2e-3])
    axs1.axis([0, 7, 2, 5.5])

axs1.plot(x, ratio1, color=col_ratio1) 
axs1.plot(x, ratio2, color=col_ratio2) 

plt.tight_layout()
plt.subplots_adjust(hspace=0.0, wspace=0.2) 
plt.savefig(file_lg_dens)
plt.clf()
plt.cla()
