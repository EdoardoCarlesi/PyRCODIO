import numpy as np
import matplotlib.pylab as plt
from libio.read_ascii import *
from libcosmo.halos import *
from scipy import stats

ahf_path='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'


halos = read_ahf(ahf_path)

v = []
m = []

for hl in halos:
    v.append(np.log10(hl.vmax))
    m.append(np.log10(hl.m))

slope, intercept, r_value, p_value, std_err = stats.linregress(m, v)

print(slope, intercept)

nstep = 100
x0 = 8.0
x1 = 13.0
dx = (x1 - x0) / nstep 

mx = []; vy = []
for im in range(0, nstep):
    x = x0 + dx * im
    y = intercept + slope * x
    vy.append(y)
    mx.append(x)

plt.plot(m, v)
plt.plot(mx, vy)
#plt.show()
plt.savefig('ProjectsPlots/zzz_vmax.png')

