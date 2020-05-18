import numpy as np

def mah(a, b, z):
    fac = - a * (np.sqrt(1.0 + z) - 1.0)
    mz_m0 = np.power(1.0 + z, b) * np.exp(fac)
    
    return fac * mz_m0


z0 = 0.0
z1 = 6.0
nz = 10000

dz = (z1 - z0) / nz

#a0 = 5.54; b0 = 2.71
a0 = 6.33; b0 = 2.91
#a0 = 6.24; b0 = 2.71
#a0 = 6.01; b0 = 2.66

a1 = 5.54; b1 = 2.71
#a1 = 4.90; b1 = 2.23
#a1 = 4.50; b1 = 2.23

diff = 0.0

for iz in range(0, nz):
    z = z0 + iz * dz
    v0 = mah(a0, b0, z)
    v1 = mah(a1, b1, z)

    diff = diff + np.sqrt(np.power(v0 - v1, 2.0))

print('RMS: ', diff / nz, dz)
