import matplotlib.pyplot as plt
import numpy as np
import os

from libio.read_ascii import *
from config import *
from libcosmo.utils import *
from libcosmo.halos import *
from libcosmo.find_halo import *
from libcosmo.lg_plot import *
import pickle


fn_all_subs = 'saved/all_subs__2048.pkl'
fn_all_host = 'saved/all_host__2048.pkl'

f_subs = open(fn_all_subs, 'rb')
subs = pickle.load(f_subs)

f_host = open(fn_all_host, 'rb')
host = pickle.load(f_host)

m0 = 1.e+12

lg_name = ['M31', 'MW']

plt.rc({'text.usetex': True})
plt.axis([0.5, 4.7, 1, 150])
plt.xlabel(r'$ \frac{M} {10 ^ {12} M_{\odot}}$')
plt.ylabel(r'$ N (> 5 \times 10^8 M_{\odot})$')
plt.scatter(host[:, :, :]/m0, subs, c='black', s=20)
#axes = plt.gca()
#axes.set_xlim([0.1, 500000000000000000000000000000.0])
#axes.set_ylim([10, 150])
plt.title('Full sample')
fname = 'subs_host_all.png'

plt.tight_layout()
plt.savefig(fname)
plt.clf()
plt.cla()

for ilg in range(0, 2):
    n_err = np.zeros((2, 10), dtype='float')
    m_err = np.zeros((2, 10), dtype='float')
    n_med = np.zeros(10, dtype='float')
    m_med = np.zeros(10, dtype='float')

    for ihs in range(0, 10):
        p0 = 10; p1 = 50; p2 = 90
        min_n = np.percentile(subs[ihs, ilg, :], p0)
        med_n = np.percentile(subs[ihs, ilg, :], p1)
        max_n = np.percentile(subs[ihs, ilg, :], p2)

        n_sig_min = abs(min_n-med_n)/med_n
        n_sig_max = abs(max_n-med_n)/med_n
        n_err[0, ihs] = abs(med_n - min_n) 
        n_err[1, ihs] = abs(max_n - med_n) 
        n_med[ihs] = med_n

        p0 = 20; p1 = 50; p2 = 80
        min_m = np.percentile(host[ihs, ilg, :], p0)
        med_m = np.percentile(host[ihs, ilg, :], p1)
        max_m = np.percentile(host[ihs, ilg, :], p2)
        m_sig_min = (abs(min_m-med_m)/med_m)
        m_sig_max = (abs(max_m-med_m)/med_m)
        m_err[0, ihs] = (abs(np.log10(med_m) - np.log10(min_m)))
        m_err[1, ihs] = (abs(np.log10(max_m) - np.log10(med_m)))
        m_med[ihs] = np.log10(med_m)

    plt.rc({'text.usetex': True})
    plt.xlabel(r'$\log_{10} \frac{M} {M_{\odot}}$')
    plt.ylabel(r'$ N (> 5 \times 10^8 M_{\odot})$')
    plt.errorbar(m_med, n_med, xerr = m_err , yerr=n_err, fmt='o', color='black')
    plt.title(lg_name[ilg])
    fname = 'subs_' + lg_name[ilg] + '.png'

    plt.tight_layout()
    plt.savefig(fname)
    plt.clf()
    plt.cla()

