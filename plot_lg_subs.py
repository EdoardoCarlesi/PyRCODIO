import matplotlib.pyplot as plt
import scipy.stats as sp
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

print(host)
print(subs)

plt.scatter(host, subs, s=0.5, color='black')
plt.show()

plt.savefig('sats.png')
