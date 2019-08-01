from config import *
import random
import time
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Save all the selected merger trees to these files
out_ev1 = 'saved/rand_ev1.pkl'
out_ev2 = 'saved/rand_ev2.pkl'
out_ev3 = 'saved/rand_ev3.pkl'
#evs_lg = 'saved/lgs_evs.pkl'
evs_lg = 'saved/lgs_evs_all.pkl'
in_lgs = 'saved/rand_web_lgs_'

fileMMT_m31='saved/all_trees_lgf_m31.pkl'
f_mmt_m31 = open(fileMMT_m31, 'rb')

mmt_m31 = pickle.load(f_mmt_m31)
print(len(mmt_m31))
print(mmt_m31[4].last_major_merger(False) * (0.25))
print(mmt_m31[4].nPart[0])

f_ev_lg = open(evs_lg, 'rb')
lg_evs = pickle.load(f_ev_lg)
print(len(lg_evs[1]))

id2mmt = dict()



