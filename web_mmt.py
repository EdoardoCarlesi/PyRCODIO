from config import *
import random
import time
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from libcosmo.halos import *

# Save all the selected merger trees to these files
out_ev1 = 'saved/rand_ev1.pkl'
out_ev2 = 'saved/rand_ev2.pkl'
out_ev3 = 'saved/rand_ev3.pkl'

#evs_lg = 'saved/lgs_evs.pkl'
evs_lg = 'saved/lgs_evs_all.pkl'
in_lgs = 'saved/rand_web_lgs_'

base_lgs = 'saved/lgs_'
base_web = 'saved/lgf_web_'

out_dict_mmt = 'saved/dict_lgs_mmt.pkl'
out_dict_code = 'saved/dict_lgs_code.pkl'

fileMMT_m31='saved/all_trees_lgf_m31.pkl'
fileMMT_mw='saved/all_trees_lgf_mw.pkl'

f_mmt_mw = open(fileMMT_mw, 'rb')
f_mmt_m31 = open(fileMMT_m31, 'rb')

mmt_mw = pickle.load(f_mmt_mw)
mmt_m31 = pickle.load(f_mmt_m31)

id2mmt = dict()
code2pos = dict()

iIni = 0
iEnd = 80
gIni = 0
gEnd = 10

offset = 20
grid = 64
box = 100.0

kpc = 1.0e-3

cell_unit = float(grid / box)

for this_mw in mmt_mw:
    this_mmt = this_mw.last_major_merger(False) * 0.25
    this_id = this_mw.IDs[0]
    id_str = str(this_id)
    id2mmt[id_str] = this_mmt

for this_m31 in mmt_m31:
    this_mmt = this_m31.last_major_merger(False) * 0.25
    this_id = this_m31.IDs[0]
    id_str = str(this_id)
    id2mmt[id_str] = this_mmt

f_dict = open(out_dict_mmt, 'wb')
pickle.dump(id2mmt, f_dict)

#f_dict = open(out_dict, 'rb')
#id2mmt_bis = pickle.load(f_dict)

for i in range(iIni, iEnd):
    i_str = '%02d' % i

    for g in range(gIni, gEnd):
        g_str = '%02d' % g

        sub_dir = i_str + '_' + g_str

        f_web = base_web + sub_dir + '.pkl'
        f_lgs = base_lgs + sub_dir + '.pkl'

        existWEB = os.path.isfile(f_web)
        existLGs = os.path.isfile(f_lgs)

        if existWEB and existLGs:
            file_lgs = open(f_lgs, 'rb') 
            file_web = open(f_web, 'rb') 

            this_lgs = pickle.load(file_lgs)
            #this_web = pickle.load(file_web)
            this_com = this_lgs[0].geo_com()
            this_id_mw = this_lgs[0].LG1.ID
            this_id_m31 = this_lgs[0].LG2.ID

            i_x = int(this_com[0] * cell_unit * kpc) - offset
            i_y = int(this_com[1] * cell_unit * kpc) - offset
            i_z = int(this_com[2] * cell_unit * kpc) - offset
    
            coord = [i_x, i_y, i_z, this_id_mw, this_id_m31]
            code2pos[sub_dir] = coord
            
f_dict_code = open(out_dict_code, 'wb')
pickle.dump(code2pos, f_dict_code)

#print(code2pos)
            
#print(i_x, i_y, i_z)        
# print()
#print(id2mmt_bis)
#webname = 'saved/lgf_web_79_08.pkl'
#f_web = open(webname, 'rb')
#vweb = pickle.load(f_web)
#print(vweb[0, 12, 12, 12] * 512.0) 





