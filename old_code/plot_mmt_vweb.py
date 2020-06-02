import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd

path_dict_lgs = 'saved/dict_lgs_code.pkl'
path_dict_mmt = 'saved/dict_lgs_mmt.pkl'
path_web = 'saved/lgf_web_'

f_dict_lgs = open(path_dict_lgs, 'rb')
f_dict_mmt = open(path_dict_mmt, 'rb')

all_lgs = pickle.load(f_dict_lgs)
all_mmt = pickle.load(f_dict_mmt)

l1 = []; l2 = []; l3 = []; mmt = [];

for lg_code in all_lgs:
    ids = []
    xyz = []
    
    this_web = path_web + lg_code + '.pkl'
    f_web = open(this_web, 'rb')
    v_web = pickle.load(f_web, encoding="latin1")

    this_data = all_lgs[lg_code]
    xyz.append(int(this_data[0]))
    xyz.append(int(this_data[1]))
    xyz.append(int(this_data[2]))
    ids.append(str(this_data[3]))
    ids.append(str(this_data[4]))

    try:
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]

        this_mmt = all_mmt[ids[1]]
        mmt.append(float(this_mmt))
        evs = v_web[:, x, y, z]

        if abs(evs[0]) < 0.01:
            fact = 256.0
        else:
            fact = 1.0

        l1.append(float(evs[0]) * fact)
        l2.append(float(evs[1]) * fact)
        l3.append(float(evs[2]) * fact)

    except:
        'Do not print anything'
        #print(ids[0], ' not found.')

    f_web.close()

n_bins = 25
l1_bins = np.linspace(0.0, 1.0, n_bins)
mmt_bins = np.zeros((n_bins), dtype='float') 
err_mmt_bins = np.zeros((n_bins), dtype='float') 
n_mmt_bins = np.zeros((n_bins), dtype='int')

n_tot = len(l1)

step = 1.0 / n_bins

for il in range(0, n_tot):

    this_bin = int((l1[il] + l2[il])/ step)

    if this_bin < n_bins and this_bin > 0:
       mmt_bins[this_bin] = mmt_bins[this_bin] + mmt[il]
       n_mmt_bins[this_bin] = n_mmt_bins[this_bin] + 1


for im in range(0, n_bins):
    this_n = n_mmt_bins[im]

    if this_n > 1:
        mmt_bins[im] = mmt_bins[im] / float(this_n)
        err_mmt_bins[im] = mmt_bins[im] / np.sqrt(this_n)

plt.axis([0.0, 0.75, 1.0, 13.0])
plt.scatter(l1_bins, mmt_bins)
plt.savefig('test_evs.png')

'''
#    print(ids, xyz)

# Put all the evs and mmts in a dataframe
df = { 'Ev1' : l1,
        'Ev2' : l2,
        'Ev3' : l3, 
        'MMT' : mmt }

df = pd.DataFrame(df, columns=['Ev1', 'Ev2', 'Ev3', 'MMT'])
#bins = [0.0, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4]
#np.searchsorted(bins, df['percentage'].values)
#df['binned'] = pd.cut(df['Ev1'], bins)
#df['binned'] = np.searchsorted(bins, df['MMT'].values)
#print(df)
#s = df.groupby(pd.cut(df['percentage'], bins=bins)).size()
#s = df.groupby(pd.cut(df['Ev1'], bins=bins)).size()
#df1['binned'] = pd.cut(df1['Score'], bins)
#print(s)
#print(df)

#print(mmt)
#print(l1)
#print(l2)
#print(l3)

plt.axis([0.01, 1.0, 0.0, 13.0])
plt.scatter(l1, mmt)
plt.savefig('test_evs.png')
'''



