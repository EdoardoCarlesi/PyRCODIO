'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_lg_extract.py: Analyze a previously produced csv file containing LG properties
'''

import matplotlib.pyplot as plt
import read_files as rf
import halo_utils as hu
import seaborn as sns
import numpy as np
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn import metrics


# Output file containing all the LG properties 
out_all_lgs_csv = 'output/lg_pairs_2048.csv'

data_lg = pd.read_csv(out_all_lgs_csv)
data_lg.head()

print(data_lg.columns)

data_lg['Mtot'] = data_lg['M_M31'] + data_lg['M_MW']
data_lg['Mratio'] = data_lg['M_M31'] / data_lg['M_MW']
data_lg['Nsubtot'] = data_lg['Nsub_MW'] + data_lg['Nsub_M31']

cols = ['R', 'Vtan', 'Vrad', 'Nsubtot']
#cols = ['R', 'Vtan', 'Vrad', 'Mratio', 'Nsubtot']

X_df = data_lg[cols]
y_df = data_lg['Mtot']

model = LinearRegression(normalize=True)
X_train, X_test, y_train, y_test = train_test_split(X_df, y_df, test_size = 0.3, random_state=69)
model.fit(X_train, y_train)

predictions = model.predict(X_df)

plt.scatter(y_df, predictions)
#sns.distplot(y_df-predictions, bins=10)

mae = metrics.mean_absolute_error(y_df, predictions)
msq = metrics.mean_squared_error(y_df, predictions)
rme = np.sqrt(msq)

print('RME: ', mae/1.e+12, ' MAE: ', mae/1.e+12)

plt.show()
