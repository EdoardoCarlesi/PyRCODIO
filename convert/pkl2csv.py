'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    pkl2csv.py: convert pickle files to csv
'''

import pickle as pkl
import pandas as pd
import config as cfg


def convert_fb_pkl(vweb=False):
    runs = cfg.gen_runs(0, 5)

    for run in runs:
        this_pkl = 'output/lg_fullbox_' + run + '.pkl'
        f_pkl = open(this_pkl, 'rb')
        all_lgs = pkl.load(f_pkl)
        f_pkl.close()

        cols = all_lgs[0].header(dump=False)
        lg_df = pd.DataFrame(columns = cols)
        print(run, ') N LG Pairs: ', len(all_lgs))

        for i, lg in enumerate(all_lgs):
            this_row = lg.info(dump=False)
            this_series = pd.Series(this_row, index = lg_df.columns)

            lg_df = lg_df.append(this_series, ignore_index=True)

        lg_df.drop_duplicates(inplace=True)
        print('N unique elements: ', len(lg_df['Xc_LG'].unique()))

        print(lg_df.head())
        this_csv = 'output/lg_fullbox_' + run + '.csv'
        lg_df.to_csv(this_csv)



print('Converting FullBox CSV to .pkl')
convert_fb_pkl()
print('Done.')
