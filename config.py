'''
    Python Routines for COsmology and Data I/O (PyRCODIO) v0.2
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    config.py = store local configuration variables
'''

import halo_utils as hu
import numpy as np
import os


def gen_runs(ini, end):
    ''' Very easy '''

    runs = []
    for i in range(ini, end):
        run = '%02d' % i
        runs.append(run)

    return runs


def gen_all_runs(i0, i1, g0, g1):
    runs = []

    main = gen_runs(i0, i1)
    subs = gen_runs(g0, g1)

    for m in main:
        for s in subs:
            r = m + '_' + s
            runs.append(r)

    return runs


def sub_runs():
    sub_runs = []

    sub_runs.append('00')
    sub_runs.append('01')
    sub_runs.append('02')
    sub_runs.append('03')
    sub_runs.append('04')
    sub_runs.append('05')
    sub_runs.append('06')
    sub_runs.append('07')
    sub_runs.append('08')
    sub_runs.append('09')

    return sub_runs


def simu_runs():
    all_runs = [] 

    all_runs.append('00_06')	
    all_runs.append('01_12') 	
    all_runs.append('17_10') 
    all_runs.append('17_13') 
    all_runs.append('34_13')
    all_runs.append('45_17')
    all_runs.append('55_02') 
    all_runs.append('09_18')
    all_runs.append('64_14')
    all_runs.append('37_11')
    all_runs.append('62_14')

    return all_runs


def lg_models():
        all_lg_models = []
        model_index = dict()
        model_count = 0

        this_model = '00_06'
        r_iso = 1700.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -1.0 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '01_12'
        r_iso = 1500.; r_max = 1500.; r_min = 200.; m_min = 5.0e+11; m_max = 8.5e+12; ratio_max = 3.5; vrad_max = 100.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '17_10'
        r_iso = 1700.; r_max = 1600.; r_min = 200.; m_min = 9.0e+11; m_max = 10.e+12; ratio_max = 3.5; vrad_max = 100.0	
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '17_13'
        r_iso = 1000.; r_max = 1900.; r_min = 150.; m_min = 5.0e+11; m_max = 10.e+12; ratio_max = 5; vrad_max = 500.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '34_13'; r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 9.0e+11; 
        m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0	
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '45_17'; 
        r_iso = 1700.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -1.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '55_02'
        r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 9.0e+11; m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '64_14'
        r_iso = 1300.; r_max = 1900.; r_min = 200.; m_min = 8.0e+11; m_max = 4.5e+12; ratio_max = 2.5; vrad_max = 100.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '09_18'
        r_iso = 1700.; r_max = 1900.; r_min = 200.; m_min = 0.7e+12; m_max = 7.0e+12; ratio_max = 4.5; vrad_max = 100.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '62_14'
        r_iso = 1500.; r_max = 1500.; r_min = 300.; m_min = 2.5e+12; m_max = 7.0e+12; ratio_max = 4.5; vrad_max = 500.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '37_11'
        r_iso = 1000.; r_max = 1500.; r_min = 100.; m_min = 5.0e+11; m_max = 5.0e+12; ratio_max = 2.5; vrad_max = 100.0
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '08_11'
        r_iso = 1500.; r_max = 2000.; r_min = 200.; m_min = 5.0e+11; m_max = 3.0e+12; ratio_max = 2; vrad_max = 100.0	
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = 'GENERIC'
        r_iso = 2000.; r_max = 1500.; r_min = 250.; m_min = 4.0e+11; m_max = 4.0e+12; ratio_max = 5.0; vrad_max = 100.0	# Exremely generic LG Model 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1
    
        ##################################################################
        # Different LG models used for the "How common is the LG?" paper #
        ##################################################################

        this_model = 'M1'
        r_iso = 2000.; r_max = 1500.; r_min = 250.; m_min = 4.0e+11; m_max = 5.0e+12; ratio_max = 4.0; vrad_max = 50.0	# Generic LG Model 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = 'M2'
        r_iso = 2000.; r_max = 1300.; r_min = 300.; m_min = 4.5e+11; m_max = 4.0e+12; ratio_max = 3.0; vrad_max = 0.0	# Generic LG Model 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = 'M3'
        r_iso = 2000.; r_max = 1000.; r_min = 350.; m_min = 5.0e+11; m_max = 3.0e+12; ratio_max = 3.0; vrad_max = -20.0	# Generic LG Model 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = 'M4'
        r_iso = 2000.; r_max = 900.; r_min = 400.; m_min = 5.5e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -40.0	# Generic LG Model 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = 'M5'
        r_iso = 2000.; r_max = 800.; r_min = 450.; m_min = 6.0e+11; m_max = 2.0e+12; ratio_max = 2.5; vrad_max = -60.0	# Generic LG Model 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = 'M6'
        r_iso = 2000.; r_max = 700.; r_min = 450.; m_min = 6.5e+11; m_max = 1.75e+12; ratio_max = 2.0; vrad_max = -80.0    # Generic LG Model 
        lg_model = hu.LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        return (all_lg_models, model_index)

