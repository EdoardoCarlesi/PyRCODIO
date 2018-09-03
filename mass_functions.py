#!/usr/bin/python

import matplotlib.pyplot as plt
import scipy.stats as sp
import numpy as np
import os

from libio.read_ascii import *
from config import *
from libcosmo.utils import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.lg_plot import *
import pickle

simuruns = simu_runs()
n_subrun = 10

for simu in simuruns:

	for i_sub in (0, n_subrun):
		subrun = '%02d' % i_sub
		m_fname='saved/'+simu+'_'+subrun+'_mains_all.pkl'
		s_fname='saved/'+simu+'_'+subrun+'_sats.pkl'

		try:
			hand_main = open(m_fname, 'r')
			hand_sats = open(s_fname, 'r')
			main = pickle.load(hand_main)
			sats = pickle.load(hand_sats)

			n_main = len(main)
			r_sub = 1.5
			n_lg = 2
			n_simu += 1

			print 'Found: ', s_fname

		except:	
			n_lg = 0

		# The first two are the main LG members
		#for i_lg in range(0, n_lg):


