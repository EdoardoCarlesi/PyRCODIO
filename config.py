#!/usr/bin/python
from libcosmo.find_halos import *
from libcosmo.halo import *
import numpy as np
import os

def simu_runs():
	all_runs = [] 

	all_runs.append('00_06')	
	all_runs.append('01_12') 	
	all_runs.append('17_10') 
	all_runs.append('17_13') 
	all_runs.append('34_13')
	all_runs.append('45_17')
	all_runs.append('55_02') 
	all_runs.append('64_14')
	all_runs.append('09_18')
	all_runs.append('62_14')

	return all_runs

# 2048 runs
def lg_models():
	all_lg_models = []
	model_index = dict()
	model_count = 0

	this_model = '00_06'
	d_max = 5000.; r_iso = 1700.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 00_06 
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '01_12'
	d_max = 8000.; r_iso = 1550.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 8.5e+12; ratio_max = 3.5; vrad_max = -1.0	# 01_12
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '17_10'
	d_max = 8000.; r_iso = 1900.; r_max = 1700.; r_min = 200.; m_min = 9.0e+11; m_max = 10.e+12; ratio_max = 3.5; vrad_max = -1.0	# 17_10
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '17_13'
	d_max = 8000.; r_iso = 1900.; r_max = 1700.; r_min = 200.; m_min = 9.0e+11; m_max = 10.e+12; ratio_max = 3.5; vrad_max = -1.0	# 17_10
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '34_13'
	d_max = 7000.; r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 9.0e+11; m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 34_13
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '45_17'
	d_max = 5000.; r_iso = 1700.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 45_17
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({lg_model:model_count}) ; model_count += 1

	this_model = '55_02'
	d_max = 7000.; r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 9.0e+11; m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 55_02
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '64_14'
	d_max = 7000.; r_iso = 1800.; r_max = 1400.; r_min = 200.; m_min = 8.0e+11; m_max = 4.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 64_14
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '09_18'
	d_max = 6000.; r_iso = 1700.; r_max = 1200.; r_min = 200.; m_min = 0.9e+12; m_max = 7.0e+12; ratio_max = 4.5; vrad_max = -1.0	# 64_14
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	this_model = '62_14'
	d_max = 9000.; r_iso = 1000.; r_max = 1200.; r_min = 100.; m_min = 4.0e+11; m_max = 10.5e+12; ratio_max = 4.5; vrad_max = 100.0	# 64_14
	lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	lg_model.model_code = this_model
	all_lg_models.append(lg_model)
	model_index.update({this_model:model_count}) ; model_count += 1

	return (all_lg_models, model_index)

class Settings:
	home_dir = ''
	data_dir = ''
	outp_dir = ''	
	env_type = None
	data_dir = None 
	ahf_path = ''
	sub_path = ''
	base_path = ''
	plot_out = ''
	base_run = None
	run_num = None
	sub_path_ahf = ''
	sub_path_ics = ''
	sub_path_snap = ''
	resolution = None
	n_ic = None
	n_z0 = None
	file_z0 = ''
	base_out = '' 

	snapshot=''

	file_lg_name = ''
	file_lgall_name = ''
	file_lare_name = ''
	file_codes_name = ''

	file_ahf_in = ''
	file_ic_in = ''
	file_z0_in = ''
	file_ahf = ''
	file_ic = ''
	file_z0 = ''

	box_center = []

	def __init__(self, home_dir, outp_dir, env_type, resolution, snapshot):
		self.snapshot = snapshot
		self.home_dir = home_dir
		self.outp_dir = outp_dir
		self.env_type = env_type
		self.resolution = resolution

		self.base_run = ''
		self.info()
		self.init_paths()
		self.init_ascii()

	def info(self):
		print 'Home dir: ', self.home_dir
		print 'Data dir: ', self.data_dir
		print 'Lare out: ', self.file_lare_name
		print 'Env type: ', self.env_type

	def re_init(self):
		self.init_paths()
		self.init_ascii()
		#self.init_files()
	
	# Initialize local paths - where to look for all the files and stuff
	def init_paths(self):
		if self.env_type == "HESTIA" :
			self.data_dir = 'HESTIA/'
			self.base_path = self.home_dir + self.data_dir
			self.sub_path = 'DM_ONLY/' 
			self.sub_path_ahf = 'AHF_output/'
			self.sub_path_ics = 'ICs/'
			self.sub_path_snap = 'output/snapdir_127/'

		elif self.env_type == "512box100" :
			self.data_dir = '512box100/'
			self.base_path = self.home_dir + self.data_dir
			self.sub_path = ''
			self.sub_path_ahf = 'AHF/'
			self.sub_path_snap = 'SNAPS/'
			self.sub_path_ics = 'ICs/'

		elif self.env_type == "zoom" :
			self.data_dir = 'DATA/' 
			self.base_path = self.home_dir + self.data_dir
			self.sub_path = ''
			self.sub_path_ahf = 'AHF/'
			self.sub_path_snap = 'SNAPS/'
			self.sub_path_ics = 'ICs/'

		else:
			self.data_dir = 'DATA/' 
			self.base_path = self.home_dir + self.data_dir 
			self.sub_path = ''
			self.sub_path_ahf = ''
			self.sub_path_ics = ''
			self.sub_path_snap = ''

	# Global files storing all the output information in txt format
	def init_ascii(self):
		if self.env_type == "HESTIA" :
			self.file_lg_name = self.base_path + self.outp_dir + 'lg_candidates_HESTIA_'+self.resolution+'.txt'
			self.file_lare_name = self.base_path + self.outp_dir + 'lg_candidates_LaReHESTIA_'+self.resolution+'.txt'
			self.file_lgall_name = self.base_path + self.outp_dir + 'lg_candidates_HESTIAall_'+self.resolution+'.txt'
			self.file_codes_name = self.base_path + self.outp_dir + 'lg_codes_HESTIA_'+self.resolution+'.txt'

		elif self.env_type == "512box100" :
			self.file_lg_name = self.base_path + self.outp_dir + 'lg_candidates_512box100.txt'
			self.file_lare_name = self.base_path + self.outp_dir + 'lg_candidates_LaRe512box100.txt'
			self.file_lgall_name = self.base_path + self.outp_dir + 'lg_candidates_LaRe512box100all.txt'
			self.file_codes_name = self.base_path + self.outp_dir + 'lg_codes_512box100.txt'

		elif self.env_type == "zoom" :
			self.file_lg_name = self.base_path + self.outp_dir + 'lg_candidates_' + self.resolution + '_' + self.base_run + '.txt'
			self.file_lgsub_name = self.base_path + self.outp_dir + 'lg_candidates_' + self.resolution + '_' + self.base_run + '_subhalos.txt'
			self.file_lare_name = self.base_path + self.outp_dir + 'lg_candidates_LaRe_' + self.resolution + '.txt'
			self.file_lgall_name = self.base_path + self.outp_dir + 'lg_candidates_LaReAll_' + self.resolution + '.txt'
			self.file_codes_name = self.base_path + self.outp_dir + 'lg_codes_' + self.resolution + '.txt'
		else:
			self.file_lg_name = self.base_path + self.outp_dir + 'lg_candidates_LGF_'+self.resolution+'.txt'
			self.file_lgall_name = self.base_path + self.outp_dir + 'lg_candidates_LGFall_'+self.resolution+'.txt'
			self.file_lare_name = self.base_path + self.outp_dir + 'lg_candidates_LaReLGF_'+self.resolution+'.txt'
			self.file_codes_name = self.base_path + self.outp_dir + 'lg_codes_'+self.resolution+'.txt'


	def init_files(self, base_run, run_num):
		self.base_run = base_run
		self.run_num = run_num
	
		if self.resolution == "1024" or self.resolution == "2048":
			self.base_run = self.base_run + '/' + self.run_num

		if self.env_type == "HESTIA" :
			self.snapshot = 'HESTIA_100Mpc_512_'+self.base_run+'.127.z0.000.AHF_halos'
			self.file_z0 = 'snapshot_127'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.home_dir + '/HESTIA/LaRe/' + self.base_out
			self.n_ic = 2
			self.n_z0 = 8
			self.file_ahf_in = self.base_path + '/'+ self.sub_path + self.base_run + '/' + self.sub_path_ahf + self.snapshot
			self.file_z0_in = self.base_path + '/' + self.sub_path + self.base_run + '/' + self.sub_path_snap + self.file_z0
			self.file_ic_in = self.base_path + '/' + self.sub_path + self.base_run + '/' + self.sub_path_ics + self.file_ic

		elif self.env_type == "512box100" :
			self.file_z0 = 'snapshot_054'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.base_path + self.resolution + '/' + self.base_run + '/' + self.base_out
			self.n_ic = 1
			self.n_z0 = 1
			self.file_ahf_in = self.base_path + self.sub_path_ahf  + self.base_run + '/' + self.snapshot
			self.file_z0_in = self.base_path +  self.sub_path_snap + self.base_run + '/' + self.file_z0
			self.file_ic_in = self.base_path +  self.sub_path_ics  + self.base_run + '/' + self.file_ic

		elif self.env_type == "zoom" :
			self.file_z0 = 'snapshot_054'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.base_path + self.resolution + '/' + self.base_run + '/' + self.base_out
			self.n_ic = 1
			self.n_z0 = 1
			self.ahf_path = self.base_path + self.resolution + '/' + self.base_run + '/'
			self.file_ahf_in = self.base_path + self.sub_path_ahf  + self.base_run + '/' + self.snapshot
			self.file_z0_in = self.base_path +  self.sub_path_snap + self.base_run + '/' + self.file_z0
			self.file_ic_in = self.base_path +  self.sub_path_ics  + self.base_run + '/' + self.file_ic
	
		else:
			self.file_z0 = 'snapshot_054'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.base_path + self.resolution + '/' + self.base_run + '/' + self.base_out
			self.n_ic = 2
			self.n_z0 = 1
			self.file_ahf_in = self.base_path + self.resolution + '/' + self.base_run + '/' + self.snapshot
			self.file_z0_in = self.base_path + self.resolution + '/' + self.base_run + '/' + self.file_z0
			self.file_ic_in = self.base_path + self.resolution + '/ICs/' + self.file_ic

		# This is the same for all
		self.plot_out =	self.base_path + self.outp_dir + self.base_run + '_particles_LG_LV.png'

	def get_zoom_output(self):
		return (self.file_lg_name, self.file_lgsub_name)

	def get_files_input(self):
		return (self.file_z0, self.base_out, self.n_ic, self.n_z0)












