from libcosmo.halos import *
from libcosmo.units import *
from libcosmo.utils import *
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
	all_runs.append('09_18')
	all_runs.append('64_14')
	all_runs.append('37_11')
	all_runs.append('62_14')

	return all_runs

# 2048 runs
def lg_models():
        all_lg_models = []
        model_index = dict()
        model_count = 0

        this_model = '00_06'
        r_iso = 1700.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 00_06 
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '01_12'
        r_iso = 1500.; r_max = 1500.; r_min = 200.; m_min = 5.0e+11; m_max = 8.5e+12; ratio_max = 3.5; vrad_max = 100.0	# 01_12
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '17_10'
        r_iso = 1700.; r_max = 1600.; r_min = 200.; m_min = 9.0e+11; m_max = 10.e+12; ratio_max = 3.5; vrad_max = 100.0	# 17_10
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '17_13'
        r_iso = 1200.; r_max = 1900.; r_min = 200.; m_min = 5.0e+11; m_max = 10.e+12; ratio_max = 5; vrad_max = 500.0	# 17_13
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '34_13'; r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 9.0e+11; 
        m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 34_13
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '45_17'; 
        r_iso = 1700.; r_max = 1200.; r_min = 200.; m_min = 5.0e+11; m_max = 2.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 45_17
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '55_02'
        r_iso = 1900.; r_max = 1500.; r_min = 200.; m_min = 9.0e+11; m_max = 8.5e+12; ratio_max = 2.5; vrad_max = -1.0	# 55_02
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '64_14'
        r_iso = 1300.; r_max = 1900.; r_min = 200.; m_min = 8.0e+11; m_max = 4.5e+12; ratio_max = 2.5; vrad_max = 100.0	# 64_14
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '09_18'
        r_iso = 1700.; r_max = 1900.; r_min = 200.; m_min = 0.7e+12; m_max = 7.0e+12; ratio_max = 4.5; vrad_max = 100.0	# 64_14
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '62_14'
        r_iso = 1000.; r_max = 1200.; r_min = 100.; m_min = 4.0e+11; m_max = 10.5e+12; ratio_max = 4.5; vrad_max = 100.0	# 64_14
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        this_model = '37_11'
        r_iso = 1000.; r_max = 1500.; r_min = 100.; m_min = 5.0e+11; m_max = 5.0e+12; ratio_max = 2.5; vrad_max = 100.0	# 34_13
        lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
        lg_model.model_code = this_model
        all_lg_models.append(lg_model)
        model_index.update({this_model:model_count}) ; model_count += 1

        return (all_lg_models, model_index)

class Settings:
	# Standard paths where the files are kept
	home_dir = None
	data_dir = None
	outp_dir = None
	env_type = None
	data_dir = None 
	ahf_path = None
	sub_path = None
	base_path = None
	plot_out = None
	base_run = None
	run_num = None
	sub_path_ahf = None
	sub_path_ics = None
	sub_path_snap = None
	base_out = None

	# Properties of the files and simulation
	resolution = None
	n_ic = None
	n_z0 = None
	n_ahf = 1	
	file_z0 = None
	box_center = [50, 50, 50]
	
	# Redshift values and snapshot numbers for the simulation
	redshifts = []
	snapshots = []

	# Suffixes
	ahf_file = ''
	ahf_halos_str = '.AHF_halos'
	ahf_parts_str = '.AHF_particles'
	snapshot_str = 'snapshot_'

	# Output file names
	file_png_name = None
	file_lg_name = None
	file_lgall_name = None
	file_lare_name = None
	file_codes_name = None

	# Input file names
	file_ahf_in = None
	file_ic_in = None
	file_z0_in = None
	file_ahf = None
	file_ic = None
	file_z0 = None

	def __init__(self, home_dir, outp_dir, env_type, resolution, snap_str):
		self.home_dir = home_dir
		self.outp_dir = outp_dir
		self.env_type = env_type
		self.resolution = resolution
		self.snapshot_str = snap_str

		self.redshifts = []
		self.snapshots = []

		self.base_run = ''
		self.info()
		self.init_paths()
		self.init_ascii()

	def info(self):
		print('Home dir: ', self.home_dir)
		print('Data dir: ', self.data_dir)
		print('Lare out: ', self.file_lare_name)
		print('Env type: ', self.env_type)

	def re_init(self):
		self.init_paths()
		self.init_ascii()

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

	def get_chunk_files(self, isnap, n_files):
		snap_str = '%03d' % isnap
		#self.redshifts = ahf_redshifts(self.ahf_path)
		self.snapshots = ahf_snapshots(self.ahf_path)

		this_snap = self.base_file_chunk + snap_str + '.'
		this_suff_halo = '.z' + self.redshifts[isnap * n_files] + self.ahf_halos_str
		this_suff_part = '.z' + self.redshifts[isnap * n_files] + self.ahf_parts_str

		return (this_snap, this_suff_halo, this_suff_part)


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

	# These files change at different steps, are "dynamically" allocated
	def init_files(self, base_run, run_num):
		self.base_run = base_run
		self.run_num = run_num

		# TODO change all the numbers!!!!! add a use_snap variable to substitute the 054 and 127
		if self.resolution == "1024" or self.resolution == "2048":
			self.base_run = self.base_run + '/' + self.run_num

		if self.env_type == "HESTIA" :
			self.ahf_file = 'HESTIA_100Mpc_512_'+self.base_run+'.127.z0.000.AHF_halos'
			self.file_z0 = 'snapshot_127'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.home_dir + '/HESTIA/LaRe/' + self.base_out
			self.n_ic = 2
			self.n_z0 = 8
			self.file_ahf_in = self.base_path + '/'+ self.sub_path + self.base_run + '/' + self.sub_path_ahf + self.ahf_file
			self.file_z0_in = self.base_path + '/' + self.sub_path + self.base_run + '/' + self.sub_path_snap + self.file_z0
			#self.file_ic_in = self.base_path + '/' + self.sub_path + self.base_run + '/' + self.sub_path_ics + self.file_ic

		elif self.env_type == "512box100" :
			self.ahf_file = self.snapshot_str 
			self.file_z0 = 'snapshot_054'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.base_path + self.resolution + '/' + self.base_run + '/' + self.base_out
			self.n_ic = 1
			self.n_z0 = 1
			self.file_ahf_in = self.base_path + self.sub_path_ahf  + self.base_run + '/' + self.ahf_file
			self.file_z0_in = self.base_path +  self.sub_path_snap + self.base_run + '/' + self.file_z0
			self.file_ic_in = self.base_path +  self.sub_path_ics  + self.base_run + '/' + self.file_ic

		elif self.env_type == "zoom" :
			#self.ahf_file = 'snapshot_054.0000.z0.000.AHF_halos'
			self.ahf_file = self.snapshot_str 
			self.file_z0 = 'snapshot_054'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.base_path + self.resolution + '/' + self.base_run + '/' + self.base_out
			self.n_ic = 1
			self.n_z0 = 1
			self.ahf_path = self.base_path + self.resolution + '/' + self.base_run + '/'
			self.file_ahf_in = self.base_path + self.sub_path_ahf  + self.base_run + '/' + self.ahf_file
			self.file_z0_in = self.base_path +  self.sub_path_snap + self.base_run + '/' + self.file_z0
			#self.file_ic_in = self.base_path +  self.sub_path_ics  + self.base_run + '/' + self.file_ic	#FIXME
	
			# Initialize 
			#self.redshifts = ahf_redshifts(self.ahf_path)
			#self.snapshots = ahf_snapshots(self.ahf_path)

		else:
			self.ahf_file = self.snapshot_str 
			#self.ahf_file = 'snapshot_054.z0.000.AHF_halos'
			self.file_z0 = 'snapshot_054'
			self.base_out = 'lare_z0_' + self.base_run + '.dat'
			self.lare_out = self.base_path + self.resolution + '/' + self.base_run + '/' + self.base_out
			self.n_ic = 2
			self.n_z0 = 1
			self.file_ahf_in = self.base_path + self.resolution + '/' + self.base_run + '/' + self.ahf_file
			self.file_z0_in = self.base_path + self.resolution + '/' + self.base_run + '/' + self.file_z0
			#self.file_ic_in = self.base_path + self.resolution + '/ICs/' + self.file_ic	#FIXME

		# This is the same for all
		self.plot_out =	self.base_path + self.outp_dir + self.base_run + '_particles_LG_LV.png'

	def get_ahf_files(self, i_snap, mpi_task):
		mpi_str = '%04d' % mpi_task

		if self.env_type == "HESTIA" :	# FIXME this needs to be set straight
			this_part_file = self.ahf_path+self.snapshot_str+self.snapshots[i_snap]+'.'+mpi_str+'.z'+self.redshifts[i_snap]+self.ahf_parts_str
			this_halo_file = self.ahf_path+self.snapshot_str+self.snapshots[i_snap]+'.'+mpi_str+'.z'+self.redshifts[i_snap]+self.ahf_halos_str
		else:
			this_part_file = self.ahf_path+self.snapshot_str+self.snapshots[i_snap]+'.'+mpi_str+'.z'+self.redshifts[i_snap]+self.ahf_parts_str
			this_halo_file = self.ahf_path+self.snapshot_str+self.snapshots[i_snap]+'.'+mpi_str+'.z'+self.redshifts[i_snap]+self.ahf_halos_str

		return (this_part_file, this_halo_file)

	def get_png_output(self, info):
		self.file_png_name = self.base_path + self.outp_dir + info + '_' + self.resolution + '_' + self.base_run + '.png'
		return self.file_png_name

	def get_zoom_output(self):
		return (self.file_lg_name, self.file_lgsub_name)

	def get_files_input(self):
		return (self.file_z0, self.base_out, self.n_ic, self.n_z0)
