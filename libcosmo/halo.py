import math 
from find_halos import *
from utils import *

class Halo:
	ID = 1234567890123456789
	m  = 0.0
	x  = [0.0, 0.0, 0.0]
	v  = [0.0, 0.0, 0.0]
	r = 0.0
	nsub = 0
	npart = 0

	def __init__(self):
		self.ID = 1234567890123456789
		self.m  = 0.0
		self.x  = [0.0, 0.0, 0.0]
		self.v  = [0.0, 0.0, 0.0]
		self.r = 0.0
		self.nsub = 0
		self.npart = 0

	def initialize(self, ind, mass, pos, vel, rad, n_sub, n_part):
		self.ID = ind
		self.m  = mass 
		self.x  = pos 
		self.v  = vel
		self.r  = rad 
		self.nsub = n_sub
		self.npart = n_part

	def distance(self, rad):
		dist = math.sqrt( pow((self.x[0] - rad[0]), 2) + pow((self.x[1] - rad[1]), 2) + pow((self.x[2] - rad[2]), 2) )
		return dist

	def info(self):
		return "ID: %ld, M: %.3e, X: (%7.3f, %7.3f, %7.3f), Vel: (%7.3f, %7.3f, %7.3f)" % \
			(self.ID, self.m, self.x[0], self.x[1], self.x[2], self.v[0], self.v[1], self.v[2])



class HaloThroughZ:
	n_steps = 0		# Number of steps available for this halo
	t_step = []		# Time Step in Myrs
	z_step = []		# Time steps in redshift
	n_mergers = []		# Number of mergers happening at that timestep
	halo = []		# List storing a Halo structure per each timestep
	subhalos = []		# List storing a Halo structure per each timestep
	last_major_merger = 0.0
	formation_time = 0.0
	is_smooth = False

	def __init__(self):
		z = []
		n_mergers = []	
		halo_z = []
		LastMajorMerger = 0.0

	def v_t(self):
		vel_t = np.zeros((3, self.n_steps))

		for ixy in range(0, self.n_steps):
			vel_t[:][ixy] = self.halo[ixy].v[:]

		return vel_t


	def x_t(self):
		pos_t = np.zeros((3, self.n_steps))

		for ixy in range(0, self.n_steps):
			pos_t[:][ixy] = self.halo[ixy].x[:]

		return pos_t

	def add_halo(self, halo):
		self.halo.append(halo)

	def add_subhalos(self, subhalos):
		self.subhalos.append(subhalos)

	# TODO do all mmt computation
	def last_major_merger(self):
		return self.last_major_merger
	
	# TODO compute halo formation time
	def formation_time(self):
		return self.formation_time

	# TODO Get rid of flybys and spurious events polluting the halo merger history
	def smooth_history(self):
		self.is_smooth = True
		# TODO do the actual smoothing
		# TODO add a library with Gravity & Hubble expansion etc. to find out grav potentials and bound objects



class LocalGroupModel:
	d_max = 5000. # Box center distance
	d_iso = 2000.
	r_max = 1250.
	r_min = 250.
	m_max = 5.e+12
	m_min = 5.e+11
	mtot_max = 8.e+12
	mratio_max = 4.
	vrad_max = -1.0
	center = [50000.0] * 3
	model_type = '00'

	def __init__(self, d_max, d_iso, r_max, r_min, m_max, m_min, mratio_max, vrad_max):
		self.d_max = d_max
		self.d_iso = d_iso
		self.r_max = r_max
		self.r_min = r_min
		self.m_max = m_max
		self.m_min = m_min
		self.mratio_max = mratio_max
		self.vrad_max = vrad_max

	def info(self):
		lg_mod_info = "D_max: %.3f, D_iso: %.3f, Rh_max: %.3f, Rh_min: %.3f, M_min: %.3f\n" % \
				(self.d_max, self.d_iso, self.r_max, self.r_min, self.m_min)


class SubHalos():
	host = Halo()
	header = ''
	host_name = 'LG'
	hubble = 1.0
	sub = []
	code = '00_00_00'
	n_sub = 0
	
	# This one holds the selected subhalos above a given threshold	
	select_subs = []
	host_coords = []
	n_select_subs = 0
	
	moi_evals = np.zeros((3))
	moi_evecs = np.zeros((3, 3))
	moi_red_evals = np.zeros((3))
	moi_red_evecs = np.zeros((3, 3))

	def __new__(cls, *args, **kwargs):
		instance = super(SubHalos, cls).__new__(cls, *args, **kwargs)
		return instance

	def __init__(self, host, subs, code, host_name):
		self.sub = []
		self.host = host
		self.code = code
		self.host_name = host_name
		self.init_subs(subs)

	def init_subs(self, subs):
		self.n_sub = len(subs)
		self.sub = []

		for ih in range(0, self.n_sub):
			self.sub.append(subs[ih])

	def sub_over_m(self, m_min):
		n_sub_min = 0
		subs_min = []

		for ih in range(0, self.n_sub):
			if self.sub[ih].m > m_min:
				n_sub_min += 1
				subs_min.append(self.sub[ih])

		return subs_min

	def mass_function(self):
		masses = []

		for im in range(0, self.n_sub):
			masses.append(self.sub[im].m)

		(x_m, y_n) = mass_function(masses)
		return (x_m, y_n)

	def sub_over_n(self, n_min):
		n_sub_min = 0
		subs_min = []

		for ih in range(0, self.n_sub):
			if self.sub[ih].npart > n_min:
				n_sub_min += 1
				subs_min.append(self.sub[ih])

		return subs_min

	def anisotropy(self, type_cut, value_cut):
		do_tensor = True

		if type_cut == "mass":
			these_subs = self.sub_over_m(value_cut)
		elif type_cut == "part":
			these_subs = self.sub_over_n(value_cut)
		else:
			do_tensor = False
			print('Wrong kind of cutoff. Only "mass" and "part" tags are possible.')

		if do_tensor == True:
			# First of all rescale all coordinates in the system of the host halo
			#subs_x = []
			subs_x = np.zeros((1, 3))
			this_x = np.zeros((3))
			subs_n = len(these_subs)
			subs_m = np.zeros((subs_n))

			self.select_subs = these_subs
			self.n_select_subs = subs_n
			
			for ih in range(0, subs_n):
				subs_m[ih] = 1.0	
				#subs_m[ih] = these_subs[ih].m

				for ix in range(0, 3):
					this_x[ix] = these_subs[ih].x[ix] - self.host.x[ix]
				
				if ih == 0:
					subs_x[0] = this_x
				else:
					subs_x.resize((ih+1, 3))
					subs_x[ih] = this_x
			self.host_coords = subs_x

			# Automatically compute all types of tensors
			print 'Inertia tensor computed using %d subhalos.' % subs_n
			(self.moi_evals, self.moi_evecs) = moment_inertia(subs_x, subs_m)
			
			print 'Reduced inertia tensor computed using %d subhalos.' % subs_n
			(self.moi_red_evals, self.moi_red_evecs) = moment_inertia_reduced(subs_x, subs_m)

	def basis_eigenvectors(self, evec_type):
		change_basis = True

		if evec_type == "inertia":
			use_evec = self.moi_evecs

		elif evec_type == "inertia_reduced":
			use_evec = self.moi_red_evecs

		else:
			change_basis = False
			print 'Error. Select "inertia" or "inertia_reduced". '

		if change_basis == True:
			print 'Computing subhalo positions in the %s tensor eigenvector basis.' % evec_type
			new_coords = np.zeros((self.n_select_subs, 3))

			for i_s in range(0, self.n_select_subs):
				for i_x in range(0, 3):
					prod = dot_prod(use_evec[i_x], self.host_coords[i_s])
					new_coords[i_s][i_x] = prod 
		return new_coords

	def header(self):
		self.header = '# ID_sub(1)    M_sub(2)[Msun]  R_sub(3)[Mpc]   Host(4)         SimuCode(5)\n'
		return self.header

	def all_info(self, type_cut, value_cut): 
		file_sub_line = ''
		print_lines = True
		h0 = self.hubble

		if type_cut == "mass":
			these_subs = self.sub_over_m(value_cut)
		elif type_cut == "part":
			these_subs = self.sub_over_n(value_cut)
		else:
			print_lines = False
			print('Wrong kind of cutoff. Only "mass" and "part" tags are possible.')

		if print_lines == True:
			n_print = len(these_subs)
			
			for il in range(1, n_print):
				this_r = self.host.distance(these_subs[il].x)

				line = '%ld    %.2e    %7.2f   %s      %s\n' % \
 	      			 (these_subs[il].ID,  these_subs[il].m/h0, this_r/h0, self.host_name, self.code)
				file_sub_line += line

		return file_sub_line



class LocalGroup:
	code = '00_00_00'
	c_box = [50000.0] * 3
	ahf_file = 'this_file.AHF_halos'
	vrad = -100.
	r = 770.
	virgo_x = [47500., 61000., 49500.]
	d_cbox = 0.0
	d_virgo = 0.0
	rating = 0.0
	com = [0.0] * 3
	hubble = 0.67

	LG1 = Halo()
	LG2 = Halo()

	def __new__(cls, *args, **kwargs):
		instance = super(LocalGroup, cls).__new__(cls, *args, **kwargs)
		return instance

	def __init__(self, code):
		self.code = code

	def init_halos(self, lg1, lg2):

		if lg1.m > lg2.m:
			self.LG1 = lg1
			self.LG2 = lg2
		else:
			self.LG1 = lg2
			self.LG2 = lg1

		if lg1.x[0] != 0 and lg2.x != 0:
			self.r = self.r_halos()
			self.vrad = self.v_radial()
			self.com = self.com()

	def rating(self):
		self.rating = rate_lg_pair(self.LG1, self.LG2, self.c_box)

	def com(self):
		self.com = center_of_mass([self.LG1.m, self.LG2.m], [self.LG1.x, self.LG2.x])
		return self.com 

	def lg_member(self, n_member):
		if n_member == 0:
			return self.LG1
		if n_member == 1:
			return self.LG2

	def r_halos(self):
		self.r = distance(self.LG1.x, self.LG2.x)
		return self.r

	def v_radial(self):
		self.vrad = vel_radial(self.LG1.x, self.LG2.x, self.LG1.v, self.LG2.v)
		self.vrad += self.hubble * self.r * 0.1
		return self.vrad

	def m_ratio(self):
		if self.LG1.m > self.LG2.m:
			return (self.LG1.m / self.LG2.m)
		else:
			return (self.LG1.m / self.LG2.m)
	def header(self):
		header = '# ID_M31(1) ID_MW(2) M_M31(3)[Msun] M_MW(4)[Msun] R_MwM31(5)[Mpc] Vrad_M31(6)[phys, km/s] Nsub_M31(7) Nsub_MW(8) SimuCode(9) Xcom(10) Ycom(11) Zcom(12)\n'
		return header

	def info(self):
		h0 = self.hubble 
		kpc = 1000.
		file_lg_line = '%ld   %ld   %7.2e   %7.2e   %7.2f   %7.2f   %5d   %5d  %s  %5.2f  %5.2f  %5.2f\n' % \
			(self.LG1.ID, self.LG2.ID, self.LG1.m/h0, self.LG2.m/h0, self.r/h0, self.vrad, \
				self.LG1.nsub, self.LG2.nsub, self.code, self.com[0]/kpc, self.com[1]/kpc, self.com[2]/kpc)

		return file_lg_line

