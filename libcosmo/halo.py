import math 
import pickle
import numpy as np
import find_halos as fh
from scipy import interpolate
from operator import *
from utils import *
from units import *

class Halo:
	ID = 1234567890123456789
	m  = 0.0
	x  = [0.0, 0.0, 0.0]
	v  = [0.0, 0.0, 0.0]
	l  = [0.0, 0.0, 0.0]
	r = 0.0
	nsub = 0
	rsub = 0.0	# This is supposed to be Rvir but it can be slightly larger, to look for nearby bound halos
	subhalos = [] 	# This one only contains all sub-halo IDs
	npart = 0

	def __init__(self):
		self.ID = 1234567890123456789
		self.m  = 0.0
		self.x  = [0.0, 0.0, 0.0]
		self.v  = [0.0, 0.0, 0.0]
		self.r = 0.0
		self.nsub = 0
		self.npart = 0

	def __getitem__(self, key):
		if key == 0:
			return self.ID
		elif key == 1:
			return self.m
		elif key == 2:
			return self.r
		else:
			return self.nsub

	def initialize(self, ind, mass, pos, vel, rad, n_sub, n_part):
		self.ID = ind
		self.m  = mass 
		self.x  = pos 
		self.v  = vel
		self.r  = rad 
		self.npart = n_part

	def distance(self, rad):
		dist = math.sqrt( pow((self.x[0] - rad[0]), 2) + pow((self.x[1] - rad[1]), 2) + pow((self.x[2] - rad[2]), 2) )
		return dist

	def sub_halos(self, all_halos, rsub):
		self.rsub = rsub
		subs = fh.find_halos_point(self.x, all_halos, self.rsub)
		self.nsub = len(subs) - 1

		# WARNING TODO this assumes that the first halo is the HOST halo
		for isub in range(1, self.nsub+1):
			self.subhalos.append(subs[isub])

	def sort_sub_mass(self):
		sort_mass = sorted(self.subhalos, key=itemgetter(1), reverse=True)

		#for isub in range(0, self.nsub):
		#	print isub, sort_mass[isub].m, sort_mass[isub].x

		return sort_mass

#	def sort_sub_dist(self):


	def info(self):
		#return "ID: %ld, M: %.3e, X: (%7.3f, %7.3f, %7.3f), Vel: (%7.3f, %7.3f, %7.3f)" % \
		#	(self.ID, self.m, self.x[0], self.x[1], self.x[2], self.v[0], self.v[1], self.v[2])
		return "%ld %.3e %.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d" % \
			(self.ID, self.m, self.r, self.x[0], self.x[1], self.x[2], self.v[0], self.v[1], self.v[2], self.nsub, self.npart)



class HaloThroughZ:
	halo = []		# List storing a Halo structure per each timestep
	t_step = []		# Time Step in Myrs
	z_step = []		# Time steps in redshift
	subhalos = []		# List storing a Halo structure per each timestep
	n_mergers = []		# Number of mergers happening at that timestep

	n_steps = 0		# Number of steps available for this halo
	is_smooth = False
	formation_time = 0.0
	last_major_merger = 0.0

	def __init__(self, n_steps):
		self.n_steps = n_steps
		self.halo = []
		self.t_step = []		
		self.z_step = []		
		self.subhalos = []		
		self.n_mergers = []	
		self.formation_time = 0.0
		self.last_major_merger = 0.0

	def m_t(self):
		mass_t = np.zeros((self.n_steps3))

		for ixy in range(0, self.n_steps):
			mass_t[ixy] = self.halo[ixy].m

		return mass_t


	def v_t(self):
		vel_t = np.zeros((3, self.n_steps))

		for ixy in range(0, self.n_steps):
			for iv in range(0, 3):
				vel_t[iv][ixy] = self.halo[ixy].v[iv]

		return vel_t

	def x_t(self):
		pos_t = np.zeros((3, self.n_steps))

		for ixy in range(0, self.n_steps):
			for ip in range(0, 3):
				pos_t[ip][ixy] = self.halo[ixy].x[ip]

		return pos_t

	def add_step(self, halo, t, z):
		self.halo.append(halo)
		self.t_step.append(t)
		self.z_step.append(z)

	def add_subhalos(self, subhalos):
		self.subhalos.append(subhalos)

	# TODO do all mmt computation
	def last_major_merger(self):
		m_merg = 0.1 

		#for im in range(self.n_steps+1, 0, -1):
		for jm in range(0, self.n_steps-1):
			im = self.n_steps -jm -1

			m0 = self.halo[im].m
			m1 = self.halo[im-1].m
			z0 = self.z_step[im]
			z1 = self.z_step[im-1]
			
			dM = abs(m0 - m1) / m0
			dZ = 0.5 * (z0 + z1)
			#print im, im-1, z0, z1, m0, m1, dM
			
			if dM > m_merg:
				lmmz = dZ 
				itime = im

		self.last_major_merger = z2Myr(lmmz)

		#print im, lmmz, self.last_major_merger
		return self.last_major_merger
	
	# Use Z in computations - will be helpful if will be needed to use mass(z) interpolation functions
	def formation_time(self):
		m_half = 0.5 * self.halo[0].m
		im = 0
		m0 = 1.e+16

		while m0 > m_half:
			m0 = self.halo[im].m
			z0 = self.z_step[im]
			# Only add a step if the object is above the threshold
			if m0 > m_half:
				im += 1

		form_z = 0.5 * (self.z_step[im] + self.z_step[im-1])
		self.formation_time = z2Myr(form_z)

		#print m0, z0, self.formation_time
		return self.formation_time

	def dump_history(self, f_name):
		f_out = open(f_name, 'wb')
		header = "# z(0) ID(1) M(2) r(3) x(4) y(5) z(6) vx(7) vy(8) vz(9) nsub(10) npart(11)\n"
		f_out.write(header)

		for ih in range(0, self.n_steps):
			line_z = "%5.3f " % (self.z_step[ih])
			#print ih, line_z, self.halo[ih].info(), "\n"
			line = line_z + self.halo[ih].info() + "\n"
			f_out.write(line)

	def load_file(self, f_name):
		f_in = open(f_name, 'r')
		line = f_in.readline()

		while line:
			line = f_in.readline()
			line = line.strip()
			column = line.split()
			halo = Halo()
			n_col = len(column)
			
			if n_col > 0:
				z_step = float(column[0])
				hid = long(column[1])
				t_step = z2Myr(z_step)
				mass = float(column[2])
				pos = [float(column[4]), float(column[5]), float(column[6])]
				vel = [float(column[7]), float(column[8]), float(column[9])]
				rvir = float(column[3])
				nsub = int(column[10])
				npart = int(column[11])
				halo.initialize(hid, mass, pos, vel, rvir, nsub, npart)
				self.add_step(halo, t_step, z_step)

	# TODO Get rid of flybys and spurious events polluting the halo merger history
	def smooth_history(self):
		self.is_smooth = True
		# TODO do the actual smoothing
		# TODO add a library with Gravity & Hubble expansion etc. to find out grav potentials and bound objects


class SubHaloThroughZ(HaloThroughZ):
	
	host = HaloThroughZ(0)

	# In principle we allow for an halo to cross the viral radius several times - there might be more accretion times/steps
	acc_time = []
	acc_step = []

	def __init__(self, n_steps):
		self.n_steps = n_steps
		self.halo = []
		self.t_step = []		
		self.z_step = []		
		self.subhalos = []		
		self.n_mergers = []	
		self.formation_time = 0.0
		self.last_major_merger = 0.0

	def x_t_host_center(self):
		pos_t = np.zeros((3, self.n_steps))

		for ixy in range(0, self.n_steps):
			for ip in range(0, 3):
				pos_t[ip][ixy] = self.halo[ixy].x[ip] - self.host.halo[ixy].x[ip]

		return pos_t


	def accretion_time(self):
		d_old = 0.0
		r_old = 0.0
		n_passage = 1

		for iat in range(0, self.n_steps):
			try:	
				d_host = self.halo[iat].distance(self.host.halo[iat].x) 
				r_host = self.host.halo[iat].r 
				m_sub = self.halo[iat].m 
				
				if iat > 0:
					acc_z = 0.5 * (self.z_step[iat] + self.z_step[iat-1])

				if iat > 0 and (d_host < r_host and d_old > r_old) or (d_host > r_host and d_old < r_old):
					self.acc_step.append(iat)
					self.acc_time.append(z2Myr(acc_z))

				d_old = d_host		
				r_old = r_host 

			except:
				print 'Error at step ', iat

		return self.acc_time

	# Mass at accretion time - n_time should be 0 but if the halo is accreted and ejected several times n_time could be larger
	def max_mass(self, n_time):
		self.accretion_time()
		return self.halo[self.acc_step[n_time]].m



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

# Class for bulk properties of sub halo groups of a given halo - mass functions, spatial distribution, etc.
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

