#from libpynbody.read_ahf import *
from libcosmo.halos import *
#from libio.read_ascii import *
from pygadgetreader import *
#import gadfly as gdf

import numpy as np
import pynbody as pyn
import pynbody.plot.sph as sph
import matplotlib.pylab as plt
import pickle

#snap_path='/home/eduardo/CLUES/DATA/CF2P5/57252/snapshot_019'
#snap_path='/home/eduardo/CLUES/DATA/01/snapshot_035'
snap_path='/home/eduardo/CLUES/DATA/HESTIA/2048/37_11/snapshot_127'
#snap_path='/home/eduardo/CLUES/DATA/1024/00/snapshot_054'
#ahf_path='/home/eduardo/CLUES/DATA/1024/00/snapshot_054.z0.000.AHF_halos'
#sim = gdf.Simulation(snap_path)
#snap = sim.load_snapshot('_054')

s = pyn.load(snap_path)
s.physical_units()
h = s.halos() 

#print(s._block_list)

#mw = h[2]
mw = h[151]

#print(mw['pos'][0])
print(mw.properties['n_gas'])
print(mw.properties['n_star'])
print(mw.properties['mass']/1.e+12)
print(mw.properties['M_star']/1.e+12)

plt.ioff()
#pyn.analysis.angmom.faceon(mw)
#cen = pyn.analysis.halo.center(mw, mode='hyb', retcen=True)
#sph.image(mw.g,qty="rho",units="g cm^-3",width=100,cmap="Greys")
#pyn.plot.image(s.d[pyn.filt.Sphere('10 Mpc')], width='30 Mpc', units = 'Msol kpc^-2', cmap='Greys');
#pyn.plot.image(s.d[pyn.filt.Sphere('2000 kpc')], width='10000 kpc', units = 'Msol kpc^-2', cmap='Greys');

#print(mw.g)
#print(mw.g[0].properties)
#pyn.plot.image(mw.g, width=100, cmap='Blues')

plt.draw()
plt.savefig('test.png')


#f, axs = plt.subplots(1,2,figsize=(14,6))
#sph.velocity_image(mw.g, vector_color="cyan", qty="temp",width=50,cmap="YlOrRd",
#                   denoise=True,approximate_fast=False, subplot=axs[0], show_cbar = False)
#sph.image(mw.g,qty="rho",width=50,cmap="YlOrRd", denoise=True,approximate_fast=False)
#pyn.analysis.angmom.sideon(mw)
#pyn.plot.stars.render(s,width='20 kpc')
#pyn.plot.image(mw.g, width=100, cmap='Blues', threaded=False)
#pyn.plot.image(s.d[pynbody.filt.Sphere('10 Mpc')], width='10 Mpc', units = 'Msol kpc^-2', cmap='Greys');
#snap_path='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054'
#snap_path='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054'
#snap_path='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_h78/70_00/snapshot_127'
#ahf_path='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
#ahf_path='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_'
#ahf_path='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.AHF_halos'
#snap_path='/home/eduardo/CLUES/pynbody/testdata/g15784.lr.01024.gz'
#ahf_file = pyn.load(ahf_path)
#snap_path='/home/eduardo/CLUES/DATA/CF2P5/59050/snapshot_019'
#s.physical_units()
#print(s.families())
#for i in range(0, 5):
#    print(i, h[i].properties['halo_id'], h[i].properties['Mvir'])
#pyn.analysis.halo.center(h[1], mode='hyb')
#pyn.plot.sph.image(s, qty="rho", units="rho", width=100, cmap="Greys")
#halos = pyn.halo.ahf.AHFCatalogue._load_ahf_halos(s, ahf_path)
#halos = pyn.halo.ahf.AHFCatalogue(s, pathname=ahf_path)
#halos = pyn.halo.AHFCatalogue(s)
#halos._load_ahf_halos(ahf_path+'halos')
#h = s.halos(ahf_basename=ahf_path)
#h = pyn.load(ahf_path)
#print('gas=%d, dark=%d' % (len(s.gas), len(s.dark)))
#h = pyn.halo.ahf.AHFCatalogue('sim', ahf_basename=ahf_path)
#h = pyn.halo.ahf.AHFCatalogue(s, ahf_basename=ahf_path)
#h = pyn.halo.ahf.AHFCatalogue('sim', ahf_path)



'''
# THIS ARE THE OLD PLOTTING ROUTINES THAT USE SEABORN
elif oldLoad:

    doBar=False
    #doBar=True
    n_levels = 10

    # Small box settings
    #contour = np.logspace(-6.0, -1.0, num=20)

    # Full box
    contour = np.logspace(-6.5, -4.0, num=8)

    print(snap)
    print('Loading files: ', slab_x_fname)
    x_out = open(slab_x_fname, 'rb')
    y_out = open(slab_y_fname, 'rb')
    all_data_x = pickle.load(x_out)
    all_data_y = pickle.load(y_out)

    data_x = [];     data_y = []

    for ix in range(0, len(all_data_x)):
        dx = all_data_x[ix]
        dy = all_data_y[ix]

        if abs(dx) < side and abs(dy) < side:
            data_x.append(dx)
            data_y.append(dy)

    print('Done. Selected ', len(data_x), ' particles out of ', len(all_data_x))
    figname=snapname+simname+'_'+code+'_' + str(gridSize) + '_sm' + str(smooth) + '_fac' + str(reduce_fac) + '_side' + str(side) + '.png'

    plt.figure(figsize=(20, 20))
    matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel("SGY (Mpc/h)")
    plt.xlabel("SGX (Mpc/h)")

    print('Smoothing out plot on a ', gridSize, ' grid with a ', smooth, ' Gaussian kernel.')
    sns.kdeplot(data_x, data_y, cmap="coolwarm", shade=True, shade_lowest=True, gridsize=gridSize,
        bw=smooth, levels=n_levels, cbar=doBar)
    # SCATTERPLOT just in case
    #plt.plot(data_x, data_y, linestyle='', marker='o', markersize=0.1)
    #plt.savefig('test_scatter.png', dpi=100)
    plt.tight_layout()
    print('Saving output to file: ', figname)
    plt.savefig(figname, dpi=300)

'''



'''
f_ahf='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
halos = read_ahf(f_ahf)
for ih in range (0, len(halos)):
	mass = halos[ih].m/1e+8
	print  mass, ' ', halos[ih].vmax
#fname = 'saved/dyn_stats_MW_00_06_01.pkl'
#fname='saved/web_00_06_01.pkl'
#filesub = open(fname, 'r')
#sats = pickle.load(filesub)
#print sats
#random_table_triaxialities(10, 1000, False)

for npts in range(3, 200):
	#rand = random_table_triaxialities(10, 1000, True)
	random_table_triaxialities(npts, 10000, False)

#print rand

#print rand_points_sphere(16)

#(evs, disps, perc) = random_triaxialities(10, 1000, 0.25)

#fopen = open('saved/angles_align_satellite_anisotropy_00_LG.pkl', 'r')
#angs = pickle.load(fopen)
#print angs

#print evs, perc


vec = [0.0, 1.0, 0.0]

(e, d, c) = random_triaxialities_and_angles(20, 1000, vec)

print c

basedir='/home/eduardo/CLUES/DATA/512/'
snapfile=basedir+'70_00/snapshot_054'; ahffile=basedir+'70_00/snapshot_054.AHF_halos'; rho_out='rho_70_00_EC.png'
#snapfile=basedir+'60_00/snapshot_054'; ahffile=basedir+'60_00/snapshot_054.AHF_halos'; rho_out='rho_60_00_EC.png'
#snapfile=basedir+'NewParams/60_01/snapshot_054'; ahffile=basedir+'NewParams/60_01/snapshot_054.0000.z0.000.AHF_halos'; rho_out='rho_60_01_EC_new.png'
#snapfile=basedir+'/NOAM/60_00/output/snapdir_127/snapshot_127'; ahffile=basedir+'NOAM/60_00/snapshot_127.0000.z0.000.AHF_halos'; rho_out='rho_60_00_NIL.png'
#snapfile=basedir+'/NOAM/60_01/output/snapdir_127/snapshot_127'; ahffile=basedir+'NOAM/60_01/AHF_output/HESTIA_100Mpc_512_60_01.127.z0.000.AHF_halos'; rho_out='rho_60_01_NIL.png'
#basedir='/home/eduardo/CLUES/DATA/2048/'
 
halo_all = read_ahf(ahffile)

rad = 15000.
rescale = 2

#virgo_x = [61000.0, 49500.0, 47000]
virgo_x = [47000.0, 61000.0, 49500.0]
#virgo_x = [49000.0, 47000.0, 61500.0]
center = [50000.0, 50000.0, 50000.0]
#find_halos_point(x_c, halo_all, radius):
halos_c = find_halos_point(center, halo_all, rad)

lg1 = halos_c[0]
lg2 = halos_c[1]

n_lgs = len(halos_c)

#for i_lg in range(0, n_lgs):
#	print halos_c[i_lg].info()


#print lg1.info()
#print lg2.info()

#plot_lglv(settings.file_z0_in, ahf_all, settings.plot_out, best_lg.LG1, best_lg.LG2, x0, rescale, 2)
plot_lglv(snapfile, halo_all, rho_out, lg1, lg2, virgo_x, rescale, 2)


webfile='/home/eduardo/CLUES/DATA/2048/00_06/00/zoom_2048_054.000064.Vweb-ascii'
ahffile='/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'

halos = read_ahf(ahffile)

index='5260986610347328281'

print halos[0].id_index[index]

#grid = read_vweb(webfile, 64, 100)

pos = [50, 50, 50]

print grid.rhoPos(pos)
print grid.velPos(pos)
print grid.evalPos(pos)
print grid.evecPos(pos)

	
grid = Grid(256, 100)


ijk = grid.phys2grid(pos)

print ijk
print grid.grid2phys(ijk)


substart = 0
subend = 10

satfname='/home/eduardo/CLUES/PyRCODIO/saved/00_06_'+subrun+'_sats.pkl'
mainfname='/home/eduardo/CLUES/PyRCODIO/saved/00_06_'+subrun+'_mains.pkl'
fhandsat=open(satfname, 'r')
fhandmain=open(mainfname, 'r')

mains = pickle.load(fhandmain)
ssats = pickle.load(fhandsat)

print len(ssats[1][:])#.anisotropy("part", 20)

for i in range(0, 40):
	#print ssats[0][i].anisotropy('part', 10)
	
	(vals, red_vals, it_vals) = ssats[1][i].anisotropy('part', 30)

	#print 'Step ', i
	#print 'Vals   : ', vals[0]/vals[2], vals[1]/vals[2]
	#print 'RedVals: ', red_vals[0]/red_vals[2], red_vals[1]/red_vals[2]
	#print 'ItVals : ', it_vals[0]/it_vals[2], it_vals[1]/it_vals[2]
	#print ''	

#print mains
#print len(ssats[1])
#print (ssats[0])
#print ssats[0].host[0].info()
#print ssats
#print mains[0].halo[0].info()


ids = [0] * 6
#ids = np.zeros((6))
ids[0] = 10283893
ids[1] = 123456789
ids[2] = 123456789
ids[3] = 18923838
ids[4] = 123456789
ids[5] = 123456789


#ids = [10283747492921, 123456789, 1903784749321]
this = 123456789

#ids = [0, 1, 2, 3, 4, 5, 6, 7]

pts = np.where(ids == this)
#pts = np.where(ids > 3 and ids < 5)

print pts[0], pts

ahf_file='/store/clues/HESTIA/RE_SIMS/4096/GAL_FOR/55_02/AHF_output/HESTIA_100Mpc_4096_55_02.127.z0.000.AHF_halos'

(lg_model, lg_dict) = lg_models()

this_index = int(lg_dict['55_02'])
this_model = lg_model[this_index]

print lg_model[0].info()
print lg_model[1].info()
print lg_model[2].info()
print lg_model[4].info()

print lg_dict
print this_index, this_model.info(), this_model

ahf = read_ahf(ahf_file)
hs = find_lg(ahf, this_model)

print hs[0].info()
print hs[1].info()

print hs[2].info()
print hs[3].info()

#base_path = '/home/eduardo/CLUES/DATA/1024/00_06/LV/'
#snap_file = base_path + 'snapshot_020'
#ahf_file = base_path + 'snapshot_020.0000.z0.000.AHF_halos.pttype1'
ahf_snap='/home/eduardo/CLUES/DATA/QL/1024/00_00/snapshot_054.AHF_halos'
center = [250000.0, 250000.0, 250000.0]
this_halos = read_ahf(ahf_snap)
(pos, mass, names) = locate_clusters(this_halos, center)
#print pos
#print mass
#print names

n = len(pos)

for i in range(0, n):
	pos[i][0] -= 250000.
	pos[i][1] -= 250000.
	pos[i][2] -= 250000.

	#pos[i][0] /= 0.677
	#pos[i][1] /= 0.677
	#pos[i][2] /= 0.677

	mass[i] /= 0.677
	print 'Name:%s,  Mass:%e, Pos: %s ' % (names[i], mass[i], pos[i])


ahf_part1 = '/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_particles'
ahf_part2 = '/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_053.0000.z0.017.AHF_particles'

ahf_halo1 = '/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
ahf_halo2 = '/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_053.0000.z0.017.AHF_halos'

(ids1, parts1) = read_particles(ahf_part1)
#(ids2, parts2) = read_particles(ahf_part2)

ahf1 = read_ahf(ahf_halo1)
#ahf2 = read_ahf(ahf_halo2)

base_path = '/home/eduardo/CLUES/DATA/2048/00_06/00/'
root_file = 'snapshot_'
suff_ahf = '.AHF_halos'
suff_part = '.AHF_particles'

h1 = []
p1 = []
id1 = []

h1.append(ahf1[3])
p1.append(parts1[3])
id1.append(ids1[3])

ini_snap = 50
end_snap = 54
min_common = 15

merger_tree(end_snap, ini_snap, base_path, root_file, suff_part, suff_ahf, min_common, h1, p1, id1)

#find_progenitors
#id1 = parts1[0]
#id2 = parts2[0]
#common = compare_particles(id1, id2)
#print 'found particles: ', common






ahf_file = '/home/eduardo/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
ini_snap = 0
end_snap = 54


out_path='output/'
num_run = '00_06_00'

(n_host, n_subs) = n_mah_files(out_path, num_run)

print n_host
print n_subs

all_halo = []
all_sub = []
all_sub_x = []
all_sub_y = []

for i in range(0, 1):
	this_halo = HaloThroughZ(end_snap - ini_snap)
	halo_name = mah_main_file_name(out_path, num_run, i)
	this_halo.load_file(halo_name)
	all_halo.append(this_halo)
	tmp_subs = []

	for s in range(0, 10):
		this_sub = SubHaloThroughZ(end_snap - ini_snap)
		sub_name = mah_sub_file_name(out_path, num_run, i, s)
		this_sub.host = this_halo
		this_sub.load_file(sub_name)
		#print i, s, this_sub.accretion_time()
		sub_x = this_sub.x_t_host_center()
		all_sub_x.append(sub_x[0][:])
		all_sub_y.append(sub_x[1][:])
		tmp_subs.append(this_sub)

	all_sub.append(tmp_subs)


plot_trajectory(all_sub_x, all_sub_y, 'x', 'y', out_path + 'test.png')

'''

'''
halo = []

halo.append(ahf[0])
halo.append(ahf[1])
halo.append(ahf[2])
halo.append(ahf[3])
rsub = 1.33

(all_halo, all_sub) = halos_and_subhalos_through_z(end_snap, ini_snap, base_path, root_file, ahf_suff, halo,  rsub)

n_halo = len(all_halo)

for ih in range(0, n_halo):
	file_name = mah_main_file_name(out_path, num_run, ih)
	
	all_halo[ih].dump_history(file_name)

	n_sub = len(all_sub[ih])

	for isb in range(0, n_sub):
		this_isb = '%02d' % isb
		file_name = mah_sub_file_name(out_path, num_run, ih, isb)
		all_sub[ih][isb].dump_history(file_name)
'''


'''
#halo_z = HaloThroughZ(53)
#halo_z.load_file('test.txt')
base_path = '/home/edoardo/CLUES/DATA/SIMULATIONS/LGF/2048/00_06/00/'
base_path='/home/eduardo/CLUES/DATA/2048/00_06/00/'
#base_path = '/home/eduardo/CLUES/DATA/'
#ahf_file = '/home/eduardo/CLUES/DATA/1024/00_06/02/snapshot_054.0000.z0.000.AHF_halos'
#ahf_file = '/home/eduardo/CLUES/DATA/1024/00_06/02/snapshot_054.0000.z0.000.AHF_halos'
#ahf_file = base_path + '/512/17_01/snapshot_054.0000.z0.000.AHF_halos'
#snap_file = base_path + '/512/17_01/snapshot_054'
#snap_file = base_path + '512/00_06_LV/snapshot_019'
#ahf_file = base_path + '512/00_06_LV/snapshot_019.0000.z0.000.AHF_halos'
#ahf_file = base_path + '/2048/00_06/0/snapshot_054.0000.z0.000.AHF_halos'
#snap_file = base_path + '/2048/00_06/01/snapshot_054'
root_file = 'snapshot_'
suff_file = '.AHF_halos'
#ascii_file = '/home/eduardo/CLUES/DATA/output/lg_candidates_2048_00_06.txt'
#ahf_file = base_path + '/512/17_01/snapshot_054.AHF_halos'
#ahf_file = base_path + '/512/17_00/snapshot_054.z0.000.AHF_halos'
#ahf_file = base_path + 'snapshot_020.0000.z0.000.AHF_halos'

ini_snap = 0
end_snap = 54
counter = 0
timeStep = 250.


print halo_z.v_t()
halo_z.dump_history('test.txt')
acc_time = sub_z.accretion_time()

print acc_time

all_x = []
all_y = []

halo_x = halo_z.x_t()
sub_x = sub_z.x_t()

all_x.append(halo_x[0][:])
all_x.append(sub_x[0][:])
all_y.append(halo_x[1][:])
all_y.append(sub_x[1][:])

#print all_x

plot_trajectory(all_x, all_y, 'x', 'y', 'test.png')
'''	
#		print i_snap, old_halo.info()
#		expe_x = backward_x(this_halo, this_halo, 250.0)
#		print halo_z.x_t()

#print 'FT: ', halo_z.formation_time() / 1000.
#print 'MT: ', halo_z.last_major_merger() / 1000.

'''
	if counter == 0:
		old_halo = this_halos[0]
	else:
		print i_snap, expe_x
		#print 'D_true: %.3f, D_expe: %.3f' % (distance(old_halo.x, this_halos[0].x), distance(expe_x, this_halos[0].x))
		#print angle(old_halo.v, this_halos[0].v)
		print i_snap, 'Old: ', old_halo.info()
		print i_snap, 'New: ', this_halos[0].info()
		print vec_subt(expe_x, old_halo.x)
		old_halo = this_halos[0]
	
	counter += 1



x = [50000., 50000., 50000.]
v = [400., -200., 100]

dMyrs = 1000.0

nx = new_x(x, v, dMyrs)
d = distance(nx, x)

print(nx)
print(d)

ahf = read_ahf(ahf_file)

center = [50000., 50000., 50000.]
d_max = 7000.
r_iso = 2000.
r_max = 1500.
r_min = 250.
ratio_max = 2.5
vrad_max = 0.0
m_min = 3.e+11
m_max = 5.e+12


reduce_fac = 8
lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)

#hs = []
#hs.append(ahf[0])
#hs.append(ahf[2])

#hs = find_halos_point(center, ahf, radius)
hs = find_lg(ahf, lg_model)
subs = find_halos(hs[0], ahf, radius)


(x0, m0, mt) = locate_virgo(ahf)

print x0, m0, mt

for ih in range(0, len(hs)):
	print hs[ih].info()

plot_lglv(snap_file, ahf, 'output_00_06_LV.png', hs[2], hs[3], x0, 8, 1)

center = [50000.0] * 3
rads = 7000.

halos = find_halos_point(center, ahf, rads)

nh = len(halos)

for ih in range(0, nh):
	if halos[ih].m > 4.e+11 :
		print halos[ih].info()

n_pts = 10

points = rand_points_sphere(n_pts, [0.0, 0.0, 0.0], 5.0)
masses = np.zeros((n_pts))

for im in range(0, n_pts):
	masses[im] = 1.0

(evals, evecs) = moment_inertia(points, masses)

print evals
nh = len(hs)

nlg = nh / 2

print 'Nhalos in center: ', nh
m512 = simu_pmass(100, 2048)
print 'hires pmass: %e' % m512

for ih in range(0, nlg):
	lg1 = hs[2 * ih]
	lg2 = hs[2 * ih + 1]
	lgcom = center_of_mass([lg1.m, lg2.m], [lg1.x, lg2.x])
	d12 = distance(lg1.x, lg2.x)
	dcom = distance(center, lgcom)
	rate = rate_lg_pair(lg1, lg2, center)
	print lg1.info()
	print lg2.info()

	print 'Com     : ', lgcom
	print 'Centerbo: ', dcom
	print 'Distance: ', d12 / 0.67
	print 'V_radial: ',  vel_radial(lg1.x, lg2.x, lg1.v, lg2.v) + d12 * 0.067
	print '%.3e, %d, %.5e, %.3f' % (lg1.m, lg1.npart, (lg1.m/lg1.npart), lg1.distance(center))
	print '%.3e, %d, %.5e, %.3f' % (lg2.m, lg2.npart, (lg2.m/lg2.npart), lg2.distance(center))
	print 'Rated: %f\n' % rate

plot_pos = "false"
ptype = 1

lgs = read_lgs(ascii_file)

(bin_m, bin_n) = bin_lg_sub(lgs)

print bin_m
print bin_n

#plot_lg(snap_file, 'output_lg.png', hs[0], hs[1], reduce_fac, ptype, plot_pos)
#lgcom = center_of_mass([hs[0].m, hs[1].m], [hs[0].x, hs[1].x])

#print hs[0].info()

#ptypes = 4

#(xv, mv, mvtot) = locate_virgo(ahf)

#plot_lglv(snap_file, ahf, "output.png", hs[0], hs[1], xv, reduce_fac, ptypes)

#red = readhead(snap_file, 'redshift')
#time = readhead(snap_file, 'time')
#dmPart1 = readsnap(snap_file, 'pos', 1)
#dmPart2 = readsnap(snap_file, 'pos', 2)

array = np.zeros((3, 100), dtype='float')

rType1 = 5000.
rType2 = 15000.
shellZ = 2000.

print len(dmPart1)

npymp = 12

ps1 = select_particles(lgcom, dmPart1, rType1, facMpc)

#ps2 = select_particles(lgcom, dmPart2, rType2, facMpc)

(xa, xb) = find_slice(ps1, 2, hs[0].x[2], shellZ)
t3 = time.clock()
print 'Time past: ', (t3-t2)

xamax = np.amax(xa)
xamin = np.amin(xa)
xbmax = np.amax(xb)
xbmin = np.amin(xb)

shellXY = 1500.

xamax = lgcom[0] + shellXY
xamin = lgcom[0] - shellXY
xbmax = lgcom[1] + shellXY
xbmin = lgcom[1] - shellXY

facX = 10000.
x_min = facX * math.floor(xamin / facX)
x_max = facX * math.ceil(xamax / facX)
y_min = facX * math.floor(xbmin / facX)
y_max = facX * math.ceil(xbmax / facX)

x_min = xamin 
x_max = xamax 
y_min = xbmin
y_max = xbmax

print xamax, xamin
print xbmax, xbmin

ptsize = 0.1

plt.axis([x_min, x_max, y_min, y_max])
plt.scatter(xa, xb, s=ptsize, c=2)
plt.show()


#print ps1
#h0 = ahf[0]
#radius = h0.r * 2.
#hs = find_halos(h0, ahf, radius)
#hs = find_lg(ahf, center, radius, m_min, r_min, r_max)
n_m = 6
masses = [1.0] * n_m

coords = [[0] * 3] * n_m

#print coords
#print masses

for ip in range(0, n_m):
	coords[ip][0] = ip
	coords[ip][1] = ip * 2
	coords[ip][2] = ip * 3
	

print coords

moment_inertia(coords, masses)
#print masses

#nh = len(hs)
#for i in range(0, nh):
#	print hs[i].info()

#print 'found subhalos: %d' % len(hs)
#nu.particles.ellipticities(0,0)

#r = [50000, 50000, 50000]
#print ahf[0].distance(r)
#print "Halo mass: %.3f \n" % h.m
'''
