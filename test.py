'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    test.py: this file is used to test new routines and functionalities
'''

import read_files as rf
import numpy as np
import pickle

run = '00_06'
r_iso = 2e+3
r_max = 1.2e+3
r_min = 0.4e+3
ratio_max = 4.0
vrad_max = 10.0
radius = 8.0e+3
m_min = 0.5e+12
m_max = 2.5e+12
d_max = 10.e+3

#lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
#lg_model = LocalGroupModel(r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)

center = [50.e+3, 50.e+3, 50.e+3]


f_web = '/home/edoardo/CLUES/DATA/ICs/ic_web_00_00.000064.Vweb-csv'

web = rf.extract_vweb(file_name=f_web, center=center, radius=radius)


'''
ahf_center = find_halos_point(center, ahf, radius)
hs = find_lg(ahf_center, lg_model)

print(hs)


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


'''
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

(lg_model, lg_dict) = lg_models()


'/z/carlesi/STORE/512box100/AHF/00_06/merged_054.AHF_halos'
print('# Run\tM(1.e+14)\td(True)\t\tSGx,\t\tSGy,\t\tSGz')
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
=======
    import read_files as rf
import tools as t
import halo_utils as hu
import pandas as pd
import units as u

print(u.particle_density(1024, 100.0))


file_ahf = '/media/edoardo/data1/DATA/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
file_mah = '/media/edoardo/data1/DATA/00_06/00/halo_'
form_mah = '.allinfo'
data_path = '/home/edoardo/CLUES/PyRCODIO/data/'

time = rf.read_time(data_path)
all_halo_mah = []

halos = rf.read_ahf_halo(file_ahf)

c = [50.e+3, 50.e+3, 50.e+3]
r = 2.e+3

id_list = hu.halo_ids_around_center(halos, c, r)
mahs = rf.read_mah_halo(id_list, file_mah, time)
#print(mahs)
#mahs[1].host = mahs[0]
#print(mahs[1].trajectory_around_host())
#print(mahs[0].m_t())
#print(mahs[10].t_m_max())
#print(mahs[0].x_t())
#print(mahs[0].formation_time())
#print(mahs[0].last_major_merger())
#print(mahs[0].m_t())
#print(mahs[0].last_major_merger())

#print(mahs[0].head())
#print(mahs[0]['ID'])
#cols = str(mahs[0].columns)
#row = str(mahs[0].iloc[[1]])
#print(row.split())
#print(cols.split())
#mahs[0].columns = cols.split()
#print(mahs[0].head()) 
#print(mahs[0]['#ID(1)']) 
#print(mahs[0]['Mvir(4)']) 
#print(mahs[0]['HostHalo(2)']) 

"""
    TEST READ IN ROUTINES FOR MAH and AHF

file_halo = '/media/edoardo/Elements/CLUES/DATA/2048/00_06/00/snapshot_054.0000.z0.000.AHF_halos'
file_tree = '/media/edoardo/Elements/CLUES/DATA/trees/2048/00_06/00/df_all_ids.csv'
halo = rf.read_ahf_halo(file_halo)
tree = rf.read_csv_tree(file_tree)

id0 = halo['ID'].loc[0]
hh = halo[halo['ID'] == id0]
this_halo = hu.Halo(hh)

this_halo.assign_subhalos(halo)

hs = hu.HaloHistory(10)

hs.halos.append(hh)
hs.trajectory_around_host()

print(tree[tree])

#print(this_halo.info())
#print(this_halo.distance(pos))
'''



