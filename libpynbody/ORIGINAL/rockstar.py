import glob
import os.path
import sys

import numpy as np

from . import HaloCatalogue, DummyHalo, Halo
from .. import util


class RockstarCatalogue(HaloCatalogue):
    def __init__(self, sim, dummy=False, pathname=None, format_revision=None,
                 filenames=None, sort=False, **kwargs):
        """Initialize a RockstarCatalogue.

        **kwargs** :


        *dummy*: if True, the particle file is not loaded, and all
                 halos returned are just dummies (with the correct
                 properties dictionary loaded). Use load_copy to get
                 the actual data in this case.

        *sort*: if True, resort the halos into descending order of
                particle number. Otherwise, leave in RockStar output order.

        *filenames*: a list of filenames of each of the RockStar outputs.
                     You probably want to use pathname instead, which specifies
                     the path to the output folder.

        *pathname*: the path of the output folder with the individual RockStar outputs

        *format_revision*: Override the header's format revision information. Specify
                    1, 2, or 'caterpillar' for Rockstar prior to 2014, post 2014, and
                    customized for the caterpillar project respectively

        """
        HaloCatalogue.__init__(self, sim)

        self._dummy = dummy

        if filenames is not None:
            self._files = filenames
        else:
            if pathname is None:
                pathname = os.path.dirname(sim.filename)
            self._files = glob.glob(os.path.join(pathname,'halos*.bin'))
            if len(self._files)==0 :
                self._files = glob.glob(os.path.join(pathname, 'halos*.boundbin'))
            self._files.sort()

        if len(self._files)==0:
            raise IOError("Could not find any Rockstar output. Try specifying pathname='/path/to/rockstar/outputfolder'")

        self._cpus = [RockstarCatalogueOneCpu(sim,dummy,file_i, format_revision=format_revision) for file_i in self._files]
        self._prune_files_from_wrong_scalefactor()
        self._cpus[0]._init_iord_to_fpos()
        for cpu in self._cpus:
            cpu._iord_to_fpos = self._cpus[0]._iord_to_fpos

        self._init_index_ar()
        if sort:
            self._sort_index_ar()

    def _prune_files_from_wrong_scalefactor(self):
        new_cpus = []
        new_files = []
        for file,cpu in zip(self._files,self._cpus):
            if abs(self.base.properties['a']-cpu._head['scale'][0])<1e-6:
                new_cpus.append(cpu)
                new_files.append(file)
        self._cpus = new_cpus
        self._files = new_files

    def _pass_on(self, function, *args, **kwargs):
        if self._index_ar is not None:
            cpui, hi = self._cpus[self._index_ar[args[0],0]], self._index_ar[args[0],1]
            return function(cpui,hi,*args[1:],**kwargs)
        for i in self._cpus:
            try:
                return function(i,*args,**kwargs)
            except KeyError:
                pass

    def __getitem__(self, k):
        return self._pass_on(RockstarCatalogueOneCpu.__getitem__, k)

    def load_copy(self, k):
        return self._pass_on(RockstarCatalogueOneCpu.load_copy, k)

    def get_group_array(self):
        ar = np.zeros(len(self.base), dtype=int)-1
        for cpu_i in self._cpus:
            cpu_i._update_grp_array(ar)
        return ar

    def make_grp(self, name='grp'):
        """
        Creates a 'grp' array which labels each particle according to
        its parent halo.
        """
        self.base[name]= self.get_group_array()

    def __len__(self):
        return sum([len(x) for x in self._cpus])

    def _init_index_ar(self):
        index_ar = np.empty((len(self),2),dtype=np.int32)

        for cpu_id, cpu in enumerate(self._cpus):
            i0 = cpu._halo_min
            i1 = cpu._halo_max
            index_ar[i0:i1,0]=cpu_id
            index_ar[i0:i1,1]=np.arange(i0,i1,dtype=np.int32)

        self._index_ar = index_ar


    def _sort_index_ar(self):
        num_ar = np.empty(len(self))

        for cpu_id, cpu in enumerate(self._cpus):
            i0 = cpu._halo_min
            i1 = cpu._halo_max
            num_ar[i0:i1]=self._cpus[cpu_id]._halo_lens

        num_ar = np.argsort(num_ar)[::-1]
        self._index_ar = self._index_ar[num_ar]


    @staticmethod
    def _can_run(sim):
        return False

    @staticmethod
    def _can_load(sim, **kwargs):
        for file in glob.glob(os.path.join(os.path.dirname(sim.filename), 'halos*.bin')):
            if os.path.exists(file):
                return True
        return False


class RockstarCatalogueOneCpu(HaloCatalogue):
    """
    Class to handle sub-catalogues produced by Rockstar. Users should normally not use this class,
    rather using RockstarCatalogue which collates the multiple sub-files that Rockstar produces.
    """

    head_type = np.dtype([('magic',np.uint64),('snap',np.int64),
                          ('chunk',np.int64),('scale','f'),
                          ('Om','f'),('Ol','f'),('h0','f'),
                          ('bounds','f',6),('num_halos',np.int64),
                          ('num_particles',np.int64),('box_size','f'),
                          ('particle_mass','f'),('particle_type',np.int64),
                          ('format_revision',np.int32),
                          ('rockstar_version',np.str_,12)])

    halo_types = {
        1: np.dtype([('id', np.int64), ('pos', 'f', 3), ('vel', 'f', 3),
                  ('corevel', 'f', 3), ('bulkvel', 'f', 3), ('m', 'f'),
                  ('r', 'f'),
                  ('child_r', 'f'), ('vmax_r', 'f'), ('mgrav', 'f'),
                  ('vmax', 'f'), ('rvmax', 'f'), ('rs', 'f'),
                  ('klypin_rs', 'f'), ('vrms', 'f'), ('J', 'f', 3),
                  ('energy', 'f'), ('spin', 'f'), ('alt_m', 'f', 4),
                  ('Xoff', 'f'), ('Voff', 'f'), ('b_to_a', 'f'),
                  ('c_to_a', 'f'), ('A', 'f', 3), ('b_to_a2', 'f'),
                  ('c_to_a2', 'f'), ('A2', 'f', 3), ('bullock_spin', 'f'),
                  ('kin_to_pot', 'f'), ('m_pe_b', 'f'), ('m_pe_d', 'f'),
                  ('num_p', np.int64), ('num_child_particles', np.int64),
                  ('p_start', np.int64), ('desc', np.int64),
                  ('flags', np.int64), ('n_core', np.int64),
                  ('min_pos_err', 'f'), ('min_vel_err', 'f'),
                  ('min_bulkvel_err', 'f')], align=True)  # Rockstar format v1
        ,
        2: np.dtype([('id',np.int64),('pos','f',3),('vel','f',3),
                          ('corevel','f',3),('bulkvel','f',3),('m','f'),
                          ('r','f'),
                          ('child_r','f'),('vmax_r','f'),('mgrav','f'),
                          ('vmax','f'),('rvmax','f'),('rs','f'),
                          ('klypin_rs','f'),('vrms','f'),('J','f',3),
                          ('energy','f'),('spin','f'),('alt_m','f',4),
                          ('Xoff','f'),('Voff','f'),('b_to_a','f'),
                          ('c_to_a','f'),('A','f',3),('b_to_a2','f'),
                          ('c_to_a2','f'),('A2','f',3),('bullock_spin','f'),
                          ('kin_to_pot','f'),('m_pe_b','f'),('m_pe_d','f'),
                          ('halfmass_radius','f'),
                          ('num_p',np.int64),('num_child_particles',np.int64),
                          ('p_start',np.int64),('desc',np.int64),
                          ('flags',np.int64),('n_core',np.int64),
                          ('min_pos_err','f'),('min_vel_err','f'),
                          ('min_bulkvel_err','f')], align=True), # Rockstar format v2, includes halfmass_radius

        'caterpillar': np.dtype([('id',np.int64),
                                 ('pos','f',3),('vel','f',3),
                          ('corevel','f',3),('bulkvel','f',3),('m','f'),
                          ('r','f'),
                          ('child_r','f'),('vmax_r','f'),('mgrav','f'),
                          ('vmax','f'),('rvmax','f'),('rs','f'),
                          ('klypin_rs','f'),('vrms','f'),('J','f',3),
                          ('energy','f'),('spin','f'),('alt_m','f',4),
                          ('Xoff','f'),('Voff','f'),('b_to_a','f'),
                          ('c_to_a','f'),('A','f',3),('b_to_a2','f'),
                          ('c_to_a2','f'),('A2','f',3),('bullock_spin','f'),
                          ('kin_to_pot','f'),('m_pe_b','f'),('m_pe_d','f'),
                          ('halfmass_radius','f'),
                          ('num_p',np.int64),('num_child_particles',np.int64),

                          ('p_start',np.int64),('desc',np.int64),
                          ('flags',np.int64),('n_core',np.int64),
                          ('min_pos_err','f'),('min_vel_err','f'),
                          ('min_bulkvel_err','f'),
                          ('num_bound', 'i8'), ('num_iter', 'i8')]
                                , align=True) # Hacked rockstar from caterpillar project
    }


    def __init__(self, sim, dummy=False, filename=None, format_revision=None, **kwargs):
        """Initialize a RockstarCatalogue.

        **kwargs** :


        *dummy*: if True, the particle file is not loaded, and all
                 halos returned are just dummies (with the correct
                 properties dictionary loaded). Use load_copy to get
                 the actual data in this case.


        """

        HaloCatalogue.__init__(self,sim)

        self._dummy = dummy

        if filename is not None: self._rsFilename = filename
        else:
            self._rsFilename = util.cutgz(glob.glob('halos*.bin')[0])

        try :
            f = util.open_(self._rsFilename, 'rb')
        except IOError:
            raise IOError("Halo catalogue not found -- check the file name of catalogue data or try specifying a catalogue using the filename keyword")

        with f:
            self._head = np.fromstring(f.read(self.head_type.itemsize),
                                       dtype=self.head_type)
            if format_revision is None:
                format_revision = self._head['format_revision'][0]

            self.halo_type = self.halo_types[format_revision]
            unused = f.read(256 - self._head.itemsize)

            self._nhalos = self._head['num_halos'][0]

            self._load_rs_halos(f,sim)




    def __len__(self):
        return len(self._halo_lens)

    def _load_all(self):
        for i in range(self._halo_min, self._halo_max):
            self._load_rs_halo_if_required(i)
            self._load_rs_particles_for_halo_if_required(i)

    def calc_item(self, i):
        if self.base is None:
            raise RuntimeError("Parent SimSnap has been deleted")

        self._load_rs_halo_if_required(i)
        self._load_rs_particles_for_halo_if_required(i)

        return self._halos[i]

    def load_copy(self, i):
        """Load a fresh SimSnap with only the particles in halo i"""
        if i<self._halo_min or i>=self._halo_max:
            raise KeyError("No such halo")

        from . import load
        return load(self.base.filename, take=self._get_particles_for_halo(i))

    def _update_grp_array(self, ar):
        """Insert into an existing grp array halos from this CPU"""
        self._load_all()
        for halo in list(self._halos.values()):
            ar[halo.get_index_list(self.base)] = halo._halo_id

    def _setup_children(self):
        """
        Creates a 'children' array inside each halo's 'properties'
        listing the halo IDs of its children. Used in case the reading
        of substructure data from the AHF-supplied _substructure file
        fails for some reason.
        """

        for i in range(self._nhalos):
            self._halos[i+1].properties['children'] = []

        for i in range(self._nhalos):
            host = self._halos[i+1].properties.get('hostHalo', -2)
            if host > -1:
                try:
                    self._halos[host+1].properties['children'].append(i+1)
                except KeyError:
                    pass





    def _load_rs_halo_if_required(self, i):
        if i not in self._halos:
            self._halos[i] = self._get_dummy_for_halo(i)

    def _get_dummy_for_halo(self, n):
        if n<self._halo_min or n>=self._halo_max:
            raise KeyError("No such halo")

        with util.open_(self._rsFilename, 'rb') as f:
            f.seek(self._haloprops_offset+(n-self._halo_min)*self.halo_type.itemsize)
            halo_data = np.fromfile(f, dtype=self.halo_type, count=1)

        hn = DummyHalo()
        # TODO: properties are in Msun / h, Mpc / h
        hn.properties = dict(list(zip(halo_data.dtype.names,halo_data[0])))
        return hn


    def _load_rs_halos(self, f, sim):
        self._haloprops_offset = f.tell()
        self._halo_offsets = np.empty(self._head['num_halos'][0],dtype=np.int64)
        self._halo_lens = np.empty(self._head['num_halos'][0],dtype=np.int64)
        offset = self._haloprops_offset+self.halo_type.itemsize*self._head['num_halos'][0]


        self._halo_min = int(np.fromfile(f, dtype=self.halo_type, count=1)['id'])
        self._halo_max = int(self._halo_min+self._head['num_halos'][0])

        f.seek(self._haloprops_offset)

        this_id = self._halo_min

        for h in range(self._head['num_halos'][0]):
            halo_data =np.fromfile(f, dtype=self.halo_type, count=1)
            assert halo_data['id']==this_id
            self._halo_offsets[this_id-self._halo_min] = int(offset)
            if 'num_bound' in self.halo_type.names:
                num_ptcls = int(halo_data['num_bound'])
            else:
                num_ptcls = int(halo_data['num_p'])
            self._halo_lens[this_id-self._halo_min] = num_ptcls
            offset+=num_ptcls*np.dtype('int64').itemsize
            this_id+=1


    def _get_particles_for_halo(self, num):
        self._init_iord_to_fpos()
        with util.open_(self._rsFilename, 'rb') as f:
            f.seek(self._halo_offsets[num-self._halo_min])
            halo_ptcls=np.fromfile(f,dtype=np.int64,count=self._halo_lens[num-self._halo_min])
            halo_ptcls = self._iord_to_fpos[halo_ptcls]
            halo_ptcls.sort()

        return halo_ptcls

    def _load_rs_particles_for_halo(self, num):
        halo_ptcls = self._get_particles_for_halo(num)

        properties_from_proxy = self._halos[num].properties

        self._halos[num]=Halo(num, self, self.base, halo_ptcls)
        self._halos[num]._descriptor = "halo_"+str(num)

        self._halos[num].properties.update(properties_from_proxy)

    def _load_rs_particles_for_halo_if_required(self, num):
        if isinstance(self._halos[num], DummyHalo) and not self._dummy:
            self._load_rs_particles_for_halo(num)



    def _load_ahf_substructure(self, filename):
        f = util.open_(filename)

        for i in range(len(self._halos)):

            haloid, nsubhalos = [int(x) for x in f.readline().split()]
            self._halos[haloid+1].properties['children'] = [
                int(x)+1 for x in f.readline().split()]

        f.close()




    def writestat(self, outfile=None, hubble=None):
        """
        write a condensed skid.stat style ascii file from ahf_halos
        file.  header + 1 halo per line. should reproduce `Alyson's
        idl script' except does not do last 2 columns (Is it a
        satellite?) and (Is central halo is `false'ly split?).  output
        units are set to Mpc Msun, km/s.

        user can specify own hubble constant hubble=(H0/(100
        km/s/Mpc)), ignoring the snaphot arg for hubble constant
        (which sometimes has a large roundoff error).
        """
        s = self._base()
        mindarkmass = min(s.dark['mass'])

        if hubble is None:
            hubble = s.properties['h']

        if outfile is None: outfile = self._base().filename+'.stat'
        print("write stat file to ", outfile)
        fpout = open(outfile, "w")
        header = "#Grp  N_tot     N_gas      N_star    N_dark    Mvir(M_sol)       Rvir(kpc)       GasMass(M_sol) StarMass(M_sol)  DarkMass(M_sol)  V_max  R@V_max  VelDisp    Xc   Yc   Zc   VXc   VYc   VZc   Contam   Satellite?   False?   ID_A"
        print(header, file=fpout)
        for ii in np.arange(self._nhalos)+1:
            print('%d '%ii, end=' ')
            sys.stdout.flush()
            h = self[ii].properties  # halo index starts with 1 not 0
##  'Contaminated'? means multiple dark matter particle masses in halo)"
            icontam = np.where(self[ii].dark['mass'] > mindarkmass)
            if (len(icontam[0]) > 0):
                contam = "contam"
            else:
                contam = "clean"
## may want to add implement satellite test and false central breakup test.
            ss = "     "  # can adjust column spacing
            outstring = str(ii)+ss
            outstring += str(len(self[ii]))+ss+str(len(self[ii].g))+ss
            outstring += str(len(self[ii].s)) + ss+str(len(self[ii].dark))+ss
            outstring += str(h['m']/hubble)+ss+str(h['r']/hubble)+ss
            outstring += str(self[ii].g['mass'].in_units('Msol').sum())+ss
            outstring += str(self[ii].s['mass'].in_units('Msol').sum())+ss
            outstring += str(self[ii].d['mass'].in_units('Msol').sum())+ss
            outstring += str(h['vmax'])+ss+str(h['vmax_r']/hubble)+ss
            outstring += str(h['vrms'])+ss
        ## pos: convert kpc/h to mpc (no h).
            outstring += str(h['pos'][0][0]/hubble)+ss
            outstring += str(h['pos'][0][1]/hubble)+ss
            outstring += str(h['pos'][0][2]/hubble)+ss
            outstring += str(h['vel'][0][0])+ss+str(h['vel'][0][1])+ss
            outstring += str(h['vel'][0][2])+ss
            outstring += contam+ss
            outstring += "unknown" + \
                ss  # unknown means sat. test not implemented.
            outstring += "unknown"+ss  # false central breakup.
            print(outstring, file=fpout)
        fpout.close()

    def writetipsy(self, outfile=None, hubble=None):
        """
        write halos to tipsy file (write as stars) from ahf_halos
        file.  returns a shapshot where each halo is a star particle.

        user can specify own hubble constant hubble=(H0/(100
        km/s/Mpc)), ignoring the snaphot arg for hubble constant
        (which sometimes has a large roundoff error).
        """
        from . import tipsy
        from .analysis import cosmology
        from snapshot import _new as new
        import math
        s = self._base()
        if outfile is None: outfile = s.filename+'.gtp'
        print("write tipsy file to ", outfile)
        sout = new(star=self._nhalos)  # create new tipsy snapshot written as halos.
        sout.properties['a'] = s.properties['a']
        sout.properties['z'] = s.properties['z']
        sout.properties['boxsize'] = s.properties['boxsize']
        if hubble is None: hubble = s.properties['h']
        sout.properties['h'] = hubble
    ### ! dangerous -- rho_crit function and unit conversions needs simplifying
        rhocrithhco = cosmology.rho_crit(s, z=0, unit="Msol Mpc^-3 h^2")
        lboxkpc = sout.properties['boxsize'].ratio("kpc a")
        lboxkpch = lboxkpc*sout.properties['h']
        lboxmpch = lboxkpc*sout.properties['h']/1000.
        tipsyvunitkms = lboxmpch * 100. / (math.pi * 8./3.)**.5
        tipsymunitmsun = rhocrithhco * lboxmpch**3 / sout.properties['h']

        print("transforming ", self._nhalos, " halos into tipsy star particles")
        for ii in range(self._nhalos):
            h = self[ii+1].properties
            sout.star[ii]['mass'] = h['m']/hubble / tipsymunitmsun
            ## tipsy units: box centered at 0. (assume 0<=x<=1)
            sout.star[ii]['x'] = h['pos'][0][0]/lboxmpch - 0.5
            sout.star[ii]['y'] = h['pos'][0][1]/lboxmpch - 0.5
            sout.star[ii]['z'] = h['pos'][0][2]/lboxmpch - 0.5
            sout.star[ii]['vx'] = h['vel'][0][0]/tipsyvunitkms
            sout.star[ii]['vy'] = h['vel'][0][1]/tipsyvunitkms
            sout.star[ii]['vz'] = h['vel'][0][2]/tipsyvunitkms
            sout.star[ii]['eps'] = h['r']/lboxkpch
            sout.star[ii]['metals'] = 0.
            sout.star[ii]['phi'] = 0.
            sout.star[ii]['tform'] = 0.
        print("writing tipsy outfile %s"%outfile)
        sout.write(fmt=tipsy.TipsySnap, filename=outfile)
        return sout

    def writehalos(self, hubble=None, outfile=None):
        s = self._base()
        if outfile is None:
            statoutfile = s.filename+".rockstar.stat"
            tipsyoutfile = s.filename+".rockstar.gtp"
        else:
            statoutfile = outfile+'.stat'
            gtpoutfile = outfile+'.gtp'
        self.make_grp()
        self.writestat(statoutfile, hubble=hubble)
        shalos = self.writetipsy(gtpoutfile, hubble=hubble)
        return shalos