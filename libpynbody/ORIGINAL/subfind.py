import os.path
import weakref

import numpy as np

from . import HaloCatalogue, Halo
from .. import units
from ..array import SimArray


class SubfindCatalogue(HaloCatalogue):

    """
        Class to handle catalogues produced by the SubFind halo finder.
        Can import FoF groups (default) or subhalos (by setting subs=True).
        Groups are sorted by mass (descending), most massive group is halo[0]. Subhalos within each group are also sorted by mass.
        IDs within groups are sorted by binding energy, UNLESS the IDs in the snapshots are not ordered starting from 0. In that case,
        setting order=False will generate the correct halo catalogue but the ordering by binding energy within each halo will be lost.
        Additional properties calculated by SubFind can be accessed via e.g. 'halo[0].properties'.
        make_grp=True generates a particle array with the Fof group number, and -1 if the particles is not in a group. *This will take some time! v=True enables a crude progress report.
    """


    def __init__(self, sim, subs=False, order=True, make_grp=None, v=False, **kwargs):
        self._base = weakref.ref(sim)
        self._order=order
        self._subs=subs
        if self._order is False:
            if (self.base['iord'][0] != 0 or self.base['iord'][1] != 1):
                raise ValueError("IDs are not ordered. Use argument order=False to load halos for this simulation.")

        self._halos = {}
        HaloCatalogue.__init__(self,sim)
        self.dtype_int = sim['iord'].dtype
        self.dtype_flt='float32' #SUBFIND data is apparently always single precision???
        self.halodir = self._name_of_catalogue(sim)
        self.header = self._readheader()
        if subs is True:
            if self.header[6]==0:
                raise ValueError("This file does not contain subhalos")
            if make_grp is True:
                raise ValueError("subs=True and make_grp=True are not compatible.")
        self._tasks = self.header[4]
        self.ids = self._read_ids()
        self._keys={}
        self._halodat, self._subhalodat=self._read_groups()
        #self.data_len, self.data_off = self._read_groups()
        if make_grp:
            self.make_grp(v=v)

    def make_grp(self, name='grp', v=False):
        """
        #Creates a 'grp' array which labels each particle according to
        #its parent halo. This can take quite some time!
        Option: v=True prints out 'progress' in terms of total number of groups.
        #"""
        self.base[name] = self.get_group_array(v=v) #np.zeros(len(self.base), dtype=int)#self.get_group_array()

    def get_group_array(self, v=False):
        ar = np.zeros(len(self.base), dtype=int)-1
        for i in range(0, self.header[1]): #total number of groups
            if v:
                print("Halo #", i , "of", self.header[1])
            halo=self[i]
            ar[halo.get_index_list(self.base)] = halo._halo_id
        return ar

    def get_halo_properties(self, i, with_unit=True):
        if with_unit:
            extract = units.get_item_with_unit
        else:
            extract = lambda array, element: array[element]
        properties = {}
        if self._subs is False:
            for key in self._keys:
                properties[key] = extract(self._halodat[key], i)
            if self.header[6] > 0:
                properties['children'] = np.where(self._subhalodat['sub_groupNr'] == i)[
                    0]  # this is the FIRST level of substructure, sub-subhalos (etc) can be accessed via the subs=True output (below)
        else:
            for key in self._keys:
                properties[key] = extract(self._subhalodat[key], i)
            properties['children'] = np.where(self._subhalodat['sub_parent'] == i)[
                0]  # this goes down one level in the hierarchy, i.e. a subhalo will have all its sub-subhalos listed, but not its sub-sub-subhalos (those will be listed in each sub-subhalo)
        return properties

    def _get_halo(self,i):
        """This also works if the IDs are not sorted, but it will break the ordering by binding energy which is not desirable. We do however save the group's mostboundID"""
        if self._order is False:
            if self._subs is True:
                #this needs to be tested again on a snapshot that is not ordered!
                x = Halo(i, self, self.base, np.where(np.in1d(self.base['iord'], self.ids[self._subhalodat['sub_off'][i]:self._subhalodat['sub_off'][i]+self._subhalodat['sub_len'][i]]  )))
            else:
                x = Halo(i, self, self.base, np.where(np.in1d(self.base['iord'], self.ids[self._halodat['group_off'][i]:self._halodat['group_off'][i]+self._halodat['group_len'][i]]  )))

        else:
            if self._subs is False: #to use groups as halos:
                x = Halo(i, self, self.base, self.ids[self._halodat['group_off'][i]:self._halodat['group_off'][i]+self._halodat['group_len'][i]] )
            else:
                x=Halo(i, self, self.base, self.ids[self._subhalodat['sub_off'][i]:self._subhalodat['sub_off'][i]+self._subhalodat['sub_len'][i]] )

        x._descriptor = "halo_"+str(i)
        x.properties.update(self.get_halo_properties(i))
        return x

    def _readheader(self):
        header = np.array([], dtype='int32')
        filename = self.halodir + "/subhalo_tab_" + \
            self.halodir.split("_")[-1] + ".0"
        fd = open(filename, "rb")
        # read header: this is strange but it works: there is an extra value in
        # header which we delete in the next step
        header1 = np.fromfile(fd, dtype='int32', sep="", count=8)
        header = np.delete(header1, 4)
        fd.close()
        return header  # [4]

    def _read_ids(self):
        data_ids = np.array([], dtype=self.dtype_int)
        for n in range(0, self._tasks):
            filename = self.halodir + "/subhalo_ids_" + \
                self.halodir.split("_")[-1] + "." + str(n)
            fd = open(filename, "rb")
            # for some reason there is an extra value in header which we delete
            # in the next step
            header1 = np.fromfile(fd, dtype='int32', sep="", count=7)
            header = np.delete(header1, 4)
            # TODO: include a check if both headers agree (they better)
            ids = np.fromfile(fd, dtype=self.dtype_int, sep="", count=-1)
            fd.close()
            data_ids = np.append(data_ids, ids)
        return data_ids

    def _read_groups(self):
        halodat={}
        keys_flt=['mass', 'pos', 'mmean_200', 'rmean_200', 'mcrit_200', 'rcrit_200', 'mtop_200', 'rtop_200', 'contmass']
        keys_int=['group_len', 'group_off',  'first_sub', 'Nsubs', 'cont_count', 'mostboundID']
        for key in keys_flt:
            halodat[key]=np.array([], dtype=self.dtype_flt)
        for key in keys_int:
            halodat[key]=np.array([], dtype='int32')

        subhalodat={}
        subkeys_int=['sub_len', 'sub_off', 'sub_parent', 'sub_mostboundID', 'sub_groupNr']
        subkeys_flt=['sub_pos', 'sub_vel', 'sub_CM', 'sub_mass', 'sub_spin', 'sub_veldisp', 'sub_VMax', 'sub_VMaxRad', 'sub_HalfMassRad', ]
        for key in subkeys_int:
            subhalodat[key]=np.array([], dtype='int32')
        subhalodat['sub_mostboundID']=np.array([], dtype=self.dtype_int)
        #subhalodat['sub_groupNr']=np.array([], dtype=self.dtype_int) #these are special
        for key in subkeys_flt:
            subhalodat[key]=np.array([], dtype=self.dtype_flt)

        self._keys=keys_flt+keys_int
        if self._subs is True:
            self._keys=subkeys_flt+subkeys_int

        for n in range(0,self._tasks):
            filename=self.halodir+"/subhalo_tab_"+self.halodir.split("_")[-1]+"." +str(n)
            fd=open(filename, "rb")
            header1=np.fromfile(fd, dtype='int32', sep="", count=8)
            header=np.delete(header1,4)
            #read groups
            if header[0]>0:
                halodat['group_len']=np.append(halodat['group_len'], np.fromfile(fd, dtype='int32', sep="", count=header[0]))
                halodat['group_off']=np.append(halodat['group_off'], np.fromfile(fd, dtype='int32', sep="", count=header[0]))
                halodat['mass']=np.append(halodat['mass'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['pos']=np.append(halodat['pos'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=3*header[0]) )
                halodat['mmean_200']=np.append(halodat['mmean_200'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['rmean_200']=np.append(halodat['rmean_200'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['mcrit_200']=np.append(halodat['mcrit_200'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['rcrit_200']=np.append(halodat['rcrit_200'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['mtop_200']=np.append(halodat['mtop_200'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['rtop_200']=np.append(halodat['rtop_200'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['cont_count']=np.append(halodat['cont_count'], np.fromfile(fd, dtype='int32', sep="", count=header[0]))
                halodat['contmass']=np.append(halodat['contmass'],np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[0]))
                halodat['Nsubs']=np.append(halodat['Nsubs'],np.fromfile(fd, dtype='int32', sep="", count=header[0]))
                halodat['first_sub']=np.append(halodat['first_sub'],np.fromfile(fd, dtype='int32', sep="", count=header[0]))
            #read subhalos only if expected to exist from header
            if header[5]>0:
                subhalodat['sub_len']=np.append(subhalodat['sub_len'], np.fromfile(fd, dtype='int32', sep="", count=header[5]))
                subhalodat['sub_off']=np.append(subhalodat['sub_off'], np.fromfile(fd, dtype='int32', sep="", count=header[5]))
                subhalodat['sub_parent']=np.append(subhalodat['sub_parent'], np.fromfile(fd, dtype='int32', sep="", count=header[5]))
                subhalodat['sub_mass']=np.append(subhalodat['sub_mass'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[5]))
                subhalodat['sub_pos']=np.append(subhalodat['sub_pos'],np.fromfile(fd, dtype=self.dtype_flt, sep="", count=3*header[5]))
                subhalodat['sub_vel']=np.append(subhalodat['sub_vel'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=3*header[5]))
                subhalodat['sub_CM']=np.append(subhalodat['sub_CM'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=3*header[5]))
                subhalodat['sub_spin']=np.append(subhalodat['sub_spin'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=3*header[5]))
                subhalodat['sub_veldisp']=np.append(subhalodat['sub_veldisp'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[5]))
                subhalodat['sub_VMax']=np.append(subhalodat['sub_VMax'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[5]))
                subhalodat['sub_VMaxRad']=np.append(subhalodat['sub_VMaxRad'],np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[5]))
                subhalodat['sub_HalfMassRad']=np.append(subhalodat['sub_HalfMassRad'], np.fromfile(fd, dtype=self.dtype_flt, sep="", count=header[5]))
                subhalodat['sub_mostboundID']=np.append(subhalodat['sub_mostboundID'], np.fromfile(fd, dtype=self.dtype_int, sep="", count=header[5]))
                subhalodat['sub_groupNr']=np.append(subhalodat['sub_groupNr'], np.fromfile(fd, dtype='int32', sep="", count=header[5]))
            fd.close()

        halodat['pos']=np.reshape(halodat['pos'], (header[1],3))
        if header[6]>0:
            #some voodoo because some SubFind files may have (at least?) one extra entry which is not really a subhalo
            real_ones=np.where(halodat['first_sub']<header[6])[0]
            fake_ones=np.where(halodat['first_sub']>=header[6])[0]
            halodat['mostboundID']=np.zeros(len(halodat['Nsubs']),dtype=self.dtype_int)-1
            halodat['mostboundID'][real_ones]=subhalodat['sub_mostboundID'][halodat['first_sub'][real_ones]]  #useful for the case of unordered snapshot IDs

            subhalodat['sub_pos']=np.reshape(subhalodat['sub_pos'], (header[6],3))
            subhalodat['sub_vel']=np.reshape(subhalodat['sub_vel'], (header[6],3))
            subhalodat['sub_CM']=np.reshape(subhalodat['sub_CM'], (header[6],3))
            subhalodat['sub_spin']=np.reshape(subhalodat['sub_spin'], (header[6],3))

        ar_names = 'mass', 'pos', 'mmean_200', 'rmean_200', 'mcrit_200', 'rcrit_200', 'mtop_200', 'rtop_200', \
                   'sub_mass', 'sub_pos', 'sub_vel', 'sub_CM', 'sub_veldisp', 'sub_VMax', 'sub_VMaxRad', 'sub_HalfMassRad'
        ar_dimensions = 'kg', 'm', 'kg', 'm', 'kg', 'm', 'kg', 'm', \
                    'kg', 'm', 'm s^-1', 'm', 'm s^-1', 'm s^-1', 'm', 'm'

        for name, dimension in zip(ar_names, ar_dimensions):
            if name in halodat:
                halodat[name] = SimArray(halodat[name], self.base.infer_original_units(dimension))
            if name in subhalodat:
                subhalodat[name] = SimArray(subhalodat[name], self.base.infer_original_units(dimension))

        return halodat, subhalodat

    def __len__(self):
        if self._subs:
            return len(self._subhalodat['sub_pos'])
        else:
            return len(self._halodat['pos'])

    @staticmethod
    def _name_of_catalogue(sim):
        # standard path for multiple snapshot files
        snapnum = os.path.basename(
            os.path.dirname(sim.filename)).split("_")[-1]
        parent_dir = os.path.dirname(os.path.dirname(sim.filename))
        dir_path=os.path.join(parent_dir,"groups_" + snapnum)

        if os.path.exists(dir_path):
            return dir_path
        # alternative path if snapshot is single file
        else:
            snapnum = os.path.basename(sim.filename).split("_")[-1]
            parent_dir = os.path.dirname(sim.filename)
            return os.path.join(parent_dir,"groups_" + snapnum)

    @property
    def base(self):
        return self._base()

    @staticmethod
    def _can_load(sim, **kwargs):
        if os.path.exists(SubfindCatalogue._name_of_catalogue(sim)):
            return True
        else:
            return False