import numpy as np

from . import HaloCatalogue, Halo
from .. import snapshot, config_parser, array

class SubFindHDFSubhaloCatalogue(HaloCatalogue) :
    """
    Gadget's SubFind HDF Subhalo catalogue.

    Initialized with the parent FOF group catalogue and created
    automatically when an fof group is created
    """

    def __init__(self, group_id, group_catalogue) :
        super(SubFindHDFSubhaloCatalogue,self).__init__(group_catalogue.base)

        self._group_id = group_id
        self._group_catalogue = group_catalogue



    def __len__(self):
        if self._group_id == (len(self._group_catalogue._fof_group_first_subhalo)-1) :
            return self._group_catalogue.nsubhalos - self._group_catalogue._fof_group_first_subhalo[self._group_id]
        else:
            return (self._group_catalogue._fof_group_first_subhalo[self._group_id + 1] -
                    self._group_catalogue._fof_group_first_subhalo[self._group_id])

    def _get_halo(self, i):
        if self.base is None :
            raise RuntimeError("Parent SimSnap has been deleted")

        if i > len(self)-1 :
            raise RuntimeError("FOF group %d does not have subhalo %d"%(self._group_id, i))

        # need this to index the global offset and length arrays
        absolute_id = self._group_catalogue._fof_group_first_subhalo[self._group_id] + i

        # now form the particle IDs needed for this subhalo
        type_map = self.base._family_to_group_map

        halo_lengths = self._group_catalogue._subfind_halo_lengths
        halo_offsets = self._group_catalogue._subfind_halo_offsets

        # create the particle lists
        tot_len = 0
        for g_ptypes in list(type_map.values()) :
            for g_ptype in g_ptypes:
                tot_len += halo_lengths[g_ptype][absolute_id]

        plist = np.zeros(tot_len,dtype='int64')

        npart = 0
        for ptype in self.base._families_ordered():
            # family slice in the SubFindHDFSnap
            sl = self.base._family_slice[ptype]

            for g_ptype in type_map[ptype]:
                # add the particle indices to the particle list
                offset = halo_offsets[g_ptype][absolute_id]
                length = halo_lengths[g_ptype][absolute_id]
                ind = np.arange(sl.start + offset, sl.start + offset + length)
                plist[npart:npart+length] = ind
                npart += length

        return SubFindHDFSubHalo(i, self._group_id, self, self.base, plist)


    @property
    def base(self) :
        return self._base()


class SubFindHDFSubHalo(Halo) :
    """
    SubFind subhalo class
    """

    def __init__(self,halo_id, group_id, *args) :
        super(SubFindHDFSubHalo,self).__init__(halo_id, *args)

        self._group_id = group_id
        self._descriptor = "fof_group_%d_subhalo_%d"%(group_id,halo_id)

        # need this to index the global offset and length arrays
        absolute_id = self._halo_catalogue._group_catalogue._fof_group_first_subhalo[self._group_id] + halo_id

        # load properties
        sub_props = self._halo_catalogue._group_catalogue._sub_properties
        for key in sub_props :
            self.properties[key] = array.SimArray(sub_props[key][absolute_id], sub_props[key].units)
            self.properties[key].sim = self.base


class SubFindHDFHaloCatalogue(HaloCatalogue) :
    """
    Gadget's SubFind Halo catalogue -- used in concert with :class:`~SubFindHDFSnap`
    """

    # Names of various groups and attributes in the hdf file (which seemingly may vary in different versions of SubFind?)

    _fof_name = 'FOF'
    _subfind_name = 'SUBFIND'
    _subfind_grnr_name = 'GrNr'
    _subfind_first_gr_name = 'FirstSubOfHalo'

    _numgrps_name = 'Total_Number_of_groups'
    _numsubs_name = 'Total_Number_of_subgroups'

    _grp_offset_name = 'Offset'
    _grp_len_name = 'Length'

    _sub_offset_name = 'SUB_Offset'
    _sub_len_name = 'SUB_Length'

    def __init__(self, sim) :
        super(SubFindHDFHaloCatalogue,self).__init__(sim)

        if not isinstance(sim, snapshot.gadgethdf.SubFindHDFSnap):
            raise ValueError("SubFindHDFHaloCatalogue can only work with a SubFindHDFSnap simulation")

        self.__init_halo_offset_data()
        self.__init_subhalo_relationships()
        self.__init_halo_properties()
        self.__reshape_multidimensional_properties()
        self.__reassign_properties_from_sub_to_fof()

    def __init_ignorable_keys(self):
        self.fof_ignore = list(map(str.strip,config_parser.get("SubfindHDF","FoF-ignore").split(",")))
        self.sub_ignore = list(map(str.strip,config_parser.get("SubfindHDF","Sub-ignore").split(",")))

        for t in list(self.base._family_to_group_map.values()):
            # Don't add SubFind particles ever as this list is actually spherical overdensity
            self.sub_ignore.append(t[0])
            self.fof_ignore.append(t[0])

    def __init_halo_properties(self):
        self.__init_ignorable_keys()
        self._fof_properties = self.__get_property_dictionary_from_hdf(self._fof_name)
        self._sub_properties = self.__get_property_dictionary_from_hdf(self._subfind_name)


    def __get_property_dictionary_from_hdf(self, hdf_key):
        sim = self.base
        hdf0 = sim._hdf_files.get_file0_root()

        props = {}
        for property_key in list(hdf0[hdf_key].keys()):
            if property_key not in self.fof_ignore:
                props[property_key] = np.array([])

        for h in sim._hdf_files.iterroot():
            for property_key in list(props.keys()):
                props[property_key] = np.append(props[property_key], h[hdf_key][property_key].value)

        for property_key in list(props.keys()):
            arr_units = sim._get_units_from_hdf_attr(hdf0[hdf_key][property_key].attrs)
            if property_key in props:
                props[property_key] = props[property_key].view(array.SimArray)
                props[property_key].units = arr_units
                props[property_key].sim = sim

        return props


    def __reshape_multidimensional_properties(self):
        sub_properties = self._sub_properties
        fof_properties = self._fof_properties

        for key in list(sub_properties.keys()):
            # Test if there are no remainders, i.e. array is multiple of halo length
            # then solve for the case where this is 1, 2 or 3 dimension
            if len(sub_properties[key]) % self.nsubhalos == 0:
                ndim = len(sub_properties[key]) // self.nsubhalos
                if ndim > 1:
                    sub_properties[key] = sub_properties[key].reshape(self.nsubhalos, ndim)

            try:
                # The case fof FOF
                if len(fof_properties[key]) % self.ngroups == 0:
                    ndim = len(fof_properties[key]) // self.ngroups
                    if ndim > 1:
                        fof_properties[key] = fof_properties[key].reshape(self.ngroups, ndim)
            except KeyError:
                pass

    def __reassign_properties_from_sub_to_fof(self):
        reassign = []
        for k,v in self._sub_properties.items():
            if v.shape[0]==self.ngroups:
                reassign.append(k)

        for reassign_i in reassign:
            self._fof_properties[reassign_i] = self._sub_properties[reassign_i]
            del self._sub_properties[reassign_i]

    def __init_subhalo_relationships(self):
        nsub = 0
        nfof = 0
        self._subfind_halo_parent_groups = np.empty(self.nsubhalos, dtype=int)
        self._fof_group_first_subhalo = np.empty(self.ngroups, dtype=int)
        for h in self.base._hdf_files.iterroot():
            parent_groups = h[self._subfind_name][self._subfind_grnr_name]
            self._subfind_halo_parent_groups[nsub:nsub + len(parent_groups)] = parent_groups
            nsub += len(parent_groups)

            first_groups = h[self._subfind_name][self._subfind_first_gr_name]
            self._fof_group_first_subhalo[nfof:nfof + len(first_groups)] = first_groups
            nfof += len(first_groups)

    def __init_halo_offset_data(self):
        hdf0 = self.base._hdf_files.get_file0_root()

        self._fof_group_offsets = {}
        self._fof_group_lengths = {}
        self._subfind_halo_offsets = {}
        self._subfind_halo_lengths = {}

        self.ngroups = hdf0[self._fof_name].attrs[self._numgrps_name]
        self.nsubhalos = hdf0[self._fof_name].attrs[self._numsubs_name]

        for fam in self.base._families_ordered():
            ptypes = self.base._family_to_group_map[fam]
            for ptype in ptypes:
                self._fof_group_offsets[ptype] = np.empty(self.ngroups, dtype='int64')
                self._fof_group_lengths[ptype] = np.empty(self.ngroups, dtype='int64')
                self._subfind_halo_offsets[ptype] = np.empty(self.ngroups, dtype='int64')
                self._subfind_halo_lengths[ptype] = np.empty(self.ngroups, dtype='int64')

                curr_groups = 0
                curr_subhalos = 0

                for h in self.base._hdf_files:
                    # fof groups
                    offset = h[ptype][self._grp_offset_name]
                    length = h[ptype][self._grp_len_name]
                    self._fof_group_offsets[ptype][curr_groups:curr_groups + len(offset)] = offset
                    self._fof_group_lengths[ptype][curr_groups:curr_groups + len(offset)] = length
                    curr_groups += len(offset)

                    # subfind subhalos
                    offset = h[ptype][self._sub_offset_name]
                    length = h[ptype][self._sub_len_name]
                    self._subfind_halo_offsets[ptype][curr_subhalos:curr_subhalos + len(offset)] = offset
                    self._subfind_halo_lengths[ptype][curr_subhalos:curr_subhalos + len(offset)] = length
                    curr_subhalos += len(offset)


    def _get_halo(self, i) :
        if self.base is None :
            raise RuntimeError("Parent SimSnap has been deleted")

        if i > len(self)-1 :
            raise RuntimeError("Group %d does not exist"%i)

        type_map = self.base._family_to_group_map

        # create the particle lists
        tot_len = 0
        for g_ptypes in list(type_map.values()) :
            for g_ptype in g_ptypes:
                tot_len += self._fof_group_lengths[g_ptype][i]

        plist = np.zeros(tot_len,dtype='int64')

        npart = 0
        for ptype in self.base._families_ordered():
            # family slice in the SubFindHDFSnap
            sl = self.base._family_slice[ptype]

            for g_ptype in type_map[ptype]:
                # add the particle indices to the particle list
                offset = self._fof_group_offsets[g_ptype][i]
                length = self._fof_group_lengths[g_ptype][i]
                ind = np.arange(sl.start + offset, sl.start + offset + length)
                plist[npart:npart+length] = ind
                npart += length

        return SubFindFOFGroup(i, self, self.base, plist)


    def __len__(self) :
        return self.base._hdf_files[0].attrs[self._numgrps_name]

    @property
    def base(self):
        return self._base()


class SubFindFOFGroup(Halo) :
    """
    SubFind FOF group class
    """

    def __init__(self, group_id, *args) :
        super(SubFindFOFGroup,self).__init__(group_id, *args)

        self._subhalo_catalogue = SubFindHDFSubhaloCatalogue(group_id, self._halo_catalogue)

        self._descriptor = "fof_group_"+str(group_id)

        # load properties
        for key in list(self._halo_catalogue._fof_properties.keys()) :
            self.properties[key] = array.SimArray(self._halo_catalogue._fof_properties[key][group_id],
                                            self._halo_catalogue._fof_properties[key].units)
            self.properties[key].sim = self.base


    def __getattr__(self, name):
        if name == 'sub':
            return self._subhalo_catalogue
        else :
            return super(SubFindFOFGroup,self).__getattr__(name)