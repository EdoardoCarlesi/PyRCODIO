

	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	++  PyRCODIO: Python Routines for COsmology and Data I/O  ++
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		See also: http://vixra.org/abs/1704.0001

			Edoardo Carlesi 2018
			  gatto@nanowar.it
			  ecarlesi@aip.de

Tools to analyze AREPO - GADGET simulations of the Local Group and the local environment Works with GADGET1, GADGET2 and HDF5 fortran using h5py and gadget py reader.

    libSQL has to be written yet. Maybe it will not be necessary...

    libio contains some tools to identify the Lagrangian regions in the ICs and generate the mask for zoom-in initial conditions, as well as other functions to implement C and Bash functions into python

    libcosmo contains simple algorithms for the identification of clusters and LG-like objects in Constrained Simulations, as well as simple plotting routines.

    lare.py automatically identifies LG candidates in a simulation and traces its LaRe in the ICs

    satellite_stat.py identifies subhaloes in high resolution runs and computes basic properties: mass functions, anisotropy and so on

    lg_stat.py does some basic LG-related statistics

    test.py is used to test new routines as they are written

    data/ contains tables of P(k) and z(t) to be used for interpolation

TODO (ShortTerm):
	
	


TODO (LongTerm):

    Include the C libraries and bash scripts for grid generation & LagrangianRegion identification & the Fortran code for Mask generation inside a separate library within this code to make it self-consistent

    Add SQL functionalities - gather all the LG data within a database to be updated each time a new run or series of simulations runs have been finished, write some libs to access and edit it with python

	Include mpi4py and parallelize some of the routines (especially the merger tree!)



