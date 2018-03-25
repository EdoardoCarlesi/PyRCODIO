
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	++  PyRCODIO: Python Routines for COsmology and Data I/O  ++
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			Edoardo Carlesi 2018
			  gatto@nanowar.it
			  ecarlesi@aip.de

Tools to analyze AREPO - GADGET simulations of the Local Group and the local environment
Works with GADGET1, GADGET2 and HDF5 fortran using h5py and gadget py reader.

- libics contains some tools to identify the Lagrangian regions in the ICs and generate the mask
for zoom-in initial conditions.

- libcosmo contains simple algorithms for the identification of clusters and LG-like objects in 
Constrained Simulations, as well as simple plotting routines.

- lare.py automatically identifies LG candidates in a simulation and traces its LaRe in the ICs

- satellite_stat.py identifies subhaloes in high resolution runs and computes basic properties:
mass functions, anisotropy and so on

- lg_stat.py does some basic LG-related statistics

- data/ contains tables of P(k) and z(t) to be used for interpolation

TODO: 

* Include the C libraries and bash scripts for grid generation & LagrangianRegion identification & the Fortran code for Mask generation 
inside a separate library within this code to make it self-consistent

* Add SQL functionalities - gather all the LG data within a database to be updated each time a new run or series of simulations runs have
been finished, write some libs to access and edit it with python

