RMSDAnalyze
===========

RMSD analysis designed for water in TMV, but extensible to other probjects.

1. Example scripts
interface.py: 
  a. Create an HDF5 file from a .gro "movie" file
  b. Run a grid analysis on an HDF5 file
  c. Compute the atom-wise RMSD for an atomistic visualizer


2. Configuration
Requires ConfigParser (expects  configuration file formatted as:

[HDF]
file        : File to use for HDF5 storage or reading. Written to with gro2hdf5, read after conversion 
atom_dset   : Name of dataset to use with HDF5.
[plotting]
colormap    : Colormap to use with numpy. Run through pyplot.get_cmap().



RMSDAnalyze -- The library
==========================
RMSDAnalyze.convert -- Conversion from .gro to .hdf5
RMSDAnalyze.grid    -- Calculations on a grid
