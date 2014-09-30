RMSDAnalyze
===========

Spatially gridded RMSD analysis designed for water in TMV, but extensible to other probjects.

Requires ConfigParser configuration file formatted as:

[HDF]
file        : File to use for HDF5 storage or reading. Written to with gro2hdf5, read after conversion 
atom_dset   : Name of dataset to use with HDF5.
[plotting]
colormap    : Colormap to use with numpy. Run through pyplot.get_cmap().
