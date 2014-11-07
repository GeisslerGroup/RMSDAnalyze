RMSDAnalyze
===========

To use without installation, add the current directory to your PYTHONPATH

RMSD analysis designed for water in TMV, but extensible to other probjects.



RMSDAnalyze -- The library
==========================
RMSDAnalyze.convert -- Conversion from .gro to .hdf5
RMSDAnalyze.grid    -- En-masse calculations on a grid
RMSDAnalyze.local   -- The in-depth local RMSD analysis tool



basic_scripts: Example calculations
==========================
    gro2hdf.py -- run example code from RMSDAnalyze.convert
    grid.py  -- run example code from RMSDAnalyze.grid
    local.py  -- run example code from RMSDAnalyze.local


Scripts are setup to run with certain parameters, but many can be edited in 
the configuration. You should poke around in these scripts and see what 
content is demonstrated inside of them! 

All of these example scripts are written to run with my TMV files; you will
need to edit, by hand, parameters like the box-vectors, indexes of water
molecules and ion molecules, and the coordinate system. However, this are
all certainly things that you NEED to edit anyway, so the open API provides
objects to handle many cases and is extensible for custom objects as needed.

The ConfigParser configuration:
    [HDF]
    file        : File to use for HDF5 storage or reading. Written to with gro2hdf5, read after conversion 
    atom_dset   : Name of dataset to use with HDF5.
    [plotting]
    colormap    : Colormap to use with numpy. Run through pyplot.get_cmap().
