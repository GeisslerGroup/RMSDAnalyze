import ConfigParser
import h5py
import logging
import numpy as np


def main():
    logging.basicConfig(level=logging.DEBUG)
    config = ConfigParser.RawConfigParser()
    config.read('RMSDAnalyze.cfg')
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    display_type = config.get('plotting','display_type')
    output_dir  = config.get('plotting','output_dir')

    import matplotlib
    if display_type != 'display':
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import RMSDAnalyze.grid as grid
    import RMSDAnalyze.coords as coords
    import RMSDAnalyze.op as op
    
    colormap= config.get('plotting','colormap')
    colormap= plt.cm.get_cmap(colormap)
    
    nframes = 100

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        print ds.shape
        print ds.attrs["dt"]
        rmsd_dT = 20 # Units of "frames"
        rmsd_lambda = op.RMSDLambda( b_activity = False,                 
                                       b_scaletime = False, 
                                       rmsd_delay = rmsd_dT * .1,         
                                       cutoff = .38,  
                                       sharpness = 12)
        hex = op.HexPBC([ 22.34405, 18.91217,   9.91809 ,
                             0.00000,  0.00000, -10.91895 ,
                             0.00000, -0.00000,  -0.00000 ])
        center = hex.GetPBCCenter()
        rmsd_lambda.SetTitle()
        coord_system = coords.SlabCoords(10.0, 4.0)
        #coord_system = coords.RadialCoords(10.0, 4.0)
        colorrange = None
        file_name = output_dir + "/TMV-test"
        stars = [[0, 0], 
                [3.5, -1], 
                [2.5, -1.5],
                [9.5, 3.5]]
        stars_colors = ['k', 'b', 'r', '#984ea3']

        data, density, pos = grid.GridOP(
                ds, dynamic_step=rmsd_dT, op_type='rmsd', 
                colorrange=colorrange, display_type = display_type, 
                file_name=file_name, rmsd_lambda = rmsd_lambda, 
                colormap=colormap, coord_system = coord_system, 
                pbc=hex, nframes = nframes, center = center, plot = False)
        print "Going to OP Plotter"
        atom_to_mark = []
        for monomer in xrange(34):
            for resatm in xrange(10):
                cys101 = 1467 + 2461 * monomer + resatm
                atom_to_mark.append(cys101)
            for resatm in xrange(10):
                cys157 = 2388 + 2461 * monomer + resatm
                atom_to_mark.append(cys157)
        logging.debug("Length of mark atoms: {}".format(len(atom_to_mark)))

        grid.GridOPPlotter(ds, data, density, pos,
                dynamic_step=rmsd_dT, op_type='rmsd', 
                colorrange=colorrange, display_type = display_type, 
                file_name=file_name, rmsd_lambda = rmsd_lambda, 
                colormap=colormap, coord_system = coord_system, 
                pbc=hex, nframes = nframes, center = center,
                stars = stars, stars_colors=stars_colors,
                atom_to_mark = atom_to_mark)

if __name__ == "__main__":
    main()
