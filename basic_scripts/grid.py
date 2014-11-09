import ConfigParser
import h5py
import logging
import numpy as np


#def main():
#    hex = grid.HexPBC([ 22.34405, 18.91217,   9.91809 ,
#                         0.00000,  0.00000, -10.91895 ,
#                         0.00000, -0.00000,  -0.00000 ])
#    print hex
#
#    v1 = np.array([[22.0,   0, 9.91895]])
#    v2 = np.array([[ 0.0,   0, 9.91895]])
#    print hex.NearestPeriodic(v1, v2)
#

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


    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        print ds.shape
        print ds.attrs["dt"]

        rmsd_dT = 20 # Units of "frames"

        #rmsd_lambda = grid.RMSDLambda( b_activity = True,                 
        #                               b_scaletime = True, 
        #                               rmsd_delay = rmsd_dT * .1,         
        #                               cutoff = .38,  
        #                               sharpness = 12)
        rmsd_lambda = op.RMSDLambda( b_activity = False,                 
                                       b_scaletime = False, 
                                       rmsd_delay = rmsd_dT * .1,         
                                       cutoff = .38,  
                                       sharpness = 12)

        hex = op.HexPBC([ 22.34405, 18.91217,   9.91809 ,
                             0.00000,  0.00000, -10.91895 ,
                             0.00000, -0.00000,  -0.00000 ])

        rmsd_lambda.SetTitle()

        coord_system = coords.SlabCoords(10.0, 4.0)
        coord_system = coords.RadialCoords(10.0, 4.0)

        colorrange = None
        file_name = output_dir + "/TMV-test"

        grid.GridOP(ds, dynamic_step=rmsd_dT, op_type='rmsd', 
                colorrange=colorrange, display_type = display_type, 
                file_name=file_name, rmsd_lambda = rmsd_lambda, 
                colormap=colormap, coord_system = coord_system, 
                pbc=hex, nframes = 5)


if __name__ == "__main__":
    main()
