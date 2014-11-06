import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ConfigParser
import matplotlib.pyplot as plt
import h5py
import RMSDAnalyze.grid as grid
import logging


def main():
    logging.basicConfig(level=logging.DEBUG)
    config = ConfigParser.RawConfigParser()
    config.read('RMSDAnalyze.cfg')
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    colormap= config.get('plotting','colormap')
    colormap= plt.cm.get_cmap(colormap)


    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        print ds.shape
        print ds.attrs["dt"]

        rmsd_dT = 199 # Units of "frames"

        #rmsd_lambda = grid.RMSDLambda( b_activity = True,                 
        #                               b_scaletime = True, 
        #                               rmsd_delay = rmsd_dT * .1,         
        #                               cutoff = .38,  
        #                               sharpness = 12)
        rmsd_lambda = grid.RMSDLambda( b_activity = False,                 
                                       b_scaletime = False, 
                                       rmsd_delay = rmsd_dT * .1,         
                                       cutoff = .38,  
                                       sharpness = 12)
        rmsd_lambda.SetTitle()

        coord_system = grid.SlabCoords(10.0, 4.0)

        colorrange = [0, .5]
        display_type = 'display'
        file_name ="/home/jhaberstroh/Dropbox/Physics/subgroup/2014-10-28/TMV/real_cool_file"

        grid.GridOP(ds, dynamic_step=rmsd_dT, op_type='rmsd', 
                colorrange=colorrange, display_type = display_type, 
                file_name=file_name, rmsd_lambda = rmsd_lambda, 
                colormap=colormap, coord_system = coord_system)


if __name__ == "__main__":
    main()
