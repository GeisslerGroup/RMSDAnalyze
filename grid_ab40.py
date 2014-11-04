import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ConfigParser
import matplotlib.pyplot as plt
import h5py
import RMSDAnalyze.grid as grid
import logging

def main():
    config = ConfigParser.RawConfigParser()
    config.read('.RMSDAnalyze')
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    colormap= config.get('plotting','colormap')
    colormap= plt.cm.get_cmap(colormap)


    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        print ds.shape
        print ds.attrs["dt"]

        rmsd_dT = 2 # Units of "frames"

        rmsd_lambda = grid.RMSDLambda( b_activity = False,                 
                                       b_scaletime = False, 
                                       rmsd_delay = rmsd_dT * .1,         
                                       cutoff = .38,  
                                       sharpness = 12)
        rmsd_lambda.SetTitle()


        coord_system = grid.SlabCoords(6.0, 2.0, dir1 = 2, dir2 = 1, thickness=10)
        gridsize = [30,10]

        colorrange = [0, .4]
        display_type = 'png'
        file_name ="/home/jhaberstroh/Dropbox/Physics/subgroup/2014-10-28/ab_image"

        logging.basicConfig(level=logging.DEBUG)
    
        grid.GridOP(ds, dynamic_step=rmsd_dT, op_type='rmsd', 
                colorrange=colorrange, display_type = display_type, 
                file_name=file_name, rmsd_lambda = rmsd_lambda, 
                colormap=colormap, coord_system = coord_system,
                water_pos=5761, ion_pos=58336, gridsize=gridsize)


if __name__ == "__main__":
    main()