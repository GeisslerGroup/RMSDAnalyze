import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import RMSDAnalyze.local as local
import RMSDAnalyze.grid as grid
import ConfigParser
import h5py
import numpy as np

if __name__=="__main__":
    config = ConfigParser.RawConfigParser()
    config.read('.RMSDAnalyze')
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    plt.hold(True)

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]

        my_coords = [#local.LocalSlabCoords(2.0, 4.0, -2.5, -.5, 2.0),
                local.LocalSlabCoords(9., 10., 3., 4., 1.),
                local.LocalSlabCoords(-.5, .5, -.5, .5, 1.)]
        cs_legend = [#'chromophore', 
                'external', 'central pore']
        plots = []
        for i,coords in enumerate(my_coords):

            op_bins = np.linspace(0,1,1000)
            op_dist_t = []
            op_legend = []
            times = xrange(5,51,5)
            for i in times:
                dt = "{}ps".format(i*.1)
                print("Using dt = {}".format(dt))
                rmsd_lambda = grid.RMSDLambda( b_activity = False,                 
                                               b_scaletime = False, 
                                               rmsd_delay = i)
                #                        cutoff=args.rmsd_activityparm[0],  
                #                        sharpness=args.rmsd_activityparm[1])
                op_legend.append(dt)
                op_dist = local.LocalOP(ds, dynamic_step = i, op_type='rmsd',
                        rmsd_lambda = rmsd_lambda, water_pos=83674, ion_pos = 423427,
                        coord_system=coords,
                        bins = op_bins)
                op_dist /= float(np.sum(op_dist))
                op_dist_t.append(op_dist)
    
            op_bins_ctr = (op_bins[0:-1] + op_bins[1:]) / 2
            op_t = [np.sum(op_bins_ctr * op_dist) for op_dist in op_dist_t]
            
            print op_t
            plot, = plt.plot(times, op_t)
            plots.append(plot)

#        legend = []
#        for op_dist in op_dist_t:
#            fig, = plt.plot((op_bins[0:-1] + op_bins[1:]) / 2, op_dist)
#            legend.append(fig)
        plt.legend(plots, cs_legend)
        plt.savefig("/home/jhaberstroh/Dropbox/Physics/subgroup/2014-10-28/TMV/cool_file2.png")
