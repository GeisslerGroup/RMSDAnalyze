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

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]


        op_bins = np.linspace(0,1,1000)
        op_dist_t = []
        op_legend = []
        for i in xrange(5,50,5):
            dt = "{}ps".format(i*.1)
            print("Using dt = {}".format(dt))
            rmsd_lambda = grid.RMSDLambda( b_activity = False,                 
                                           b_scaletime = True, 
                                           rmsd_delay = i)
            #                        cutoff=args.rmsd_activityparm[0],  
            #                        sharpness=args.rmsd_activityparm[1])
            op_legend.append(dt)
            op_dist = local.LocalOP(ds, dynamic_step = i, op_type='rmsd',
                    rmsd_lambda = rmsd_lambda, water_pos=83674, ion_pos = 423427,
                    coord_system=local.LocalSlabCoords(4., 5., 1., 2., 1.),
                    bins = op_bins)
            op_dist /= float(np.sum(op_dist))
            op_dist_t.append(op_dist)

        legend = []
        for op_dist in op_dist_t:
            fig, = plt.plot((op_bins[0:-1] + op_bins[1:]) / 2, op_dist)
            legend.append(fig)
        plt.legend(legend, op_legend)
        plt.savefig("/home/jhaberstroh/Dropbox/Physics/subgroup/2014-10-28/TMV/cool_file.png")
