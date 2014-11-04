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

#    with h5py.File(hdffile,'r') as h5:
#        ds = h5[hdfatom]
#
#        my_coords = [local.LocalSlabCoords(3., 4., -2.0, -1.0, 1.0),
#                local.LocalSlabCoords(9., 10., 3., 4., 1.),
#                local.LocalSlabCoords(-.5, .5, -.5, .5, 1.)]
#        cs_legend = ['chromophore', 
#                'external', 'central pore']
#        plots = []
#        for i,coords in enumerate(my_coords):
#
#            op_bins = np.linspace(0,1,1000)
#            op_dist_t = []
#            op_legend = []
#            times = [1,2,3,4] + range(5,51,5)
#            for i in times:
#                dt = "{}ps".format(i*.1)
#                print("Using dt = {}".format(dt))
#                rmsd_lambda = grid.RMSDLambda( b_activity = False,                 
#                                               b_scaletime = False, 
#                                               rmsd_delay = i)
#                #                        cutoff=args.rmsd_activityparm[0],  
#                #                        sharpness=args.rmsd_activityparm[1])
#                op_legend.append(dt)
#                op_dist = local.LocalOP(ds, dynamic_step = i, op_type='rmsd',
#                        rmsd_lambda = rmsd_lambda, water_pos=83674, ion_pos = 423427,
#                        coord_system=coords,
#                        bins = op_bins)
#                op_dist /= float(np.sum(op_dist))
#                op_dist_t.append(op_dist)
#    
#            op_bins_ctr = (op_bins[0:-1] + op_bins[1:]) / 2
#            op_t  = [0] + [np.sum(op_bins_ctr * op_dist) for op_dist in op_dist_t]
#            times = [0] + times
#            
#            print op_t
#            plot, = plt.plot([time * .1 for time in times], op_t)
#            plots.append(plot)
#
#        plt.legend(plots, cs_legend, loc=4)
#        plt.savefig("/home/jhaberstroh/Dropbox/Physics/subgroup/2014-10-28/TMV/meanRMSD.png")

# PLOT THE DISTRIBUTIONS
    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]

        my_coords = [local.LocalSlabCoords(3., 4., -2.0, -1.0, 1.0)]
        cs_legend = ['chromophore']
        plots = []
        for i,coords in enumerate(my_coords):

            op_bins = np.linspace(0,.5,50)
            op_dist_t = []
            op_legend = []
            times = range(5, 51, 10)
            for i in times:
                dt = "{}ps".format(i*.1)
                print("Using dt = {}".format(dt))
                rmsd_lambda = grid.RMSDLambda( b_activity = False,                 
                                               b_scaletime = True, 
                                               rmsd_delay = i * .1)
                #                        cutoff=args.rmsd_activityparm[0],  
                #                        sharpness=args.rmsd_activityparm[1])
                op_legend.append(dt)
                op_dist = local.LocalOP(ds, dynamic_step = i, op_type='rmsd',
                        rmsd_lambda = rmsd_lambda, water_pos=83674, ion_pos = 423427,
                        coord_system=coords,
                        bins = op_bins)
                op_dist /= float(np.sum(op_dist))
                op_dist_t.append(op_dist)

            plots = []
            leg_label = ["{}ps".format(i * .1) for i in times]
            for op_dist in op_dist_t:
                fig, = plt.plot((op_bins[0:-1] + op_bins[1:]) / 2, op_dist, '-o', linewidth=3)
                plots.append(fig)
            plt.legend(plots, leg_label)
            plt.savefig("/home/jhaberstroh/Dropbox/Physics/subgroup/2014-10-28/TMV/diffusion_dist_chromo.png")
