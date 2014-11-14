import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import RMSDAnalyze.local as local
import RMSDAnalyze.grid as grid
import RMSDAnalyze.op as op
import RMSDAnalyze.coords as coords
import ConfigParser
import h5py
import numpy as np
import logging
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('op_type', choices=['rmsd','angle','diffusion'])
    parser.add_argument('local_type', choices=['mean','dist'])
    args = parser.parse_args()

    op_type = args.op_type
    out = args.local_type
    diffusion = False
    if args.op_type == 'diffusion':
        op_type = 'rmsd'
        diffusion = True

    #logging.basicConfig(level=logging.DEBUG)
    config = ConfigParser.RawConfigParser()
    config.read('RMSDAnalyze.cfg')
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    outdir = config.get('plotting','output_dir')
    plt.hold(True)
    logging.debug('Local type: {}'.format(out))
    assert( out == 'mean' or out == 'dist' )

    hex = op.HexPBC([ 22.34405, 18.91217,   9.91809 ,
                         0.00000,  0.00000, -10.91895 ,
                         0.00000, -0.00000,  -0.00000 ])
    center = hex.GetPBCCenter()

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        ctr_1 = -1.5
        ctr_2 = -1.0
        my_coords = [coords.LocalSlabCoords(3., 4., 
                            ctr_1 - .5, ctr_1 + .5, 1.0),
                     coords.LocalSlabCoords(2., 3., 
                            ctr_2 - .5, ctr_2 + .5, 1.0),
                     coords.LocalSlabCoords(9., 10., 3., 4., 1.),
                     coords.LocalSlabCoords(-.5, .5, -.5, .5, 1.)]
        plot_color = ['b', 'r', '#984ea3', 'k']
        cs_legend = ['chromophore_inner', 
                'chromophore_outer',
                'external', 
                'central pore']
        plots = []
        for i,coords in enumerate(my_coords):
            logging.debug("Plotting for {}".format(cs_legend[i]))
            if op_type == 'rmsd':
                op_bins = np.linspace(0,1,1000)
            elif op_type == 'angle':
                op_bins = np.linspace(0,3.4,1000)
            else:
                raise NotImplementedError(
                        "Plotting range unkown for op_type: {}".
                        format(op_type))
            op_dist_t = []
            op_legend = []
            if out == 'dist':
                times = [1, 5, 20]
            if out == 'mean':
                times = [1, 2, 3, 4] + range(5, 41, 5)
            for delay in times:
                dt = "{}ps".format(delay*.1)
                print("Using dt = {}".format(dt))
                #rmsd_lambda = op.RMSDLambda( b_activity = False,                 
                #                               b_scaletime = False, 
                #                               rmsd_delay = delay)
                #                        cutoff=args.rmsd_activityparm[0],  
                #                        sharpness=args.rmsd_activityparm[1])
                rmsd_lambda = None
                op_legend.append(dt)

                op_dist = local.LocalOP(ds, dynamic_step = delay, op_type=op_type,
                        rmsd_lambda = rmsd_lambda, water_pos=83674, ion_pos = 423427,
                        coord_system=coords, bins = op_bins, pbc = hex, center = center,
                        scaletime = diffusion)
                op_dist /= float(np.sum(op_dist))
                op_dist_t.append(op_dist)

            if out == 'dist':
                plots = []
                leg_label = ["{}ps".format(delay * .1) for delay in times]
                for op_dist in op_dist_t:
                    fig, = plt.plot((op_bins[0:-1] + op_bins[1:]) / 2, op_dist, '-o', linewidth=3)
                    plots.append(fig)
                plt.legend(plots, leg_label)
                plt.savefig("{}/{}_dist_{}.png".format(outdir, args.op_type, cs_legend[i]))
                plt.clf()
    
            if out == 'mean':
                op_bins_ctr = (op_bins[0:-1] + op_bins[1:]) / 2
                op_t = [np.sum(op_bins_ctr * op_dist) for op_dist in op_dist_t]
                if not scaletime:
                    op_t  = [0] + [np.sum(op_bins_ctr * op_dist) for op_dist in op_dist_t]
                    times = [0] + times
                print op_t
                plot, = plt.plot([time * .1 for time in times], op_t, 
                        color=plot_color[i])
                plots.append(plot)

        if out == 'mean':
            plt.legend(plots, cs_legend, loc=4)
            plt.savefig("{}/TMV-mean{}.png".format(outdir, args.op_type))

