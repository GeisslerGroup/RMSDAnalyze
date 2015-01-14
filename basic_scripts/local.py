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

DT=.1
#logging.basicConfig(level=logging.DEBUG)

if __name__=="__main__":
    op_types = ['diffusion','rmsd','angle1','angle2']
    parser = argparse.ArgumentParser()
    parser.add_argument('op_type', choices=op_types)
    parser.add_argument('local_type', choices=['mean','dist'])
    parser.add_argument("-cfg", default='.RMSDAnalyze', type=str, help="Configuration file location. (Default = .RMSDAnalyze)")
    parser.add_argument("-test", action='store_true', help="Run in testing mode?")
    
    args = parser.parse_args()

    op_type = args.op_type
    out = args.local_type
    diffusion = False
    legendre = 1

    if args.op_type == 'diffusion':
        op_type = 'rmsd'
        diffusion = True
    if args.op_type == 'angle1':
        op_type = 'angle'
        legendre = 1
    if args.op_type == 'angle2':
        op_type = 'angle'
        legendre = 2

    #logging.basicConfig(level=logging.DEBUG)
    config = ConfigParser.RawConfigParser()
    config.read(args.cfg)
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    system_type = config.get('local','system_type')
    outdir = config.get('plotting','output_dir')
    txtdir = config.get('plotting','text_dir')
    plt.hold(True)
    logging.debug('Local type: {}'.format(out))
    assert( out == 'mean' or out == 'dist' )


    system_types = ['tmv','water']
    if not system_type in system_types:
        raise ValueError("system_type ({}) not found. Use one of {}".
                format(system_type, system_types))
    
    if system_type == 'tmv':
        water_pos=83674
        ion_pos = 423427
        hex = op.HexPBC([ 22.34405, 18.91217,   9.91809 ,
                            0.00000,  0.00000, -10.91895 ,
                            0.00000, -0.00000,  -0.00000 ])
        center = hex.GetPBCCenter()
        stars = [[0, 0], 
                [-2.5, -1.],
#               [-2.5, -0.],
#               [-3.5, -0.], 
#               [-4.5, -0.],
                [-5.5, -0.],
                [-9.5, 3.5]]
        my_coords = []
        for star_pos in stars:
            xm, xp = star_pos[0] - .5, star_pos[0] + .5
            ym, yp = star_pos[1] - .5, star_pos[1] + .5
            new_coord = coords.LocalSlabCoords(xm, xp, ym, yp, 1.0)
            my_coords.append(new_coord)
        plot_color = ['k', #plot_color = ['k',     
                '#e41a1c', #        '#d7191c',
#               '#377eb8', #        '#fdae61',
#               '#4daf4a', #        '#ffffbf',
#               '#984ea3', #        '#abdda4',
                '#ff7f00', #        '#2b83ba',
                '#a65628'] #        '#7570b3']
        cs_legend = ['central-pore',
                'pore-low', 
#                'pore-mid',
#                'disc-shallow',
#                'disc-mid',
                'disc-deep',
                'external' ]
        sys_prefix = 'TMV'
        if args.test:
            print "Running in Testing Mode"
            my_coords  = [my_coords[4],  my_coords[5]]
            plot_color = [plot_color[0], plot_color[1]]
            cs_legend  = [cs_legend[0], cs_legend[1]]
            times = [times[0], 23]
            print plot_color

    if system_type=='water':
        water_pos = 0
        ion_pos = 12126
        hex = op.HexPBC([ 5.05657, 5.05657, 5.05657,
                          0.00000, 0.00000, 0.00000 ,
                          0.00000, 0.00000, 0.00000 ])
        center = hex.GetPBCCenter()
        stars = [[0, 0]], 
        my_coords = [coords.LocalSlabCoords(-1.5, 1.5, -1.5, 1.5, 3.)]
        plot_color = ['k'] 
        cs_legend = ['control']
        sys_prefix = 'water'

    if out == 'dist':
        times = [1, 5, 20]
    if out == 'mean' and op_type == 'rmsd':
        times = [1, 2, 3, 4] + range(5, 111, 10)
        if not diffusion:
            times = [0] + times
    if out == 'mean' and op_type == 'angle':
        times = range(0, 5) + range(5, 51, 5) + range(55, 116, 10)

    try:
        hex
        center
        plot_color
        cs_legend
        stars
        times
        water_pos
        ion_pos
        sys_prefix
    except NameError:
        raise ValueError("system_type = water is not fully specified."+
                "Please check source and run again.")

    params = 3
    out_data = np.zeros((len(times), params * len(my_coords) + 1))
    out_data[:,0] = np.array(times) * DT
    plots = []

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        for i, coords in enumerate(my_coords):
            logging.debug("Plotting for {}".format(cs_legend[i]))
            if op_type == 'rmsd':
                op_bins = np.linspace(0,1,200)
            elif op_type == 'angle':
                op_bins = np.linspace(-1,1,200)
            else:
                raise NotImplementedError(
                        "Plotting range unkown for op_type: {}".
                        format(op_type))

            # Compute the distributions for every time slide for one 
            # coordinate system
            op_dist_t = []
            N_t = []
            op_legend = []
            for delay in times:
                dt = "{}ps".format(delay * DT)
                print("Using dt = {}".format(dt))
                op_legend.append(dt)

                if delay == 0:
                    op_dist = np.zeros(len(op_bins)-1)
                    if op_type == 'rmsd':
                        op_dist[0] = 1.
                    elif op_type == 'angle':
                        op_dist[-1] = 1.
                    op_dist_t.append(op_dist)
                    op_legend.append(dt)
                    N_t.append(0)

                else:
                    op_dist = local.LocalOP(ds, dynamic_step = delay, op_type=op_type,
                            rmsd_lambda = None, water_pos=water_pos, ion_pos = ion_pos,
                            coord_system=coords, bins = op_bins, pbc = hex, center = center,
                            scaletime = diffusion, legendre=legendre)
                    N_t.append(np.sum(op_dist))
                    op_dist /= float(np.sum(op_dist))
                    op_dist_t.append(op_dist)

            # Output the analysis for all time slices of this coordinate
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
                op_dist = np.array(op_dist)
                op_bins_ctr = np.array((op_bins[0:-1] + op_bins[1:]) / 2)
                op_t = np.array([np.sum(op_bins_ctr * op_dist) for op_dist in op_dist_t])
                err_t = np.array([np.sum(np.square(op_bins_ctr) * op_dist) for op_dist in op_dist_t])
                N_t = np.array(N_t)
                err_t -= np.square(op_t)
                for err_ind in xrange(len(err_t)):
                    print("Value {0}: {1:.3f} +/- {2:.3f}"
                            .format(err_ind, op_t[err_ind], err_t[err_ind]))
                err_t /= np.sqrt(N_t / 100)
                err_t[N_t == 0] = 0
                if diffusion: 
                    err_t *= (2 * op_t)
                    op_t  = np.square(op_t)
                for err_ind in xrange(len(err_t)):
                    print("Value {}: {:.3f} +/- {:.3f} \t N={}"
                            .format(err_ind, op_t[err_ind], err_t[err_ind], N_t[err_ind]))
                print('op value: {}'.format(op_t))
                print('err value: {}'.format(err_t))
                #plot, = plt.plot([time * DT for time in times], op_t, 
                #        color=plot_color[i])
                plot = plt.errorbar(x=[time * DT for time in times], y=op_t, 
                        yerr = err_t, elinewidth=1, color=plot_color[i])
                out_data[:,params*i+1] = op_t
                out_data[:,params*i+2] = err_t
                out_data[:,params*i+3] = N_t
                plots.append(plot)

        if out == 'mean':
            plt.legend(plots, cs_legend, loc=3)
            plt.savefig("{}/{}-mean{}.png".format(outdir, sys_prefix, args.op_type))
            np.savetxt("{}/{}-mean{}.csv".format(txtdir, sys_prefix, args.op_type),
                    out_data, header=" ".join(['time'] + cs_legend) )

