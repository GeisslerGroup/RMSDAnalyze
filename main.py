import RMSDAnalyze.grid
import RMSDAnalyze.atomic
import RMSDAnalyze.convert
import argparse
import ConfigParser
import h5py
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-trr", nargs=2, type=str, help="Load program from .trr and .tpr file (respectively), overwriting contents of configuration HDF5 file")
    parser.add_argument("-gro", type=str, help="Load program from .gro file, overwriting contents of configuration HDF5 file")
    parser.add_argument("-cfg", type=str, default=".RMSDAnalyze", help="Configuration file to load options from, default is .RMSDAnalyze")
    parser.add_argument("-rmsd_grid", type=str, nargs='+', choices=['display','png'], default=None, help="Output type for rmsd_grid")
    parser.add_argument("-rmsd_type", type=str, choices=['rmsd','angle'], default='rmsd',help='Type of RMSD calculation to perform')
    parser.add_argument("-rmsd_grid_file", type=str, default=None, help="Filename base for RMSD grid. If none passed, default file name is \'default_xxps.png\'.")
    parser.add_argument("-rmsd_dt", type=float, default=3.0, help="dt between RMSD calculations, converted to frames by rounding")
    parser.add_argument("-rmsd_activityparm", nargs=2, type=float, help="Parameters to convert RMSD to activity, namely r0 and k (cutoff distance and logistic rate, respectively)")
    parser.add_argument('-rmsd_scaletime', action='store_true', help="Scale time by sqrt(rmsd_dt)")
    parser.add_argument('-rmsd_colorrange', nargs=2, type=float, help="Min and max value to map to a color")
    parser.add_argument("-rmsd_out", type=str, help="File output destination for RMSD calculation")

    parser.add_argument("-op_grid", type=str, nargs='+', choices=['display','png'], default=None, help="Output type for op_grid")
    parser.add_argument("-op_type", type=str, choices=['q6','density'], default='density',help='Type of RMSD calculation to perform')

    parser.add_argument("-act_out", type=str, help="File output destination for activity XYZ calculation")
    parser.add_argument("-act_dt", type=float, help="dt between activity calculations, converted to frames by rounding")
    parser.add_argument("-act_frames", type=int, default=1, help="Frames (from the beginning of the file) to use for activity calculation")
    parser.add_argument("-act_openflag", choices=['w','a'], default='w', help="Open flag for writing to file. Use 'a' to append for sweeping activity")
    parser.add_argument("-act_dist", type=float, default=.5, help="Characteristic distance to measure activity")
    parser.add_argument("-act_slope", type=float, default=3., help="Rate for the exponential in the activity order parameter")
    parser.add_argument("-gro_copy", nargs=2, type=str, help="Output location for gro file")
    parser.add_argument("-gro_frames", nargs="+", type=int, help="array of frames to copy. Note: optimized for writing early frames")
    args = parser.parse_args()

    config = ConfigParser.RawConfigParser()
    config.read(args.cfg)
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    colormap= config.get('plotting','colormap')
    colormap= plt.cm.get_cmap(colormap)

    if args.gro and args.trr:
        raise ValueError("Cannot load from both .gro and .trr; select only one!")
    if args.trr != None:
        RMSDAnalyze.convert.trr2hdf5(args.trr, hdffile, hdfatom)
    if args.gro != None:
        RMSDAnalyze.convert.gro2hdf5(args.gro, hdffile, hdfatom)



    if args.gro_copy:
        print "Creating a copy of frames {} of {} with grocopy".format(args.gro_frames, args.gro_copy[0])
        RMSDAnalyze.convert.copygro(args.gro_copy[0], args.gro_copy[1], args.gro_frames)

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        print ds.shape
        print ds.attrs["dt"]
        if args.rmsd_out:
            rmsd_dT = int(round(args.rmsd_dt / ds.attrs["dt"]))
            RMSDAnalyze.atom.AtomRMSD(ds,args.rmsd_out, rmsd_dT)
        if args.act_out:
            print "activity parameters: r0 = {}, k = {}".format(args.act_dist, args.act_slope)
            activity_dT = int(round(args.act_dt / ds.attrs["dt"]))
            RMSDAnalyze.atom.AtomActivity(ds,args.act_out, args.act_openflag, args.act_frames, activity_dT, args.act_dist, args.act_slope)

        if args.rmsd_grid:
            rmsd_dT = int(round(args.rmsd_dt / ds.attrs["dt"]))

            rmsd_lambda = None
            if args.rmsd_activityparm:
                print "Creating RMSDLambda to compute activity"
                rmsd_lambda = RMSDAnalyze.grid.RMSDLambda(
                                        b_activity = True,                 \
                                        b_scaletime = args.rmsd_scaletime, \
                                        rmsd_delay = args.rmsd_dt,         \
                                        cutoff=args.rmsd_activityparm[0],  \
                                        sharpness=args.rmsd_activityparm[1])
                rmsd_lambda.SetTitle(args.rmsd_type)
            elif args.rmsd_scaletime:
                print "Creating RMSDLambda to scale time"
                rmsd_lambda = RMSDAnalyze.grid.RMSDLambda(
                                        b_activity = False,                \
                                        b_scaletime = args.rmsd_scaletime, \
                                        rmsd_delay = args.rmsd_dt)
                rmsd_lambda.SetTitle(args.rmsd_type)
            if 'png' in args.rmsd_grid and args.rmsd_grid_file == None:
                args.rmsd_grid_file = "default_{}ps".format(args.rmsd_dt)
            RMSDAnalyze.grid.GridRMSDRadial(ds,rmsd_dT, rmsd_type=args.rmsd_type, colorrange=args.rmsd_colorrange, \
                    display_type = args.rmsd_grid, file_name=args.rmsd_grid_file, rmsd_lambda = rmsd_lambda, \
                    colormap=colormap)

if __name__ == "__main__":
    main()
