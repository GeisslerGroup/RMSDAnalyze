import RMSDAnalyze.local as local
import RMSDAnalyze.grid as grid
import ConfigParser
import h5py

if __name__=="__main__":
    config = ConfigParser.RawConfigParser()
    config.read('.RMSDAnalyze')
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]

        #rmsd_lambda = grid.RMSDLambda( b_activity = True,                 
        #                        b_scaletime = args.rmsd_scaletime, 
        #                        rmsd_delay = args.rmsd_dt,         
        #                        cutoff=args.rmsd_activityparm[0],  
        #                        sharpness=args.rmsd_activityparm[1])

        op_dist, op_bins = local.LocalOP(ds, dynamic_step = 5, op_type='rmsd',
                rmsd_lambda = None, water_pos=83674, ion_pos = 423427,
                coord_system=[local.LocalSlabCoords(4., 5., 1., 2., 1.)],
                n_bins = 100)

        print op_dist, op_bins
