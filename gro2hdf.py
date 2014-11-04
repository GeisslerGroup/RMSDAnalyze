import RMSDAnalyze.convert as conv
import argparse
import ConfigParser

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gro", type=str, help="File to convert.  e.g. -gro_copy traj.gro")
    parser.add_argument("-cfg", default='.RMSDAnalyze', type=str, help="Configuration file location. (Default = .RMSDAnalyze)")
    args = parser.parse_args()


    config = ConfigParser.RawConfigParser()
    config.read(args.cfg)
    hdffile = config.get('HDF','file')
    hdfatom = config.get('HDF','atom_dset')
    
    conv.gro2hdf5(args.gro, hdffile, hdfatom)
