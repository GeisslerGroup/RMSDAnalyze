import numpy as np
import argparse
import ConfigParser
import os
import h5py
import MDAnalysis
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import scipy.sparse as mtx
import numpy.linalg as LA


def tail(f, n, offset=0):
  command = "tail -n {} {}".format(n+offset, f)
  stdin,stdout = os.popen2(command)
  stdin.close()
  lines = stdout.readlines(); stdout.close()
  return lines[0:-offset]

def copygro(gro_in, gro_out, frames):
    with open(gro_in) as f:
        i = 0
        for l in f:
            l_arr = l.split()
            if "trjconv" in l_arr:
                t0 = float(l_arr[-1])
                print "t0 = {}".format(t0)
                atoms = int(f.next())
                print "atoms = {}".format(atoms)
                break
        for atom in xrange(atoms):
            f.next()
        box = f.next().strip()
        l1 = f.next()
        t1 = float(l1.split()[-1])
        print "box = {}".format(box)
        print "t1 = {}".format(t1)
        print "dt = {}".format(t1-t0)

    lf = tail(gro_in, 1, atoms+2)[0]
    tf = float(lf.split()[-1])
    print "tf = {}".format(tf)
    N  = int(round((tf / (t1-t0)))) + 1

    if np.any(np.array(frames) >= N):
        raise ValueError("Frame value too high in copygro, frames must be less than {}".format(N))
    
    f_out = open(gro_out,'w')
    for frame in frames:
        f_out = open(gro_out,'a')
        with open(gro_in) as f:
            n = 0
            while n <= frame:
                print 'Frame ', n, "/", N
                l0 = f.next()
                l1 = f.next()
                if n == frame:
                    f_out.write(l0)
                    f_out.write(l1)
                for atom in xrange(atoms):
                    l = f.next()
                    if n == frame:
                        f_out.write(l)
                lf = f.next()
                if n == frame:
                    f_out.write(lf)
                n += 1
        f_out.close()


def trr2hdf5(trr_in, hdf_out, hdf_atom):
    print "Loading universe..."
    universe=MDAnalysis.Universe(trr_in[1], trr_in[0])
    print "Getting the dt step..."
    print universe.trajectory.dt

def gro2hdf5(gro_in, hdf_out, hdf_atom):
    with open(gro_in) as f:
        i = 0
        for l in f:
            l_arr = l.split()
            if "trjconv" in l_arr:
                t0 = float(l_arr[-1])
                print "t0 = {}".format(t0)
                atoms = int(f.next())
                print "atoms = {}".format(atoms)
                break
        for atom in xrange(atoms):
            f.next()
        box = f.next().strip()
        l1 = f.next()
        t1 = float(l1.split()[-1])
        print "box = {}".format(box)
        print "t1 = {}".format(t1)
        print "dt = {}".format(t1-t0)

    lf = tail(gro_in, 1, atoms+2)[0]
    tf = float(lf.split()[-1])
    print "tf = {}".format(tf)
    N  = int(round((tf / (t1-t0)))) + 1

    h5 = h5py.File(hdf_out,'w')
    ds = h5.create_dataset(hdf_atom, (N,atoms,3), dtype='f32')
    with open(gro_in) as f:
        for n in xrange(N):
            print 'Frame ', n, "/", N
            l0 = f.next()
            l1 = f.next()
            for atom in xrange(atoms):
                if atom % 100000 == 0:
                    print "\tatom",atom
                l = f.next()
                pos = l[20:].strip()
                pos = np.array([float(xyz) for xyz in pos.split()])
                ds[n,atom,:] = pos[:]
            lf = f.next()
    h5.close()


def AtomRMSD(data_tik, output_file, rmsd_step, time_init=0.0, time_step=0.0):
    f = open(output_file,'w') 
    for t0 in xrange(data_tik.shape[0]-rmsd_step):
        f.write("{}\n".format(data_tik.shape[1]))
        f.write("RMSD generated in XYZ by RMSDAnalyze.py @jhaberstroh: t = {}\n".format(time_init + time_step * t0))
        print "outputting time {} of {}".format(t0, data_tik.shape[0]-rmsd_step-1)
        r_ik = data_tik[t0+rmsd_step,:,:] - data_tik[t0,:,:]
        rmsd_i = np.sqrt(np.sum(np.square(r_ik), axis=1))
        rmsd_i[rmsd_i > 5.0] = 0
        #plt.hist(rmsd_i, bins=100)
        #plt.show()
        for rmsd in rmsd_i:
            f.write("C\t{}\t0.0\t0.0\n".format(rmsd))
    f.close()


def AtomActivity(data_tik, output_file, openflag, frames, activity_step, r0, k):
    f = open(output_file, openflag)
    if frames+activity_step > data_tik.shape[0]:
        raise AttributeError("Frames need for activity calculation ({}) exceeds available frames ({})".format(frames+activity_step, data_tik.shape[0]))
    print "Frames: {}".format(frames)
    for t0 in xrange(frames):
        f.write("{}\n".format(data_tik.shape[1]))
        f.write("Activity generated in XYZ by RMSDAnalyze.py @jhaberstroh: t = {}\n".format(frames))
        r_ik = data_tik[t0+activity_step,:,:] - data_tik[t0,:,:]
        rmsd_i = np.sqrt(np.sum(np.square(r_ik), axis=1))

        activity_i = np.reciprocal(1 + np.exp( -k * (rmsd_i - r0) ))
        for activity in activity_i: 
            f.write("C\t{}\t0.0\t0.0\n".format(activity))
    f.close()

class RMSDLambda:
    def __init__(self, b_scaletime, b_activity, rmsd_delay=0, cutoff=0, sharpness=0,title=None):
        self.b_scaletime = b_scaletime
        self.b_activity = b_activity
        self.rmsd_delay = rmsd_delay
        self.sharpness = sharpness
        self.cutoff = cutoff
        self.title=title
    def __call__ (self, rmsd):
        if self.b_scaletime:
            rmsd /= np.sqrt(self.rmsd_delay)
        if self.b_activity:
            return np.reciprocal(1 + np.exp( -self.sharpness * (rmsd - self.cutoff)))
        else:
            return rmsd
    def SetTitle(self, rmsd_type='rmsd'):
        print("rmsd_type={}, activity={}, scaletime={}".format(rmsd_type, self.b_activity, self.b_scaletime))
        if rmsd_type == 'angle':
            if self.b_activity:
                print "SETTING TITLE GOOD"
                self.title='Mean angular activity, cutoff={}rad, k={}rad^-1, t_delay={}ps'.format(self.cutoff, self.sharpness, self.rmsd_delay/10.0)
            else:
                print "SETTING TITLE BAD"
                self.title='Mean angular displacement, t_delay={}ps'.format(self.rmsd_delay/10.0)
        elif rmsd_type == 'rmsd':
            if self.b_activity and self.b_scaletime:
                self.title='Mean activity (scaled), cutoff={}nm/ps^(.5), k={}ps^(.5)/nm, t_delay={}ps'.format(self.cutoff, self.sharpness, self.rmsd_delay/10.0)
            elif self.b_activity:
                self.title='Mean activity, cutoff={}nm, k={}nm^-1, t_delay={}ps'.format(self.cutoff, self.sharpness, self.rmsd_delay/10.0)
            elif self.b_scaletime:
                self.title='Scaled RMSD (root diffusion), t_delay={}ps'.format(self.rmsd_delay/10.0)
            else:
                self.title='Mean RMSD, t_delay={}ps'.format(self.rmsd_delay/10.0)



def ComputeRMSD(r0_ik, r1_ik, pre_cutoff=None, post_cutoff=None, rmsd_lambda = None):
    dr_ik = r0_ik- r1_ik 
    rmsd_i = np.sqrt(np.sum(np.square(dr_ik), axis=1))
    if pre_cutoff:
        rmsd_i[rmsd_i > pre_cutoff] = None
    if rmsd_lambda:
        rmsd_i = rmsd_lambda(rmsd_i)
    if post_cutoff:
        rmsd_i[rmsd_i > post_cutoff] = None
    return rmsd_i

def ComputeRMSAngle(r0_ik, r1_ik, rmsd_lambda = None):
    r0_v = (r0_ik[0::3] - r0_ik[1::3]) + \
           (r0_ik[0::3] - r0_ik[2::3])
    r1_v = (r1_ik[0::3] - r1_ik[1::3]) + \
           (r1_ik[0::3] - r1_ik[2::3])
    r0_norm = LA.norm(r0_v, axis=1)
    r0_norm.shape = (r0_norm.shape[0],1)
    r0_v /= r0_norm
    r1_norm = LA.norm(r1_v, axis=1)
    r1_norm.shape = (r1_norm.shape[0],1)
    r1_v /= r1_norm
    theta = np.arccos(np.sum(r0_v*r1_v,axis=1))
    return rmsd_lambda(theta)

def ComputeQ6(atoms_i, cutoff_A):
    return np.ones(atoms_i.shape)


def GridOPRadial(data_tik, display_type=[], colorrange=None, rmsd_type='rmsd', \
                   file_name=None, rmsd_lambda=None, time_init=0.0, time_step=0.0,\
                   colormap=plt.cm.Spectral_r):
    if colorrange==None:
        colorrange=[None,None]
    atom_type = 'water'
    running_mean_mtx = None
    running_weight_mtx = None

    gridsize = [40,30]
    r_extent = 10.0
    z_extent = 4.0
    water_pos = 83674
    ion_pos  = 423427

    mtx_scale = max(gridsize[0] / r_extent, gridsize[1] / (2 * z_extent)) * 100

    for t0 in xrange(data_tik.shape[0]):
        print "outputting time {} of {}".format(t0, data_tik.shape[0])
        display_histogram=False
        
        # Select the atom subset to run computation on
        if atom_type == 'water':
            atoms_i = data_tik[t0,water_pos:ion_pos,:]
        else:
            raise ValueError("GridRMSDRadial passed atom_type that is not known: {}".format(atom_type))

        # Run the appropriate computation on those atoms
        if op_type == 'q6':
            op_i = ComputeQ6(atoms_i, cutoff_A=4.0)
        elif op_type == 'density':
            op_i = np.ones(atoms_i.shape)
        else:
            raise ValueError("GridOPRadial passed op_type that is not known: {}".format(op_type))

        if display_histogram:
            plt.hist(op_i)
            plt.title("RMSD Post-histogram")
            plt.show()

        # Center the initial frame to compute the histogram
        center_k = np.mean(data_tik[t0,:,:], axis=0)
        r_ik = atoms_i - center_k
        r_cyl_i = np.sqrt( np.square(r_ik[:,0]) + np.square(r_ik[:,1]))
        z_cyl_i = r_ik[:,2]
        
        # Select only oxygens from water, and only within a spatial domain specified by sub
        if atom_type == 'water':
            r = r_cyl_i[::3]
            z = z_cyl_i[::3]
            if rmsd_type=='rmsd':
                C = op_i[::3]
            elif rmsd_type=='angle':
                C = op_i
            else:
                raise ValueError("GridRMSDRadial passed rmsd_type that is not known: {}".format(rmsd_type))

        # Compute the extent of the box
        extent = [0, r_extent, -z_extent, z_extent]
        # Create a boolean vector to select atoms only within the bounds of the extents
        sub = (r < r_extent) * (np.abs(z) < z_extent)


        # Compute the mean RMSD per box and the weight of the box
        hexmeanplt   = plt.hexbin(r[sub], z[sub], C=C[sub], \
                cmap=colormap, mincnt=0, gridsize = gridsize, extent=extent) 
        cb = plt.colorbar(hexmeanplt, spacing='uniform',extend='max')
        plt.clf()
        hexweightplt = plt.hexbin(r[sub], z[sub], \
                cmap=colormap, mincnt=0, gridsize = gridsize, extent=extent) 
        cb = plt.colorbar(hexweightplt, spacing='uniform',extend='max')
        plt.clf()


        #plt.show()
        hex_mean   = hexmeanplt.get_array()
        hex_weight = hexweightplt.get_array()
        #print "mean shape, weight shape: {}, {}".format(hex_mean.shape, hex_weight.shape)
        mean_pos   =  mtx_scale * (hexmeanplt.get_offsets()   + np.array([0,z_extent]))
        weight_pos =  mtx_scale * (hexweightplt.get_offsets() + np.array([0,z_extent]))
        mean_pos.astype(int)
        weight_pos.astype(int)


        if running_mean_mtx == None and running_weight_mtx == None:
            mean_mtx   = mtx.csr_matrix( (hex_mean  , (mean_pos[:,0],   mean_pos[:,1]  )), dtype='f')
            weight_mtx = mtx.csr_matrix( (hex_weight, (weight_pos[:,0], weight_pos[:,1])) )
            running_mean_mtx   = mean_mtx.multiply(weight_mtx)
            running_weight_mtx = weight_mtx
        else:
            mean_mtx   = mtx.csr_matrix( (hex_mean  , (mean_pos[:,0],   mean_pos[:,1]  )), dtype='f')
            weight_mtx = mtx.csr_matrix( (hex_weight, (weight_pos[:,0], weight_pos[:,1])) )
            running_mean_mtx   = running_mean_mtx + mean_mtx.multiply(weight_mtx)
            running_weight_mtx = running_weight_mtx + weight_mtx


    running_mean_mtx = running_mean_mtx / running_weight_mtx

    running_mean_mtx   = running_mean_mtx.tocoo()
    running_weight_mtx = running_weight_mtx.tocoo()
    #print "RUNNING MEAN TYPE: {}".format(type(running_mean_mtx))
    #print "RUNNING MEAN DATA: {}".format(running_mean_mtx)
    vertsRMSD = np.vstack((running_mean_mtx.row,running_mean_mtx.col)).astype(float).T / mtx_scale
    vertsRMSD -= np.array([0,z_extent])
    weightedRMSD = running_mean_mtx.data
    verts_density = np.vstack((running_weight_mtx.row,running_weight_mtx.col)).astype(float).T / mtx_scale
    verts_density -= np.array([0,z_extent])
    mean_density = running_weight_mtx.data / verts_density[:,0]
    

    # Plot the protein structure
    center_k = np.mean(data_tik[:,0:water_pos,:], axis=(0,1))
    protein_r = np.sqrt(np.square(data_tik[0,0:water_pos,0] - center_k[0]) + \
                        np.square(data_tik[0,0:water_pos,1] - center_k[1]) )
    protein_z = data_tik[0,0:water_pos,2] - center_k[2]
    sub = (protein_r < r_extent) * (np.abs(protein_z) < z_extent)



    print "PLOTTING THE FINAL RESULT!!"
    counts = hexmeanplt.get_array()
    ncnts = np.count_nonzero(np.power(10,counts))
    binx = vertsRMSD[:,0]
    biny = vertsRMSD[:,1]
    plt.subplot(2,1,1)
    RMSD_out = plt.hexbin(binx, biny, C=weightedRMSD, vmin=colorrange[0], vmax=colorrange[1], \
                          cmap=colormap, gridsize=gridsize, extent=extent)
    plt.plot(protein_r[sub], protein_z[sub], '-k', alpha=.3)
    cb = plt.colorbar(RMSD_out, spacing='uniform',extend='max')

    if rmsd_lambda.title:
        title=rmsd_lambda.title
    plt.title(title)
    plt.xlabel("R, cylindrical radius from center of disc (nm)")
    plt.ylabel("Z, vertical height (nm)")

    plt.subplot(2,1,2)
    binx = verts_density[:,0]
    biny = verts_density[:,1]
    density_out = plt.hexbin(binx,biny,C=mean_density/(data_tik.shape[0]-rmsd_step), \
                             cmap=colormap, gridsize=gridsize, extent=extent)
    cb = plt.colorbar(density_out, spacing='uniform',extend='max')
    plt.plot(protein_r[sub], protein_z[sub], '-k', alpha=.3)
    plt.title("Density (atom per frame)")
    plt.xlabel("R, cylindrical radius from center of disc (nm)")
    plt.ylabel("Z, vertical height (nm)")
    plt.tight_layout()

    if 'png' in display_type:
        if file_name:
            plt.savefig(file_name + '.png')
        else:
            plt.savefig('default.png')

    if 'display' in display_type:
        plt.show()


def GridRMSDRadial(data_tik, rmsd_step, display_type=[], colorrange=None, rmsd_type='rmsd', \
                   file_name=None, rmsd_lambda=None, time_init=0.0, time_step=0.0,\
                   colormap=plt.cm.Spectral_r):
    if colorrange==None:
        colorrange=[None,None]
    atom_type = 'water'
    running_mean_mtx = None
    running_weight_mtx = None

    gridsize = [40,30]
    r_extent = 10.0
    z_extent = 4.0
    water_pos = 83674
    ion_pos  = 423427

    mtx_scale = max(gridsize[0] / r_extent, gridsize[1] / (2 * z_extent)) * 100

    for t0 in xrange(data_tik.shape[0]-rmsd_step):
        print "outputting time {} of {}".format(t0, data_tik.shape[0]-rmsd_step-1)
        display_histogram=False
        
        # Select the atom subset to run computation on
        if atom_type == 'water':
            atoms_i = data_tik[t0,water_pos:ion_pos,:]
            atoms_f = data_tik[t0+rmsd_step,water_pos:ion_pos,:]
        else:
            raise ValueError("GridRMSDRadial passed atom_type that is not known: {}".format(atom_type))

        # Runn the appropriate computation on those atoms
        if rmsd_type == 'rmsd':
            rmsd_i = ComputeRMSD(atoms_i, atoms_f, pre_cutoff=5.0, rmsd_lambda=rmsd_lambda)
        elif rmsd_type == 'angle':
            if atom_type != 'water':
                raise ValueError("atom_type \'water\' must be used with rmsd_type = \'angle\'")
            rmsd_i = ComputeRMSAngle(atoms_i, atoms_f, rmsd_lambda=rmsd_lambda)
        else:
            raise ValueError("GridRMSDRadial passed rmsd_type that is not known: {}".format(rmsd_type))

        if display_histogram:
            plt.hist(rmsd_i)
            plt.title("RMSD Post-histogram")
            plt.show()

        # Center the initial frame to compute the histogram
        center_k = np.mean(data_tik[t0,:,:], axis=0)
        r_ik = atoms_i - center_k
        r_cyl_i = np.sqrt( np.square(r_ik[:,0]) + np.square(r_ik[:,1]))
        z_cyl_i = r_ik[:,2]
        
        # Select only oxygens from water, and only within a spatial domain specified by sub
        if atom_type == 'water':
            r = r_cyl_i[::3]
            z = z_cyl_i[::3]
            if rmsd_type=='rmsd':
                C = rmsd_i[::3]
            elif rmsd_type=='angle':
                C = rmsd_i
            else:
                raise ValueError("GridRMSDRadial passed rmsd_type that is not known: {}".format(rmsd_type))

        # Compute the extent of the box
        extent = [0, r_extent, -z_extent, z_extent]
        # Create a boolean vector to select atoms only within the bounds of the extents
        sub = (r < r_extent) * (np.abs(z) < z_extent)


        # Compute the mean RMSD per box and the weight of the box
        hexmeanplt   = plt.hexbin(r[sub], z[sub], C=C[sub], \
                cmap=colormap, mincnt=0, gridsize = gridsize, extent=extent) 
        cb = plt.colorbar(hexmeanplt, spacing='uniform',extend='max')
        plt.clf()
        hexweightplt = plt.hexbin(r[sub], z[sub], \
                cmap=colormap, mincnt=0, gridsize = gridsize, extent=extent) 
        cb = plt.colorbar(hexweightplt, spacing='uniform',extend='max')
        plt.clf()


        #plt.show()
        hex_mean   = hexmeanplt.get_array()
        hex_weight = hexweightplt.get_array()
        #print "mean shape, weight shape: {}, {}".format(hex_mean.shape, hex_weight.shape)
        mean_pos   =  mtx_scale * (hexmeanplt.get_offsets()   + np.array([0,z_extent]))
        weight_pos =  mtx_scale * (hexweightplt.get_offsets() + np.array([0,z_extent]))
        mean_pos.astype(int)
        weight_pos.astype(int)


        if running_mean_mtx == None and running_weight_mtx == None:
            mean_mtx   = mtx.csr_matrix( (hex_mean  , (mean_pos[:,0],   mean_pos[:,1]  )), dtype='f')
            weight_mtx = mtx.csr_matrix( (hex_weight, (weight_pos[:,0], weight_pos[:,1])) )
            running_mean_mtx   = mean_mtx.multiply(weight_mtx)
            running_weight_mtx = weight_mtx
        else:
            mean_mtx   = mtx.csr_matrix( (hex_mean  , (mean_pos[:,0],   mean_pos[:,1]  )), dtype='f')
            weight_mtx = mtx.csr_matrix( (hex_weight, (weight_pos[:,0], weight_pos[:,1])) )
            running_mean_mtx   = running_mean_mtx + mean_mtx.multiply(weight_mtx)
            running_weight_mtx = running_weight_mtx + weight_mtx


    running_mean_mtx = running_mean_mtx / running_weight_mtx

    running_mean_mtx   = running_mean_mtx.tocoo()
    running_weight_mtx = running_weight_mtx.tocoo()
    #print "RUNNING MEAN TYPE: {}".format(type(running_mean_mtx))
    #print "RUNNING MEAN DATA: {}".format(running_mean_mtx)
    vertsRMSD = np.vstack((running_mean_mtx.row,running_mean_mtx.col)).astype(float).T / mtx_scale
    vertsRMSD -= np.array([0,z_extent])
    weightedRMSD = running_mean_mtx.data
    verts_density = np.vstack((running_weight_mtx.row,running_weight_mtx.col)).astype(float).T / mtx_scale
    verts_density -= np.array([0,z_extent])
    mean_density = running_weight_mtx.data / verts_density[:,0]
    

    # Plot the protein structure
    center_k = np.mean(data_tik[:,0:water_pos,:], axis=(0,1))
    protein_r = np.sqrt(np.square(data_tik[0,0:water_pos,0] - center_k[0]) + \
                        np.square(data_tik[0,0:water_pos,1] - center_k[1]) )
    protein_z = data_tik[0,0:water_pos,2] - center_k[2]
    sub = (protein_r < r_extent) * (np.abs(protein_z) < z_extent)



    print "PLOTTING THE FINAL RESULT!!"
    counts = hexmeanplt.get_array()
    ncnts = np.count_nonzero(np.power(10,counts))
    binx = vertsRMSD[:,0]
    biny = vertsRMSD[:,1]
    plt.subplot(2,1,1)
    RMSD_out = plt.hexbin(binx, biny, C=weightedRMSD, vmin=colorrange[0], vmax=colorrange[1], \
                          cmap=colormap, gridsize=gridsize, extent=extent)
    plt.plot(protein_r[sub], protein_z[sub], '-k', alpha=.3)
    cb = plt.colorbar(RMSD_out, spacing='uniform',extend='max')

    if rmsd_lambda.title:
        title=rmsd_lambda.title
    plt.title(title)
    plt.xlabel("R, cylindrical radius from center of disc (nm)")
    plt.ylabel("Z, vertical height (nm)")

    plt.subplot(2,1,2)
    binx = verts_density[:,0]
    biny = verts_density[:,1]
    density_out = plt.hexbin(binx,biny,C=mean_density/(data_tik.shape[0]-rmsd_step), \
                             cmap=colormap, gridsize=gridsize, extent=extent)
    cb = plt.colorbar(density_out, spacing='uniform',extend='max')
    plt.plot(protein_r[sub], protein_z[sub], '-k', alpha=.3)
    plt.title("Density (atom per frame)")
    plt.xlabel("R, cylindrical radius from center of disc (nm)")
    plt.ylabel("Z, vertical height (nm)")
    plt.tight_layout()

    if 'png' in display_type:
        if file_name:
            plt.savefig(file_name + '.png')
        else:
            plt.savefig('default.png')

    if 'display' in display_type:
        plt.show()




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
        trr2hdf5(args.trr, hdffile, hdfatom)
    if args.gro != None:
        gro2hdf5(args.gro, hdffile, hdfatom)



    if args.gro_copy:
        print "Creating a copy of frames {} of {} with grocopy".format(args.gro_frames, args.gro_copy[0])
        copygro(args.gro_copy[0], args.gro_copy[1], args.gro_frames)

    with h5py.File(hdffile,'r') as h5:
        ds = h5[hdfatom]
        print ds.shape
        print ds.attrs["dt"]
        if args.rmsd_out:
            rmsd_dT = int(round(args.rmsd_dt / ds.attrs["dt"]))
            AtomRMSD(ds,args.rmsd_out, rmsd_dT)
        if args.act_out:
            print "activity parameters: r0 = {}, k = {}".format(args.act_dist, args.act_slope)
            activity_dT = int(round(args.act_dt / ds.attrs["dt"]))
            AtomActivity(ds,args.act_out, args.act_openflag, args.act_frames, activity_dT, args.act_dist, args.act_slope)

        if args.rmsd_grid:
            rmsd_dT = int(round(args.rmsd_dt / ds.attrs["dt"]))

            rmsd_lambda = None
            if args.rmsd_activityparm:
                print "Creating RMSDLambda to compute activity"
                rmsd_lambda = RMSDLambda(b_activity = True,                 \
                                         b_scaletime = args.rmsd_scaletime, \
                                         rmsd_delay = args.rmsd_dt,         \
                                         cutoff=args.rmsd_activityparm[0],  \
                                         sharpness=args.rmsd_activityparm[1])
                rmsd_lambda.SetTitle(args.rmsd_type)
            elif args.rmsd_scaletime:
                print "Creating RMSDLambda to scale time"
                rmsd_lambda = RMSDLambda(b_activity = False,                \
                                         b_scaletime = args.rmsd_scaletime, \
                                         rmsd_delay = args.rmsd_dt)
                rmsd_lambda.SetTitle(args.rmsd_type)
            if 'png' in args.rmsd_grid and args.rmsd_grid_file == None:
                args.rmsd_grid_file = "default_{}ps".format(args.rmsd_dt)
            GridRMSDRadial(ds,rmsd_dT, rmsd_type=args.rmsd_type, colorrange=args.rmsd_colorrange, \
                    display_type = args.rmsd_grid, file_name=args.rmsd_grid_file, rmsd_lambda = rmsd_lambda, \
                    colormap=colormap)

if __name__ == "__main__":
    main()
