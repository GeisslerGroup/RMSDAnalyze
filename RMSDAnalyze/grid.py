import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import scipy.sparse as mtx
import numpy.linalg as LA

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
