import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import scipy.sparse as mtx
import numpy.linalg as LA
from sklearn.neighbors import NearestNeighbors

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
                self.title='Mean angular activity, cutoff={}rad, k={}rad^-1, t_delay={}ps'.format(self.cutoff, self.sharpness, self.rmsd_delay)
            else:
                print "SETTING TITLE BAD"
                self.title='Mean angular displacement, t_delay={}ps'.format(self.rmsd_delay)
        elif rmsd_type == 'rmsd':
            if self.b_activity and self.b_scaletime:
                self.title='Mean activity (scaled), cutoff={}nm/ps^(.5), k={}ps^(.5)/nm, t_delay={}ps'.format(self.cutoff, self.sharpness, self.rmsd_delay)
            elif self.b_activity:
                self.title='Mean activity, cutoff={}nm, k={}nm^-1, t_delay={}ps'.format(self.cutoff, self.sharpness, self.rmsd_delay)
            elif self.b_scaletime:
                self.title='Scaled RMSD (root diffusion), t_delay={}ps'.format(self.rmsd_delay)
            else:
                self.title='Mean RMSD, t_delay={}ps'.format(self.rmsd_delay)

class PlotLabeler:
    def __init__(self, title=None, xlabel=None, ylabel=None, colorrange=None, colormap=plt.cm.Spectral_r):
        self.title = title
        self.xlabel= xlabel
        self.ylabel= ylabel
        self.colorrange=colorrange
        self.colormap=colormap

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

def ComputeQ6(r_ik, cutoff_A):
    nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(r_ik)
    return np.ones(atoms_i.shape)

def OPPlotter2D(x,y, value, extent, gridsize, plotlabeler=PlotLabeler, style='hex', subplot=(1,1,1)):
    plt.subplot( subplot[0], subplot[1], subplot[2] )
    print "Gridsize: {}, extent: {}".format(gridsize, extent)
    if plotlabeler.colorrange:
        plot_out = plt.hexbin(x, y, C=value, \
                   vmin=plotlabeler.colorrange[0], \
                   vmax=plotlabeler.colorrange[1], \
                   cmap=plotlabeler.colormap, \
                   gridsize=gridsize, extent=extent)
    else:
        print x.shape
        print y.shape
        print value.shape
        plot_out = plt.hexbin(x, y, C=value, \
                   cmap=plotlabeler.colormap, \
                   gridsize=gridsize, extent=extent)
    cb = plt.colorbar(plot_out, spacing='uniform',extend='max')
    plt.title (plotlabeler.title)
    plt.xlabel(plotlabeler.xlabel)
    plt.ylabel(plotlabeler.ylabel)
    plt.tight_layout()

def UpdateRunningMean2D(running_mean_mtx, running_weight_mtx, x, y, value, extent, gridsize, mtx_scale=100, style='hex'):
    if style=='hex':
    # Compute the mean RMSD per box and the weight of the box
        meanplt   = plt.hexbin(x, y, C=value, \
                mincnt=0, gridsize = gridsize, extent=extent) 
        weightplt = plt.hexbin(x, y, \
                mincnt=0, gridsize = gridsize, extent=extent) 
    else:
        raise ValueError("Plotting style must be 'hex', received {}". format(style))
    plt.clf()
    #plt.show()
    mean   = meanplt.get_array()
    weight = weightplt.get_array()
    #print "mean shape, weight shape: {}, {}".format(hex_mean.shape, hex_weight.shape)
    mean_pos   =  mtx_scale * (meanplt.get_offsets()   + np.array([0,extent[3]]))
    weight_pos =  mtx_scale * (weightplt.get_offsets() + np.array([0,extent[3]]))
    mean_pos.astype(int)
    weight_pos.astype(int)
    mean_mtx   = mtx.csr_matrix( (mean  , (mean_pos[:,0],   mean_pos[:,1]  )), dtype='f')
    weight_mtx = mtx.csr_matrix( (weight, (weight_pos[:,0], weight_pos[:,1])) )
    if running_mean_mtx == None and running_weight_mtx == None:
        running_mean_mtx   = mean_mtx.multiply(weight_mtx)
        running_weight_mtx = weight_mtx
    else:
        running_mean_mtx   = running_mean_mtx + mean_mtx.multiply(weight_mtx)
        running_weight_mtx = running_weight_mtx + weight_mtx
    return running_mean_mtx, running_weight_mtx

def ProcessSparseRunner(running_mean, running_weight, mtx_scale, z_extent, coord='cyl'):
    #running_mean = running_mean / running_weight
    running_mean   = running_mean.tocoo()
    running_weight = running_weight.tocoo()
    #print "RUNNING MEAN TYPE: {}".format(type(running_mean_mtx))
    #print "RUNNING MEAN DATA: {}".format(running_mean_mtx)
    data_pos  = np.vstack((running_mean.row,running_mean.col)).astype(float).T / mtx_scale
    data_pos -= np.array([0,z_extent])
    weight_pos = np.vstack((running_weight.row,running_weight.col)).astype(float).T / mtx_scale
    weight_pos -= np.array([0,z_extent])
    if coord == 'cyl':
        weight_mean = running_weight.data / weight_pos[:,0]
    return (running_mean.data, data_pos), (weight_mean, weight_pos)

class RadialCoords:
    def __init__(self, r_extent, z_extent):
        self.r_extent = r_extent
        self.z_extent = z_extent
    def GetMtxScale(self, gridsize):
        return max(gridsize[0] / self.r_extent, gridsize[1] / (2 * self.z_extent)) * 100
    def ProcessSparseRunner(self, running_mean, running_weight, mtx_scale):
        running_mean = running_mean / running_weight
        running_mean   = running_mean.tocoo()
        running_weight = running_weight.tocoo()
        #print "RUNNING MEAN TYPE: {}".format(type(running_mean_mtx))
        #print "RUNNING MEAN DATA: {}".format(running_mean_mtx)
        data_pos  = np.vstack((running_mean.row,running_mean.col)).astype(float).T / mtx_scale
        data_pos -= np.array([0,self.z_extent])
        weight_pos = np.vstack((running_weight.row,running_weight.col)).astype(float).T / mtx_scale
        weight_pos -= np.array([0,self.z_extent])
        weight_mean = running_weight.data / weight_pos[:,0]
        return (running_mean.data, data_pos), (weight_mean, weight_pos)
    def GetExtent(self):
        return [0, self.r_extent, -self.z_extent, self.z_extent]
    def __call__(self, r_ik, op_i):
        r_cyl_i = np.sqrt( np.square(r_ik[:,0]) + np.square(r_ik[:,1]))
        z_cyl_i = r_ik[:,2]
        # Truncate to relevant regions of the box
        sub = (r_cyl_i < self.r_extent) * (np.abs(z_cyl_i) < self.z_extent)
        print sub.shape, op_i.shape
        r    = r_cyl_i[sub]
        z    = z_cyl_i[sub]
        op_i = op_i[sub]
        return r_cyl_i, z_cyl_i, op_i


def GridOPRadial_v2(data_tik, display_type=[], dynamic_step = 0, colorrange=[None,None], op_type='density', \
                   file_name=None, rmsd_lambda=None, colormap=plt.cm.Spectral_r, \
                   water_pos=83674, ion_pos = 423427, coord_system = RadialCoords(10.0, 4.0), gridsize=[40,30]):
    dynamic=['rmsd','angle']
    static=['q6','density']
    if dynamic_step > 0 and op_type in static:
        raise ValueError("Cannot use a dynamic step {} > 0 with an op_type {}".format(dynamic_step, op_type))
    if dynamic_step == 0 and op_type in dynamic:
        raise ValueError("Cannot use a dynamic step == 0 with an op_type {}".format(dynamic_step, op_type))

    atom_type = 'water'
    running_mean_mtx = None
    running_weight_mtx = None
    mtx_scale = coord_system.GetMtxScale(gridsize)
    for t0 in xrange(data_tik.shape[0] - dynamic_step):
        print "outputting time {} of {}".format(t0, data_tik.shape[0])
        atoms = [data_tik[t0,:,:]]
        if op_type in dynamic:
            atoms.append(data_tik[t0+dynamic_step,:,:])
        # PRE-FILTER
        if atom_type == 'water':
            for i in xrange(len(atoms)):
                atoms[i] = atoms[i][water_pos:ion_pos:,:]
        else:
            raise ValueError("GridOPRadial passed atom_type that is not known: {}".format(atom_type))
        if op_type != 'angle' and atom_type =='water':
            for i in xrange(len(atoms)):
                atoms[i] = atoms[i][::3]
        # Run the appropriate computation on those atoms
        if   op_type == 'q6':
            op_i = ComputeQ6(atoms[0], cutoff_A=4.0)
        elif op_type == 'density':
            op_i = np.ones(atoms[0].shape)
        elif op_type == 'rmsd':
            op_i = ComputeRMSD(atoms[0], atoms[1], post_cutoff=5.0, rmsd_lambda = rmsd_lambda)
        elif op_type == 'angle':
            op_i = ComputeRMSAngle(atoms[0], atoms[1], rmsd_lambda = rmsd_lambda)
        else:
            raise ValueError("GridOPRadial passed op_type that is not known: {}".format(op_type))
        # POST-FILTER
        if op_type == 'angle' and atom_type == 'water':
            atoms[0] = atoms[0][::3]
        # Convert to cylindrical coords and center
        center_k = np.mean(data_tik[t0,:,:], axis=0)
        r_ik = atoms[0] - center_k
        r,z,op_i = coord_system(r_ik, op_i)
        extent = coord_system.GetExtent()
        running_mean_mtx, running_weight_mtx = UpdateRunningMean2D(running_mean_mtx, running_weight_mtx, r, z, op_i, \
                                                                    extent, gridsize, mtx_scale=mtx_scale)
    (data, data_pos), (density,density_pos) = coord_system.ProcessSparseRunner(running_mean_mtx, running_weight_mtx, mtx_scale)
    # Build the plotlabeler
    if rmsd_lambda:
        title=rmsd_lambda.title
    elif op_type=='q6':
        title='q6 plot'
    else:
        title='Some generic order parameter plot'
    plotlabeler=PlotLabeler(title=title, \
                            xlabel="R, cylindrical radius from center of disc (nm)", \
                            ylabel="Z, vertical height (nm)", \
                            colormap = colormap, \
                            colorrange = colorrange)
    # Plot OP image
    OPPlotter2D(data_pos[:,0],data_pos[:,1], data, \
                extent, gridsize, plotlabeler=plotlabeler, subplot=(2,1,1))
    # Plot the protein structure
    center_k = np.mean(data_tik[:,0:water_pos,:], axis=(0,1))
    protein_r = np.sqrt(np.square(data_tik[0,0:water_pos,0] - center_k[0]) + \
                        np.square(data_tik[0,0:water_pos,1] - center_k[1]) )
    protein_z = data_tik[0,0:water_pos,2] - center_k[2]
    sub = (protein_r < r_extent) * (np.abs(protein_z) < z_extent)
    # Plot density image
    plotlabeler.title = "Density, no units"
    plotlabeler.colorrange = None
    OPPlotter2D(density_pos[:,0],density_pos[:,1], density, \
                extent, gridsize, plotlabeler=plotlabeler, subplot=(2,1,2))
    if 'png' in display_type:
        if file_name:
            plt.savefig(file_name + '.png')
        else:
            plt.savefig('default.png')
    if 'display' in display_type:
        plt.show()




def GridOPRadial_v1(data_tik, display_type=[], dynamic_step = 0, colorrange=[None,None], op_type='density', \
                   file_name=None, rmsd_lambda=None, colormap=plt.cm.Spectral_r, \
                   water_pos=83674, ion_pos = 423427, r_extent = 10.0, z_extent = 4.0, gridsize=[40,30]):
    dynamic=['rmsd','angle']
    static=['q6','density']
    if dynamic_step > 0 and op_type in static:
        raise ValueError("Cannot use a dynamic step {} > 0 with an op_type {}".format(dynamic_step, op_type))
    if dynamic_step == 0 and op_type in dynamic:
        raise ValueError("Cannot use a dynamic step == 0 with an op_type {}".format(dynamic_step, op_type))

    atom_type = 'water'
    running_mean_mtx = None
    running_weight_mtx = None
    mtx_scale = max(gridsize[0] / r_extent, gridsize[1] / (2 * z_extent)) * 100
    for t0 in xrange(data_tik.shape[0] - dynamic_step):
        print "outputting time {} of {}".format(t0, data_tik.shape[0])
        atoms = [data_tik[t0,:,:]]
        if op_type in dynamic:
            atoms.append(data_tik[t0+dynamic_step,:,:])

        # PRE-FILTER
        if atom_type == 'water':
            for i in xrange(len(atoms)):
                atoms[i] = atoms[i][water_pos:ion_pos:,:]
        else:
            raise ValueError("GridOPRadial passed atom_type that is not known: {}".format(atom_type))
        if op_type != 'angle' and atom_type =='water':
            for i in xrange(len(atoms)):
                atoms[i] = atoms[i][::3]
        # Run the appropriate computation on those atoms
        if   op_type == 'q6':
            op_i = ComputeQ6(atoms[0], cutoff_A=4.0)
        elif op_type == 'density':
            op_i = np.ones(atoms[0].shape)
        elif op_type == 'rmsd':
            op_i = ComputeRMSD(atoms[0], atoms[1], post_cutoff=5.0, rmsd_lambda = rmsd_lambda)
        elif op_type == 'angle':
            op_i = ComputeRMSAngle(atoms[0], atoms[1], rmsd_lambda = rmsd_lambda)
        else:
            raise ValueError("GridOPRadial passed op_type that is not known: {}".format(op_type))
        # POST-FILTER
        if op_type == 'angle' and atom_type == 'water':
            atoms[0] = atoms[0][::3]

        # Convert to cylindrical coords and center
        center_k = np.mean(data_tik[t0,:,:], axis=0)
        r_ik = atoms[0] - center_k
        r_cyl_i = np.sqrt( np.square(r_ik[:,0]) + np.square(r_ik[:,1]))
        z_cyl_i = r_ik[:,2]
        # Truncate to relevant regions of the box
        extent = [0, r_extent, -z_extent, z_extent]
        sub = (r_cyl_i < r_extent) * (np.abs(z_cyl_i) < z_extent)
        print sub.shape, op_i.shape
        r    = r_cyl_i[sub]
        z    = z_cyl_i[sub]
        op_i = op_i[sub]
        running_mean_mtx, running_weight_mtx = UpdateRunningMean2D(running_mean_mtx, running_weight_mtx, r, z, op_i, \
                                                                    extent, gridsize, mtx_scale=mtx_scale)
    (data, data_pos), (density,density_pos) = ProcessSparseRunner(running_mean_mtx, running_weight_mtx, mtx_scale, z_extent, coord='cyl')
   

    # Build the plotlabeler
    if rmsd_lambda:
        title=rmsd_lambda.title
    elif op_type=='q6':
        title='q6 plot'
    else:
        title='Some generic order parameter plot'
    plotlabeler=PlotLabeler(title=title, \
                            xlabel="R, cylindrical radius from center of disc (nm)", \
                            ylabel="Z, vertical height (nm)", \
                            colormap = colormap, \
                            colorrange = colorrange)

    # Plot OP image
    OPPlotter2D(data_pos[:,0],data_pos[:,1], data, \
                extent, gridsize, plotlabeler=plotlabeler, subplot=(2,1,1))
    # Plot the protein structure
    center_k = np.mean(data_tik[:,0:water_pos,:], axis=(0,1))
    protein_r = np.sqrt(np.square(data_tik[0,0:water_pos,0] - center_k[0]) + \
                        np.square(data_tik[0,0:water_pos,1] - center_k[1]) )
    protein_z = data_tik[0,0:water_pos,2] - center_k[2]
    sub = (protein_r < r_extent) * (np.abs(protein_z) < z_extent)
    # Plot density image
    plotlabeler.title = "Density, no units"
    plotlabeler.colorrange = None
    OPPlotter2D(density_pos[:,0],density_pos[:,1], density, \
                extent, gridsize, plotlabeler=plotlabeler, subplot=(2,1,2))

    if 'png' in display_type:
        if file_name:
            plt.savefig(file_name + '.png')
        else:
            plt.savefig('default.png')

    if 'display' in display_type:
        plt.show()
    
