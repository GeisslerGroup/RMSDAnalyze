import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import scipy.sparse as mtx
import numpy.linalg as LA
import logging
from sklearn.neighbors import NearestNeighbors

dynamic_op=['rmsd','angle']
static_op=['q6','density']

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

class SlabCoords:
    def __init__(self, dir1_extent, dir2_extent, dir1=0, dir2=2, thickness=1):
        """
        Select a region for all: abs(r[dir1]) < dir1_extent
                                 abs(r[dir2]) < dir2_extent
        Collapse the remaining direction for all atoms in a slab of dimension [thickness].
        """
        if (dir1 == dir2):
            raise ValueError("In SlabCoords: dir1 ({}) cannot be equal to dir2 ({}).".format(dir1, dir2))
        self.dir1     = dir1
        self.dir2     = dir2
        self.dir1_extent = dir1_extent
        self.dir2_extent = dir2_extent
        # dir_trunc computes the remaining direction
        self.dir_trunc= sum([0,1,2]) - dir1 - dir2 
        self.thickness = thickness
    def GetExtent(self):
        return [-self.dir1_extent, self.dir1_extent, -self.dir2_extent, self.dir2_extent]
    def __call__(self, r_ik, op_i=None):
        # Truncate to relevant regions of the box
        sub = ((np.abs(r_ik[:, self.dir1     ]) < self.dir1_extent) * 
               (np.abs(r_ik[:, self.dir2     ]) < self.dir2_extent) *
               (np.abs(r_ik[:, self.dir_trunc]) < self.thickness/2.0))
        r    = r_ik[sub,self.dir1]
        z    = r_ik[sub,self.dir2]
        if op_i != None:
            op_i = op_i[sub]
            return r, z, op_i
        else:
            return r, z

class RadialCoords:
    def __init__(self, r_extent, z_extent):
        self.r_extent = r_extent
        self.z_extent = z_extent
    def GetMtxScale(self, gridsize):
        return max(gridsize[0] / self.r_extent, gridsize[1] / (2 * self.z_extent)) * 100
    def GetExtent(self):
        return [0, self.r_extent, -self.z_extent, self.z_extent]
    def __call__(self, r_ik, op_i=None):
        r_cyl_i = np.sqrt( np.square(r_ik[:,0]) + np.square(r_ik[:,1]))
        z_cyl_i = r_ik[:,2]
        # Truncate to relevant regions of the box
        sub = (r_cyl_i < self.r_extent) * (np.abs(z_cyl_i) < self.z_extent)
        r    = r_cyl_i[sub]
        z    = z_cyl_i[sub]
        if op_i != None:
            op_i = op_i[sub]
            return r, z, op_i
        else:
            return r, z

# 22.34405  18.91217   9.91809   0.00000   0.00000 -10.91895   0.00000  -0.00000  -0.00000
class HexPBC:
    """ 
    Minimum image convention periodic boundary conditions for generic 
            triclinic systems. Constructor takes the box vectors that come
            from gromacs simulations and is output from the .gro file format.

    FROM: M. E. Tuckerman. Statistical Mechanics: Theory and Molecular 
            Simulation.  Oxford University Press, Oxford, UK, 2010

    h := [a1, a2, a3], where ai is each box vector.
    s_i = h^{-1} r_i  
    s_ij = s_i - s_j
    s_ij <-- s_ij - NINT(s_ij)        (general minimum image convention)
    r_ij = h s_ij
    """

    def __init__(self, gromacs_fmt):
        a1 = [ gromacs_fmt[0], gromacs_fmt[3], gromacs_fmt[4] ]
        a2 = [ gromacs_fmt[5], gromacs_fmt[1], gromacs_fmt[6] ]
        a3 = [ gromacs_fmt[7], gromacs_fmt[8], gromacs_fmt[2] ]
        self.h = np.vstack([a1, a2, a3])
        self.hinv = LA.inv(self.h)
    def NearestPeriodic(self, r1_ik, r2_ik):
        s1_ki = np.dot(self.hinv, r1_ik.T)
        s2_ki = np.dot(self.hinv, r2_ik.T)
        ds_ki = s1_ki - s2_ki
        ds_ki -= np.rint(ds_ki)
        dr_ik = np.dot(self.h, ds_ki).T
        return dr_ik
    def __str__(self):
        return "Rows are box vectors: \n{}".format(self.h)


def ComputeRMSD(r0_ik, r1_ik, pre_cutoff=None, post_cutoff=None, rmsd_lambda = None, pbc=None):
    """
    Compute the quantities that, when averaged, allow for a computation of RMSD.
    Returns OP_i = [|r|^2, dx, dy, dz]
    """
    if pbc:
        # Compute r1 - r0
        logging.debug("USING PERIODIC BOUNDARY CONDITIONS")
        dr_ik = pbc.NearestPeriodic(r1_ik, r0_ik)
    else:
        dr_ik = r0_ik- r1_ik 
    magsq_dr_i = np.sum(np.square(dr_ik), axis=1)
    H = np.histogram(magsq_dr_i, bins=100)
    logging.debug("Histogram of velocities: {}".format(H))
    magsq_dr_i.shape = (len(magsq_dr_i), 1)
    op_i = np.hstack([magsq_dr_i, dr_ik])
    return op_i



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
    raise NotImplementedError("q6 order parameter not implemented")
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
    plt.axis('equal')
    plt.tight_layout()

def GetMtxScale(extent, gridsize):
    L = (extent[1] - extent[0], extent[3] - extent[2])
    mtx_scale = np.array( [2.*float(gridsize[0]) / L[0],    
                           2.*float(gridsize[1]) / L[1]] )
    return mtx_scale

def ProcessSparseRunner(running_mean, running_weight, extent, gridsize):
    mtx_scale = GetMtxScale(extent, gridsize)
    nonzero = (running_weight != 0)
    running_mean   = running_mean[nonzero] / running_weight[nonzero]
    running_weight = running_weight[nonzero]

    logging.debug("Matrix scale: {}, shape = {}".format(mtx_scale, mtx_scale.shape))
    logging.debug("Grid size: {}".format(gridsize))
    x = np.arange(0, 2*gridsize[0] + 1, dtype=float) / mtx_scale[0]
    y = np.arange(0, 2*gridsize[1] + 1, dtype=float) / mtx_scale[1]
    x += extent[0]
    y += extent[2]
    yv, xv = np.meshgrid(y,x)
    nonzero.shape = (nonzero.shape[0], nonzero.shape[1])
    logging.debug("xv size, nonzero-array size: {}, {}".format(xv.shape, nonzero.shape))
    data_x = xv[nonzero].flatten()
    data_y = yv[nonzero].flatten()
    data_pos = np.vstack([data_x, data_y]).T
    logging.debug("Reconstituted data_pos: x={}, y={}".format(data_pos[:,0], data_pos[:,1]))
    logging.debug("data_pos size: {}".format( data_pos.shape))
    return (running_mean, data_pos), (running_weight, data_pos)

def CenterForMatrix(data, data_pos, extent, gridsize):
    mtx_scale = GetMtxScale(extent, gridsize)
    data_pos = mtx_scale * ( data_pos - np.array([extent[0], extent[2]]))
    data_pos = np.abs(np.round(data_pos))
    data_pos = data_pos.astype(int)

    # Include an index-column to the xy-data, and then sort on x,y
    index = np.arange(0, len(data_pos))
    index.shape = (len(data_pos), 1)
    merge = np.hstack([index,data_pos])
    merge = np.sort(merge.view('i8,i8,i8'), order=['f1','f2'])
    data = data[merge.view(np.int)[:,0]] 
    data_pos = merge.view(np.int)[:,(1,2)]
    logging.debug("Sorted data_pos: x={}, y={}".format(data_pos[:,0], data_pos[:,1]))
    return (data, data_pos)

def UpdateRunningMean2D(running_mean_mtx, running_weight_mtx, x, y, value, extent, gridsize, style='hex'):
    logging.debug("Number of entries: {}".format(x.shape))
    value_shape  = (2*gridsize[0]+1, 2*gridsize[1]+1, value.shape[1])
    weight_shape = (value_shape[0], value_shape[1], 1)
    if running_mean_mtx == None:
        running_mean_mtx = np.zeros(value_shape)
    if running_weight_mtx == None:
        running_weight_mtx = np.zeros(weight_shape)
    if style=='hex':
        mean_mtx   = np.zeros(value_shape)
        weight_mtx = np.zeros(weight_shape)
        
        logging.debug("Shape of mean_mtx: {}".format(mean_mtx.shape))
        # Compute the mean of all of the values in "value" 
        for k in xrange(value.shape[1]):
            meanplt   = plt.hexbin(x, y, C=value[:,k], mincnt=0, 
                    gridsize = gridsize, extent=extent); plt.clf()
            mean   = meanplt.get_array()
            mean, mean_pos = CenterForMatrix(mean, meanplt.get_offsets(), extent, gridsize)
            mean_mtx  [mean_pos[:,0], mean_pos[:,1], k] = mean[:]
        # Compute the statistical weight of those entries
        weightplt = plt.hexbin(x, y, mincnt=0, 
                gridsize = gridsize, extent=extent); plt.clf()
        weight = weightplt.get_array()
        weight, weight_pos = CenterForMatrix(weight, weightplt.get_offsets(), extent, gridsize)
        weight_mtx[weight_pos[:,0], weight_pos[:,1], 0] = weight[:]

        running_mean_mtx   += mean_mtx * weight_mtx
        running_weight_mtx += weight_mtx
        return running_mean_mtx, running_weight_mtx
    else:
        raise ValueError("Plotting style must be 'hex', received {}". format(style))


def OPCompute(atoms, atom_type, op_type, water_pos, ion_pos, rmsd_lambda, pbc=None):
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
        op_i = ComputeRMSD(atoms[0], atoms[1], post_cutoff=5.0, rmsd_lambda = rmsd_lambda, pbc=pbc)
    elif op_type == 'angle':
        op_i = ComputeRMSAngle(atoms[0], atoms[1], rmsd_lambda = rmsd_lambda)
    else:
        raise ValueError("GridOPRadial passed op_type that is not known: {}".format(op_type))
    # POST-FILTER
    if op_type == 'angle' and atom_type == 'water':
        atoms[0] = atoms[0][::3]
    # Convert to cylindrical coords and center
    return atoms, op_i


def GridOP(data_tik, display_type=[], dynamic_step = 0, colorrange=[None,None],
        op_type='density', file_name=None, rmsd_lambda=None, 
        colormap=plt.cm.Spectral_r, water_pos=83674, ion_pos = 423427, 
        coord_system = SlabCoords(10.0, 4.0), gridsize=[40,30],
        pbc = None, nframes = None):
    if dynamic_step > 0 and op_type in static_op:
        raise ValueError("Cannot use a dynamic step {} > 0 with an op_type {}"
                         .format(dynamic_step, op_type))
    if dynamic_step == 0 and op_type in dynamic_op:
        raise ValueError("Cannot use a dynamic step == 0 with an op_type {}"
                         .format(dynamic_step, op_type))

    if not nframes:
        nframes = data_tik.shape[0] - dynamic_step
    else:
        nframes = min(data_tik.shape[0] - dynamic_step, nframes)

    atom_type = 'water'
    running_mean_mtx = None
    running_weight_mtx = None
    for t0 in xrange(nframes):
        #print "outputting time {} of {}".format(t0, data_tik.shape[0])
        atoms = [data_tik[t0,:,:]]
        if op_type in dynamic_op:
            atoms.append(data_tik[t0+dynamic_step,:,:])
        atoms, op_i = OPCompute(atoms, atom_type, op_type, water_pos, ion_pos, rmsd_lambda, pbc)
        # Convert to cylindrical coords and center
        center_k = np.mean(data_tik[t0,:,:], axis=0)
        r_ik = atoms[0] - center_k
        extent = coord_system.GetExtent()
        r,z,op_i = coord_system(r_ik, op_i)
        running_mean_mtx, running_weight_mtx = UpdateRunningMean2D(running_mean_mtx, running_weight_mtx, r, z, op_i, \
                                                                    extent, gridsize)
    (data, data_pos), (density,density_pos) = \
            ProcessSparseRunner(running_mean_mtx, running_weight_mtx, 
                                coord_system.GetExtent(), gridsize)

    # PLOTTING FUNCTION -- should be separate function, but it shares too many arguments
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
    protein_ik = data_tik[0,0:water_pos,:] - center_k[0]
    protein_r, protein_z = coord_system(protein_ik)
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
