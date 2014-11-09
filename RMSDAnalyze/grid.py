import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import scipy.sparse as mtx
import numpy.linalg as LA
import logging
import coords
import op
from sklearn.neighbors import NearestNeighbors

class PlotLabeler:
    def __init__(self, title=None, xlabel=None, ylabel=None, colorrange=None, colormap=plt.cm.Spectral_r):
        self.title = title
        self.xlabel= xlabel
        self.ylabel= ylabel
        self.colorrange=colorrange
        self.colormap=colormap

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
        plot_out = plt.hexbin(x, y, C=value, 
                   cmap=plotlabeler.colormap, 
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
    logging.debug("PRE  nonzero shape: {}".format(nonzero.shape))
    nonzero.shape = (nonzero.shape[0], nonzero.shape[1])
    logging.debug("PRE  running_mean shape: {}".format(running_mean.shape))
    running_mean   = running_mean[nonzero, :] / running_weight[nonzero]
    logging.debug("POST running_mean shape: {}".format(running_mean.shape))
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
    return running_mean, running_weight, data_pos

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
    logging.debug("Number of objects being added: {}".format(len(x)))
    if len(value.shape) == 1:
        value.shape = (value.shape[0], 1)
    value_shape  = (2*gridsize[0]+1, 2*gridsize[1]+1, value.shape[1])
    weight_shape = (2*gridsize[0]+1, 2*gridsize[1]+1, 1)
    if running_mean_mtx is None:
        running_mean_mtx = np.zeros(value_shape)
    if running_weight_mtx is None:
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

def GridOP(data_tik, display_type=[], dynamic_step = 0, colorrange=[None,None],
        op_type='density', file_name=None, rmsd_lambda=None, 
        colormap=plt.cm.Spectral_r, water_pos=83674, ion_pos = 423427, 
        coord_system = coords.SlabCoords(10.0, 4.0), gridsize=[40,30],
        pbc = None, nframes = None, plot = True):

    if dynamic_step > 0 and op_type in op.static_op:
        raise ValueError("Cannot use a dynamic step {} > 0 with an op_type {}"
                         .format(dynamic_step, op_type))
    if dynamic_step == 0 and op_type in op.dynamic_op:
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
        if op_type in op.dynamic_op:
            atoms.append(data_tik[t0+dynamic_step,:,:])
        atoms, op_i = op.OPCompute(atoms, atom_type, op_type, 
                water_pos, ion_pos, rmsd_lambda, pbc)
        # Convert to cylindrical coords and center
        #center_k = np.mean(data_tik[t0,:,:], axis=0)
        #r_ik = atoms[0] - center_k
        r_ik = atoms[0]
        extent = coord_system.GetExtent()
        r,z,op_i = coord_system(r_ik, op_i)
        running_mean_mtx, running_weight_mtx = UpdateRunningMean2D(
                running_mean_mtx, running_weight_mtx, 
                r, z, op_i, extent, gridsize)
        print running_mean_mtx[running_mean_mtx != 0]

    data, weight, pos = ProcessSparseRunner(
            running_mean_mtx, running_weight_mtx, 
            coord_system.GetExtent(), gridsize)
    data_pos = pos

    logging.debug("Old density size: {}".format(weight.shape))
    density = coord_system.UndoJacobian(weight, pos[:, 0], pos[:, 1], gridsize)
    logging.debug("New density size: {}".format(density.shape))
    density_pos = pos
    
    if op_type == 'rmsd':
        logging.debug("Data shape: {}".format(data.shape))
        mean_r = data[:,[1,2,3]]

        data = np.sqrt(data[:,0] - np.sum(np.square(mean_r), axis=1))

    if plot:
        # PLOTTING FUNCTION -- should be separate function, but it shares too many arguments
        # Build the plotlabeler
        if rmsd_lambda:
            title=rmsd_lambda.title
        elif op_type=='q6':
            title='q6 plot'
        else:
            title='Some generic order parameter plot'
        plotlabeler=PlotLabeler(title=title, 
                xlabel="R, cylindrical radius from center of disc (nm)", 
                ylabel="Z, vertical height (nm)", 
                colormap = colormap, 
                colorrange = colorrange)
        # Plot OP image
        OPPlotter2D(data_pos[:,0],data_pos[:,1], data, 
                extent, gridsize, plotlabeler=plotlabeler, subplot=(2,1,1))
        # Plot the protein structure
        center_k = np.mean(data_tik[:,0:water_pos,:], axis=(0,1))
        protein_ik = data_tik[0,0:water_pos,:] - center_k[0]
        protein_r, protein_z = coord_system(protein_ik)
        # Plot density image
        plotlabeler.title = "Density, no units"
        logging.debug("density: {}".format(density))
        plotlabeler.colorrange = None
        OPPlotter2D(density_pos[:,0],density_pos[:,1], density, 
                extent, gridsize, plotlabeler=plotlabeler, subplot=(2,1,2))
        if 'png' in display_type:
            if file_name:
                plt.savefig(file_name + '.png')
            else:
                plt.savefig('default.png')
        if 'display' in display_type:
            plt.show()
    return data, density, pos
