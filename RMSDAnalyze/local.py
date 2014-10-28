import grid
import numpy as np

class LocalSlabCoords:
    def __init__(self, r_0, r_f, z_0, z_f, thickness=1):
        self.extent = [r_0, r_f, z_0, z_f]
        self.thickness = thickness
    def GetMtxScale(self, gridsize):
        return max(gridsize[0] / self.r_extent, gridsize[1] / (2 * self.z_extent)) * 100
    def ProcessSparseRunner(self, running_mean, running_weight, mtx_scale):
        extent = self.GetExtent()
        running_mean   = mtx.csr_matrix(running_mean/running_weight)
        running_mean   = running_mean.tocoo()
        running_weight = running_weight.tocoo()
        #print "RUNNING MEAN TYPE: {}".format(type(running_mean_mtx))
        #print "RUNNING MEAN DATA: {}".format(running_mean_mtx)
        data_pos  = np.vstack((running_mean.row,running_mean.col)).astype(float).T / mtx_scale
        data_pos += np.array([extent[0],extent[2]])
        weight_pos = np.vstack((running_weight.row,running_weight.col)).astype(float).T / mtx_scale
        weight_pos += np.array([extent[0],extent[2]])
        weight_mean = running_weight.data
        return (running_mean.data, data_pos), (weight_mean, weight_pos)
    def GetExtent(self):
        return self.extent
    def __call__(self, r_ik, op_i=None):
        # Truncate to relevant regions of the box
        sub = ((np.abs(r_ik[:, 0]) > self.extent[0]) * 
               (np.abs(r_ik[:, 0]) < self.extent[1]) *
               (np.abs(r_ik[:, 2]) > self.extent[2]) *
               (np.abs(r_ik[:, 2]) < self.extent[3]) *
               (np.abs(r_ik[:, 1]) < self.thickness/2.0))
        r    = r_ik[sub,0]
        z    = r_ik[sub,2]
        if op_i != None:
            op_i = op_i[sub]
            return r, z, op_i
        else:
            return r, z


def LocalOP(data_tik, dynamic_step = 0, op_type='density', \
                   rmsd_lambda=None, water_pos=83674, ion_pos = 423427, 
                   coord_system = [grid.RadialCoords(10.0, 4.0)], bins=np.linspace(0,1,100)):
    atom_type = 'water'

    if dynamic_step > 0 and op_type in grid.static_op:
        raise ValueError("Cannot use a dynamic step {} > 0 with an op_type {}".format(dynamic_step, op_type))
    if dynamic_step == 0 and op_type in grid.dynamic_op:
        raise ValueError("Cannot use a dynamic step == 0 with an op_type {}".format(dynamic_step, op_type))


    op_distribution = np.zeros(len(bins) - 1)
    for t0 in xrange(data_tik.shape[0] - dynamic_step):
        atoms = [data_tik[t0,:,:]]
        if op_type in grid.dynamic_op:
            atoms.append(data_tik[t0+dynamic_step,:,:])
        atoms, op_i = grid.OPCompute(atoms, atom_type, op_type, water_pos, ion_pos, rmsd_lambda)
        hist, _ = np.histogram(op_i, bins=bins)
        op_distribution[:] += hist[:]

    return op_distribution

