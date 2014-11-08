import grid
import logging
import numpy as np

class LocalSlabCoords:
    def __init__(self, r_0, r_f, z_0, z_f, thickness=1):
        self.extent = [r_0, r_f, z_0, z_f]
        self.thickness = thickness
    def GetExtent(self):
        return self.extent
    def UndoJacobian(self, value, dir1, dir2, gridsize):
        return value
    def __call__(self, r_ik, op_i=None):
        # Truncate to relevant regions of the box
        sub = ((r_ik[:, 0] > self.extent[0]) * 
               (r_ik[:, 0] < self.extent[1]) *
               (r_ik[:, 2] > self.extent[2]) *
               (r_ik[:, 2] < self.extent[3]) *
               (np.abs(r_ik[:, 1]) < self.thickness/2.0))
        r    = r_ik[sub,0]
        z    = r_ik[sub,2]
        if op_i.any():
            op_i = op_i[sub]
            return r, z, op_i
        else:
            return r, z


def LocalOP(data_tik, dynamic_step = 0, op_type='density', \
                   rmsd_lambda=None, water_pos=83674, ion_pos = 423427, 
                   coord_system = [grid.RadialCoords(10.0, 4.0)], 
                   bins=np.linspace(0,1,100), pbc=None, nframes = None):
    atom_type = 'water'
    if not nframes:
        nframes = data_tik.shape[0] - dynamic_step
    else:
        nframes = min(data_tik.shape[0] - dynamic_step, nframes)

    if dynamic_step > 0 and op_type in grid.static_op:
        raise ValueError("Cannot use a dynamic step {} > 0 with an op_type {}".
                format(dynamic_step, op_type))
    if dynamic_step == 0 and op_type in grid.dynamic_op:
        raise ValueError("Cannot use a dynamic step == 0 with an op_type {}".
                format(dynamic_step, op_type))

    op_distribution = np.zeros(len(bins) - 1)
    for t0 in xrange(nframes):
        center_k = np.mean(data_tik[t0,:,:], axis=0)
        atoms = [data_tik[t0,:,:] - center_k]
        if op_type in grid.dynamic_op:
            atoms.append(data_tik[t0+dynamic_step,:,:] - center_k)
        atoms, op_i = grid.OPCompute(atoms, atom_type, op_type, 
                water_pos, ion_pos, rmsd_lambda, pbc = None)
        logging.debug("!!!Number of ops: {}".format(op_i.shape))
        _, _, op_i = coord_system(atoms[0], op_i[:,0])
        logging.debug("!!!Number of ops: {}".format(op_i.shape))
        if len(op_i > 0):
            hist, _ = np.histogram(op_i, bins=bins)
            op_distribution[:] += hist[:]

    return op_distribution

