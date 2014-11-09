import grid
import logging
import numpy as np
import op
import coords

def LocalOP(data_tik, dynamic_step = 0, op_type='density', 
                   rmsd_lambda=None, water_pos=83674, ion_pos = 423427, 
                   coord_system = [coords.LocalSlabCoords(0.0, 1.0, 0.0, 1.0)], 
                   bins=np.linspace(0,1,100), pbc=None, nframes = None):
    atom_type = 'water'
    if not nframes:
        nframes = data_tik.shape[0] - dynamic_step
    else:
        nframes = min(data_tik.shape[0] - dynamic_step, nframes)

    if dynamic_step > 0 and op_type in op.static_op:
        raise ValueError("Cannot use a dynamic step {} > 0 with an op_type {}".
                format(dynamic_step, op_type))
    if dynamic_step == 0 and op_type in op.dynamic_op:
        raise ValueError("Cannot use a dynamic step == 0 with an op_type {}".
                format(dynamic_step, op_type))

    op_distribution = np.zeros(len(bins) - 1)
    if (op_type == 'rmsd'):
        mean_arr = []
        mean_wgt = []
        for t0 in xrange(nframes):
            atoms = [data_tik[t0,:,:]]
            if op_type in op.dynamic_op:
                atoms.append(data_tik[t0+dynamic_step,:,:])
            atoms, op_i = op.OPCompute(atoms, atom_type, op_type, 
                    water_pos, ion_pos, rmsd_lambda, pbc = None)
            _, _, op_i = coord_system(atoms[0], op_i)
            if len(op_i):
                # Remove the mean from the RMSD
                r_frame = op_i[:,[1,2,3]]
                mean_arr.append(np.mean(r_frame, axis = 0))
                mean_wgt.append(len(op_i))
        r_mean = np.average(mean_arr, weights=mean_wgt, axis=0)
    
    for t0 in xrange(nframes):
        atoms = [data_tik[t0,:,:]]
        if op_type in op.dynamic_op:
            atoms.append(data_tik[t0+dynamic_step,:,:])
        atoms, op_i = op.OPCompute(atoms, atom_type, op_type, 
                water_pos, ion_pos, rmsd_lambda, pbc = None)
        logging.debug("!!!Number of ops: {}".format(op_i.shape))
        _, _, op_i = coord_system(atoms[0], op_i)
        logging.debug("!!!Number of ops: {}".format(op_i.shape))
        if len(op_i > 0):
            if (op_type == 'rmsd'):
                # Remove the mean from the RMSD
                r_frame = op_i[:,[1,2,3]]
                rmsd = np.sqrt(np.sum(np.square(r_frame - r_mean), axis = 1))
                op_i = rmsd
            hist, _ = np.histogram(op_i, bins=bins)
            op_distribution[:] += hist[:]

    return op_distribution

