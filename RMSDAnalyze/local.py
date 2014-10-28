import grid


def LocalOP(data_tik, dynamic_step = 0, op_type='density', \
                   rmsd_lambda=None, water_pos=83674, ion_pos = 423427, 
                   coord_system = [RadialCoords(10.0, 4.0)], gridsize=[40,30], nbins):
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
    dist = np.zeros(len(coord_system),
                    data_tik.shape[0] - dynamic_step,
                    n_bins)
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
        for i,cs in enumerate(coord_system):
            r,z,op_i = coord_system(r_ik, op_i)
            #dist[i, t0, :] = np.hist(
