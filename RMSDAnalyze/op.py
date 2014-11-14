import numpy.linalg as LA
import numpy as np
import logging

dynamic_op=['rmsd','angle']
static_op=['q6','density']

# [22.34405, 18.91217,  9.91809,  0.00000,  0.00000, -10.91895,  0.00000, -0.00000, -0.00000] 
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
    def GetPBCCenter(self):
        return np.sum(self.h, axis=0) / 2.
    def __str__(self):
        return "Rows are box vectors: \n{}".format(self.h)

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
        dr_ik = r1_ik- r0_ik 
    magsq_dr_i = np.sum(np.square(dr_ik), axis=1)
    # H = np.histogram(magsq_dr_i, bins=100)
    # logging.debug("Histogram of velocities: {}".format(H))
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
    if not rmsd_lambda is None:
        return rmsd_lambda(theta)
    else:
        return theta

def ComputeQ6(r_ik, cutoff_A):
    raise NotImplementedError("q6 order parameter not implemented")
    nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(r_ik)
    return np.ones(atoms_i.shape)

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
