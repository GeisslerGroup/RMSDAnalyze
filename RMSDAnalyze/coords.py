import numpy as np

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
    def UndoJacobian(self, value, dir1, dir2, gridsize):
        return value
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
    def GetExtent(self):
        return [0, self.r_extent, -self.z_extent, self.z_extent]
    def UndoJacobian(self, value, r, z, gridsize):
        extra_rad = self.r_extent / float(gridsize[0])
        logging.debug("value shape, r shape: {}, {}".format(value.shape, r.shape))
        return value[:,0] / (r + extra_rad)
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
