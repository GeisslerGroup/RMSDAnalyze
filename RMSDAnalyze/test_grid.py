import matplotlib
matplotlib.use('Agg')
import unittest
import grid
import coords
import op
import numpy as np


class TestGrid(unittest.TestCase):
    def setUp(self):
        data_t0_ik = np.array( [[0.0, 0, 0]] )
        data_t1_ik = np.array( [[0.1, 0, 0]] )
        data_t2_ik = np.array( [[0.2, 0, 0]] )
        self.data_tik = np.array([ data_t0_ik, 
                                   data_t1_ik, 
                                   data_t2_ik ])

    def test_grid(self):
        cs = coords.SlabCoords(10,4)
        data, density, pos = grid.GridOP(self.data_tik, op_type = 'rmsd', 
                dynamic_step = 1, plot = False, water_pos = 0, ion_pos = 1, 
                coord_system = cs, gridsize=[1, 1])
        self.assertAlmostEqual(data[0], 0.)


if __name__=="__main__":
    unittest.main()
