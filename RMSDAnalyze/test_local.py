import matplotlib
matplotlib.use('Agg')
import unittest
import local
import coords
import op
import numpy as np


class TestLocal(unittest.TestCase):
    def setUp(self):
        data_tik = []
        data_tik.append( [[0.0, 0, 0]] )
        data_tik.append( [[0.5, 0, 0]] )
        data_tik.append( [[1.0, 0, 0]] )
        self.data_vel_tik = np.array( data_tik )
        data_tik.append( [[0.5, 0, 0]] )
        data_tik.append( [[0.0, 0, 0]] )
        self.data_osc_tik = np.array( data_tik )

    def test_local(self):
        cs = coords.SlabCoords(10,4)
        bins = np.linspace(0,1,100)
        op_dist = local.LocalOP(self.data_osc_tik, op_type = 'rmsd', 
                dynamic_step = 1, water_pos = 0, ion_pos = 1, 
                coord_system = cs, bins=bins)
        self.assertEqual(op_dist[0], 0)
        self.assertEqual(max(op_dist), 4)

        op_dist = local.LocalOP(self.data_vel_tik, op_type = 'rmsd', 
                dynamic_step = 1, water_pos = 0, ion_pos = 1, 
                coord_system = cs, bins=bins)
        self.assertEqual(op_dist[0], 2)


if __name__=="__main__":
    unittest.main()
