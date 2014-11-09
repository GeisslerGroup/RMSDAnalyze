import unittest
import coords
import op
import numpy as np

class TestCoords(unittest.TestCase):
    def setUp(self):
        self.radial = coords.RadialCoords(10,4)
        self.slab   = coords.SlabCoords(10,4)
        self.pbc    = op.HexPBC([22.34405, 18.91217,   9.91809,  
                                  0.00000,  0.00000, -10.91895,  
                                  0.00000, -0.00000,  -0.00000])

    def test_pbc(self):
        r1 = np.array( [[0, 0,  8.91809]] )
        r2 = np.array( [[0, 0,  0]] )
        dr = self.pbc.NearestPeriodic(r2, r1)
        self.assertAlmostEqual(dr[0,2], 1.0)

if __name__=="__main__":
    unittest.main()
