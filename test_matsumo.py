import unittest
from matsuno_c_grid import *


class TestCGrid(unittest.TestCase):
    def test_ipj(self):
        p = np.full((3, 3), 1, dtype=np.float) * units.m
        p[1, 1] = 0 * units.m

        p_ipj = ipj(p)
        self.assertEqual(0, p_ipj[1, 0].m)

    def test_ijp(self):
        p = np.full((3, 3), 1, dtype=np.float) * units.m
        p[1, 1] = 0 * units.m

        p_ipj = ijp(p)
        self.assertEqual(0, p_ipj[0, 1].m)


    def test_pgf_v(self):
        p = np.full((3, 3), 1, dtype=np.float) * units.m
        p[1, 1] = 0 * units.m

        grad = geopotential_gradient_v(p, 1 * units.m)
        self.assertEqual(G.m, grad[1, 1].m)

