import numpy as np
import unittest
from constants import units, standard_temperature, mu_air
from viscosity import *


class TestViscosity(unittest.TestCase):
    def test_laplacian(self):
        dx = 1 * units.meter
        p = np.full((5, 5), standard_temperature, dtype=np.float) * units.kelvin
        p[0] = 100 * units.kelvin + standard_temperature

        p_out = finite_laplacian_2d(p, dx)
        print(p_out.magnitude)

    def test_incompressible_part(self):
        dx = 300 * units.km
        u = np.zeros((5, 5)) * units.m / units.s
        u[2, 2] = 1 * units.m / units.s

        out = incompressible_viscosity_2d(u, mu_air, dx)
        print(out.magnitude)
        print(out.u)





if __name__ == '__main__':
    unittest.main()
