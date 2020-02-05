import unittest

from constants import standard_pressure, standard_temperature
from matsumo_temp import density_from, geopotential_from
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


class TestTemperature(unittest.TestCase):
    def test_density_units(self):
        side_len = 3
        t = np.full((side_len, side_len), standard_temperature.m) * standard_temperature.u
        p = np.zeros((side_len, side_len), dtype=np.float) * standard_pressure.u
        H = standard_pressure
        p[:] = H

        density = density_from(p, t)

        self.assertEqual((1 * units.kg / units.m ** 3).to_base_units().u, density.to_base_units().u)

    def test_geopotential_units(self):
        side_len = 3
        t = np.full((side_len, side_len), standard_temperature.m) * standard_temperature.u
        p = np.zeros((side_len, side_len), dtype=np.float) * standard_pressure.u
        H = standard_pressure
        p[:] = H

        density = density_from(p, t)
        geo = geopotential_from(density, p)

        self.assertEqual((1 * units.m).to_base_units().u, geo.to_base_units().u)
