import unittest
import numpy as np

from no_limits_2_5d import gen_geometry, plot_callback, layers, height, width

from constants import *
import temperature

from phi_port import *



class TestPhiPort(unittest.TestCase):
    def test_types(self):
        geom = gen_geometry(height, width, layers)
        p = np.full((height, width), 1) * standard_pressure - geom.ptop
        u = np.full((layers, height, width), 1) * 1.0 * units.m / units.s
        v = np.full((layers, height, width), 1) * .0 * units.m / units.s
        q = np.full((layers, height, width), 1) * 0.1 * units.dimensionless
        t = np.full((layers, height, width), 1) * temperature.to_potential_temp(standard_temperature, p)
        dx = 100 * units.m
        dt = 60 * 15 * units.s

        pt = np.transpose(p)
        tt = np.transpose(t)

        phi_t = PGF(tt, pt, geom)

        phi = np.transpose(phi_t)



