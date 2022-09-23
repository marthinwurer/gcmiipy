
import unittest
import constants
from constants import units


def to_true_temp(t, p):
    if hasattr(t, "shape"):
        assert t.shape == p.shape
    tt = t / ((constants.P0 / p) ** constants.kappa)
    # rho = p / (constants.Rd * tt)
    return tt


def to_potential_temp(tt, p):
    if hasattr(tt, "shape"):
        assert tt.shape == p.shape
    t = tt * ((constants.P0 / p) ** constants.kappa)
    return t


def to_density(tt, p):
    rho = p / (constants.Rd * tt)
    return rho





class TestTemperature(unittest.TestCase):
    def test_temperature_conversion(self):
        tt = constants.standard_temperature
        p = constants.standard_pressure

        t = to_potential_temp(tt, p)

        tt2 = to_true_temp(t, p)

        print(tt, t, tt2)

        self.assertAlmostEqual(tt, tt2)







