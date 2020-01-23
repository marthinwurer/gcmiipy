
import pint
import numpy as np

units = pint.UnitRegistry()

# from http://www.met.psu.edu/weather/other-resources/meteorological-measurements-units-and-conversions-1

# Universal Gas Constant
R = 8.3145e3 * units.J * units.K ** -1 * units.mol ** -1

# Average Molecular Weight of Dry Air
Md = 28.97 * units.g * units.mol ** -1

# Gas Constant of Dry air
Rd = 287 * units.J * units.K ** -1 * units.kg ** -1

# Density of Dry Air at 0°C and 1000mb
rd = 1.275 * units.kg * units.m ** -3

# Specific heat of dry air
Cp = 1004 * units.J * units.K ** -1 * units.kg ** -1

# potential temperature term
kappa = Rd / Cp

# standard reference pressure
P0 = 100000 * units.Pa  # 1000 hPa

# standard pressure and temperature
standard_pressure = 101325 * units.Pa
standard_temperature = 273.16 * units.K

# gravity on earth
G = 9.8 * units.m / units.s ** 2


# Dimensions for indexing state
x_dim = -1
y_dim = -2
z_dim = -3


"""
Functions for using units with numpy helper functions
"""
def unit_roll(a, shift, axis=None):
    return np.roll(a.m, shift, axis=axis) * a.units


def unit_maximum(a, b):
    return np.maximum(a.m, b.m) * a.units


def unit_minimum(a, b):
    return np.minimum(a.m, b.m) * a.units


def unit_stack(a):
    return np.stack([item.m for item in a]) * a[0].units


def get_total_variation(q):
    q_p_1 = unit_roll(q, -1, 0)
    diff = q - q_p_1
    return np.sum(np.abs(diff.m)) * q.u
