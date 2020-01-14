
import pint

units = pint.UnitRegistry()

# from http://www.met.psu.edu/weather/other-resources/meteorological-measurements-units-and-conversions-1

# Universal Gas Constant
R = 8.3145e3 * units.J * units.K ** -1 * units.mol ** -1

# Average Molecular Weight of Dry Air
Md = 28.97 * units.g * units.mol ** -1

# Gas Constant of Dry air
Rd = 287 * units.J * units.K ** -1 * units.kg ** -2

# Density of Dry Air at 0°C and 1000mb
rd = 1.275 * units.kg * units.m ** -3



# Dimensions for indexing state
x_dim = -1
y_dim = -2
z_dim = -3


"""
Functions for using units with numpy helper functions
"""
def unit_roll(a, shift, axis=None):
    return np.roll(a, shift, axis=axis) * a.units


def unit_maximum(a, b):
    return np.maximum(a, b) * a.units


def unit_minimum(a, b):
    return np.minimum(a, b) * a.units
