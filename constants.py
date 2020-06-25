
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

# non-standard reference pressure used by GCMII to simplify calculations
p_0_k = (1 * units.Pa) ** kappa

# standard pressure and temperature
standard_pressure = 101325 * units.Pa
standard_temperature = 273.16 * units.K

# temperature of mesopause
t_mesopause = 130 * units.K
p_mesopause = 0.0005 * units.kPa

# gravity on earth
G = 9.8 * units.m / units.s ** 2

# radius of earth
radius = 6.3781e6 * units.m

# dynamic viscosity (mu) of dry air at STP
mu_air = 18.5 * units.uPa * units.s

# Dimensions for indexing state
x_dim = -1
y_dim = -2
z_dim = -3


"""
Functions for using units with numpy helper functions
"""
def unit_roll(a, shift, axis=None):
    try:
        return np.roll(a.m, shift, axis=axis) * a.units
    except AttributeError:
        return np.roll(a, shift, axis=axis)



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


def courant_number(p, u, dx, dt):
    return ((np.max(u) + np.sqrt(np.mean(p) * G)) * dt / dx).to_base_units()


