



import math
import unittest

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

import constants
import low_pass
from coordinates_3d import *
from constants import *
import temperature
from ozone import ozone_at
from geometry import Geom


def solar_zenith_angle(latitude, hour_angle, declination):
    """
    from https://en.wikipedia.org/wiki/Solar_zenith_angle
    """
    return np.sin(latitude) * np.sin(declination) +\
            np.cos(latitude) * np.cos(declination) * np.cos(hour_angle)


"""
Do a basic grey gas model with o3, co2, and h2o

Do single SW pass with ozone, clouds, and water, then both LW passes with CO2 and water.
"""

# from page 51 of atmospheric dynamics chapter 2
h2o_weight = 0.125 * units.m ** 2 * units.kg ** -1

ozone_weight = 0.01 * h2o_weight.u
co2_weight = 0.05


def grey_solar(p, q, t, c, gt, utc, geom:Geom):
    """
    c: cloud fraction
    gt: ground temperature
    """


    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    rho = tp / (constants.Rd * tt)
    print("rho", rho.to_base_units())

    dp = p * geom.dsig
    oc = ozone_at(tp)

    downwelling = np.zeros((geom.layers + 1, geom.height, geom.width)) * solar_constant.u

    downwelling[-1] = solar_constant * 0.5

    dp = p * geom.dsig
    path_length = (dp / (rho * G)).to_base_units()
    print("path length")
    print(path_length)

    # Section 4.2.1 in Principles of Planetary Climate
    # gas concentrations
    gasses = [
        (oc, ozone_weight),
        (q, h2o_weight),
    ]
    absorbance = np.zeros(tp.shape)
    for gas, coefficient in gasses:
        plg = path_length * gas
        print("plg", plg.u)
        absorbance += gas * rho * path_length * coefficient

    transmittance = 10 ** -absorbance
    a_cloud = absorbance * 1.66
    t_cloud = 10 ** -a_cloud
    print("absorbance", absorbance)
    print("transmittance", transmittance)

    for layer in reversed(range(geom.layers)):
        # take the basic equation from manabe 64 21a
        previous = downwelling[layer + 1]
        trans_layer = transmittance[layer]
        t_cloud_layer = t_cloud[layer]
        non_cloud = (1 - c) * (previous * trans_layer)
        cloud = c * 0.5 * (previous * t_cloud_layer)

        downwelling[layer] = non_cloud + cloud

    print(downwelling)












