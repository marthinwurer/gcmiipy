



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


"""
Do a basic grey gas model with o3, co2, and h2o
"""

ozone_weight = 0.01
co2_weight = 0.05
# from page 51 of atmospheric dynamics chapter 2
h2o_weight = 0.125 * units.m ** 2 * units.kg ** -1


def grey_solar(p, q, t, c, gt, utc, geom:Geom):
    """
    gt: ground temperature
    """


    tp = p * geom.sig + geom.ptop
    oc = ozone_at(tp)

    downwelling = np.zeros((geom.height, geom.width, geom.layers + 1)) * solar_constant.u

    downwelling[:, :, -1] = solar_constant * 0.5

    dp = p * geom.dsig

    for layer in reversed(range(geom.layers)):

        downwelling[:,:,layer] =












