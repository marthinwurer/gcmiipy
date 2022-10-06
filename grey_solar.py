



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


def mmr_from_vmr(vmr, mmg, mma):
    """
    Using the molar mass of a gas and the molar mass of air,
    compute the mass mixing ratio from the volumetric mixing ratio
    """
    return vmr * mmg / mma


def solar_zenith_angle(latitude, hour_angle, declination):
    """
    from https://en.wikipedia.org/wiki/Solar_zenith_angle
    """
    return np.sin(latitude) * np.sin(declination) +\
            np.cos(latitude) * np.cos(declination) * np.cos(hour_angle)


def zenith_angle(longs, lats, time):
    # compute the hour angle
    hour_angle = time / (24 * units.hours) * 360 * units.degrees
    # tile the longs and lats to make the size


"""
Do a basic grey gas model with o3, co2, and h2o

Do single SW pass with ozone, clouds, and water, then both LW passes with CO2 and water.
"""

# from page 51 of atmospheric dynamics chapter 2
# also page 66
h2o_weight = 0.125 * units.m ** 2 * units.kg ** -1
liquid_weight = 5.0 * units.m ** 2 * units.kg ** -1
co2_weight = 0.3 * units.m ** 2 * units.kg ** -1

ozone_weight = 0.01 * h2o_weight.u


def compute_absorbance(gasses, rho, path_length):
    absorbance = np.zeros(rho.shape)
    for gas, coefficient in gasses:
        plg = path_length * gas
        absorbance += gas * rho * path_length * coefficient

    return absorbance


def hansen_cloud_thickness(tp, tt):
    """
    Compute the cloud optical thickness from formula 21 in hansen 1983
    """
    thickness = ((tp - 100 * units.hPa) * 0.0133 / units.hPa).to_base_units()
    thickness[tt < 258 * units.K] = 1/3
    thickness[thickness < 0] = 0
    return thickness




def grey_solar(p, q, t, c, gt, utc, dt, geom:Geom):
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

    downwelling[-1] = solar_constant * 0.25

    geopotential_depth = (dp / (rho * G)).to_base_units()
    path_length = geopotential_depth
    print("path length")
    print(path_length)

    # Section 4.2.1 in Principles of Planetary Climate
    # gas concentrations
    gasses = [
        (oc, ozone_weight),
        (q, h2o_weight),
    ]
    # absorbance = np.zeros(tp.shape)
    # for gas, coefficient in gasses:
        # plg = path_length * gas
        # print("plg", plg.u)
        # absorbance += gas * rho * path_length * coefficient

    absorbance = compute_absorbance(gasses, rho, path_length)

    transmittance = 10 ** -absorbance
    a_cloud = absorbance * 1.66  # From Manabe 
    t_cloud = 10 ** -a_cloud
    print("absorbance", absorbance)
    print("transmittance", transmittance)

    cloud_thickness = hansen_cloud_thickness(tp, tt)

    cloud_albedo = (1 - np.exp(-cloud_thickness)) * 0.7
    print(cloud_albedo)

    absorbed = np.zeros(q.shape) * solar_constant.u

    for layer in reversed(range(geom.layers)):
        # take the basic equation from manabe 64 21a
        previous = downwelling[layer + 1]
        trans_layer = transmittance[layer]
        t_cloud_layer = t_cloud[layer]
        albedo_layer = cloud_albedo[layer]
        absorbed_nc = (1 - c) * (previous * (1 - trans_layer))
        reflected = c * albedo_layer * previous
        absorbed_c = c * (1 - albedo_layer) * previous * (1 - t_cloud_layer)

        total_absorbed = absorbed_nc + absorbed_c
        transmitted = previous - total_absorbed - reflected

        downwelling[layer] = transmitted
        absorbed[layer] = total_absorbed

    print(downwelling)
    print("absorbed")
    print(absorbed)

    dT = absorbed / Cp / rho / geopotential_depth * dt
    print(dT.to_base_units())

    tt_n = tt + dT
    t_n = temperature.to_potential_temp(tt_n, tp)
    print(t_n - t)

    return t_n, downwelling







def grey_radiation(p, q, t, c, gt, utc, dt, geom:Geom):

    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    rho = tp / (constants.Rd * tt)
    dp = p * geom.dsig

    geopotential_depth = (dp / (rho * G)).to_base_units()
    path_length = geopotential_depth

    flux_shape = (geom.layers + 1, geom.height, geom.width)
    thermal_downwelling = np.zeros(flux_shape) * solar_constant.u
    solar_downwelling = np.zeros(flux_shape) * solar_constant.u

    # hour_angle = 

    
    # zenith_angle = 
    solar_downwelling[-1] = solar_constant * 0.25

    oc = ozone_at(tp)
    sw_gasses = [
        (oc, ozone_weight),
        (q, h2o_weight),
    ]
    sw_absorbance = compute_absorbance(sw_gasses, rho, path_length)

    lw_gasses = [
        (q, h2o_weight),
        # TODO CO2 distribution
        (oc, ozone_weight),
    ]
    lw_absorbance = compute_absorbance(sw_gasses, rho, path_length)

    # longwave just has emissivity at the same rate as absorbtion
    # I think for the cloud absorbance we can just convert the cloud thickness
    # to absorbance and then add them

    cloud_thickness = hansen_cloud_thickness(tp, tt)
    sw_cloud_albedo = (1 - np.exp(-cloud_thickness)) * 0.7

    lw_cloud_absorbance = cloud_thickness / np.log(10) + lw_absorbance

    lw_emissivity = 1 - 10 ** -lw_absorbance
    lw_cloud_emssivity = 1 - 10 ** -lw_cloud_absorbance

    emittance = sb_constant * tt ** 4 * lw_emissivity
    print("emittance")
    print(emittance)











