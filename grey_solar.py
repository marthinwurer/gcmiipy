



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


co2_mmr = mmr_from_vmr(300 / 1e6, M_CO2, Md)


def daily_average_irradiance(lat, declination):
    dH = np.arccos(-np.tan(lat) * np.tan(declination))
    manabe64_Sc = 2 * 41840 * units.J * units.m ** -2 * units.minute ** -1
    return manabe64_Sc / (math.pi * units.radian) * (dH * np.sin(lat) * np.sin(declination) + 
            np.cos(lat) * np.cos(declination) * np.sin(dH))



def solar_zenith_angle(latitude, hour_angle, declination):
    """
    from https://en.wikipedia.org/wiki/Solar_zenith_angle
    returns cosine(zenith_angle)
    """
    return np.sin(latitude) * np.sin(declination) +\
            np.cos(latitude) * np.cos(declination) * np.cos(hour_angle)


def zenith_angle(longs, lats, time, geom):
    # compute the hour angle (negative because sun moves west)
    hour_angle = time / (-24 * units.hours) * 360 * units.degrees
    # print(hour_angle)
    # tile the longs and lats to make the size
    t_longs = np.tile(longs, (geom.height, 1))
    # print(t_longs)
    # print(lats)
    t_lats = np.tile(lats, (1, geom.width))
    # print(t_lats)
    point_angle = t_longs + hour_angle
    # print(point_angle)
    sza = np.maximum(solar_zenith_angle(lats, point_angle, 0 * units.degrees), 0)
    # print(sza)

    return sza


"""
Do a basic grey gas model with o3, co2, and h2o

Do single SW pass with ozone, clouds, and water, then both LW passes with CO2 and water.
"""

# from page 51 of atmospheric dynamics chapter 2
# also page 66
h2o_weight = 0.125 * units.m ** 2 * units.kg ** -1
liquid_weight = 5.0 * units.m ** 2 * units.kg ** -1
co2_weight = 1 * units.m ** 2 * units.kg ** -1
co2_sw_weight = 1 * units.m ** 2 * units.kg ** -1
co2_sw_weight = co2_weight

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







def grey_radiation(p, q, tt, c, g, utc, dt, geom:Geom):

    tp = p * geom.sig + geom.ptop
    rho = tp / (constants.Rd * tt)
    dp = p * geom.dsig

    geopotential_depth = (dp / (rho * G)).to_base_units()
    path_length = geopotential_depth

    flux_shape = (geom.layers + 1, geom.height, geom.width)
    thermal_downwelling = np.zeros(flux_shape) * solar_constant.u
    thermal_upwelling = np.zeros(flux_shape) * solar_constant.u
    solar_downwelling = np.zeros(flux_shape) * solar_constant.u

    # hour_angle = 
    irradiance = daily_average_irradiance(35 * units.degrees, 0 * units.degrees)
    irradiance = manabe64_Sc = 2 * 41840 * units.J * units.m ** -2 * units.minute ** -1
    irradiance *= 0.5 * 0.5
    
    # zenith_angle = 
    solar_downwelling[-1] = irradiance

    oc = ozone_at(tp)
    sw_gasses = [
        # (oc, ozone_weight),
        (q, h2o_weight),
        (co2_mmr, co2_sw_weight),
    ]
    sw_absorbance = compute_absorbance(sw_gasses, rho, path_length)
    sw_transmittance = 10 ** -sw_absorbance
    a_cloud = sw_absorbance * 1.66  # From Manabe 
    sw_t_cloud = 10 ** -a_cloud

    lw_gasses = [
        (q, h2o_weight),
        # TODO CO2 distribution
        (co2_mmr, co2_weight),
    ]
    lw_absorbance = compute_absorbance(lw_gasses, rho, path_length)
    # co2_absorbance = co2_mmr * rho * path_length * co2_weight
    # print("co2_absorbance", co2_absorbance[0], 1-10**-co2_absorbance[0])

    # longwave just has emissivity at the same rate as absorbtion
    # I think for the cloud absorbance we can just convert the cloud thickness
    # to absorbance and then add them

    cloud_thickness = hansen_cloud_thickness(tp, tt)
    sw_cloud_albedo = (1 - np.exp(-cloud_thickness)) * 0.7

    lw_cloud_absorbance = cloud_thickness / np.log(10) + lw_absorbance

    lw_emissivity = 1 - 10 ** -lw_absorbance
    lw_cloud_emissivity = 1 - 10 ** -lw_cloud_absorbance
    print("lw_emissivity", lw_emissivity[10])

    emittance = sb_constant * tt ** 4 * ((1 - c) * lw_emissivity + c * lw_cloud_emissivity)
    # print("emittance")
    # print(emittance)
    # treat the ground as fully black
    ground_emittance = sb_constant * g.gt ** 4

    # do the downwelling
    absorbed = np.zeros(q.shape) * solar_constant.u
    reflected = np.zeros(p.shape) * solar_constant.u
    for layer in reversed(range(geom.layers)):

        # Solar Downwelling
        # take the basic equation from manabe 64 21a
        previous = solar_downwelling[layer + 1]
        trans_layer = sw_transmittance[layer]
        t_cloud_layer = sw_t_cloud[layer]
        albedo_layer = sw_cloud_albedo[layer]
        absorbed_nc = (1 - c) * (previous * (1 - trans_layer))
        sw_reflected = c * albedo_layer * previous
        absorbed_c = c * (1 - albedo_layer) * previous * (1 - t_cloud_layer)

        total_absorbed = absorbed_nc + absorbed_c
        transmitted = previous - total_absorbed - sw_reflected
        reflected += sw_reflected

        solar_downwelling[layer] = transmitted
        absorbed[layer] += total_absorbed

        # LW downwelling
        previous = thermal_downwelling[layer + 1]

        cloud_absorbtion = c * previous * lw_cloud_emissivity[layer]
        clear_absorbtion = (1 - c) * previous * lw_emissivity[layer]
        total_absorbtion = cloud_absorbtion + clear_absorbtion
        lw_transmitted = previous - total_absorbtion

        absorbed[layer] += total_absorbtion
        thermal_downwelling[layer] = lw_transmitted + emittance[layer]

    # print(absorbed)
    # print(thermal_downwelling + solar_downwelling)

    # do ground stuff
    ground_albedo = 0.1
    ground_sw_absorbtion = (1 - ground_albedo) * solar_downwelling[0]
    ground_lw_absorbtion = thermal_downwelling[0]
    ground_absorbtion = ground_sw_absorbtion + ground_lw_absorbtion

    thermal_upwelling[0] = ground_emittance

    # do upwelling
    for layer in range(geom.layers):
        previous = thermal_upwelling[layer]

        cloud_absorbtion = c * previous * lw_cloud_emissivity[layer]
        clear_absorbtion = (1 - c) * previous * lw_emissivity[layer]
        total_absorbtion = cloud_absorbtion + clear_absorbtion
        lw_transmitted = previous - total_absorbtion

        absorbed[layer] += total_absorbtion
        thermal_upwelling[layer + 1] = lw_transmitted + emittance[layer]
    
    dt_ground = (ground_absorbtion - ground_emittance) / Cg / (.1 * units.m)
    # print("ground t:", g.gt)
    print("ground dt:", dt_ground.to_base_units())
    print(ground_absorbtion, ground_emittance)

    dt_air = (absorbed - 2 * emittance) / (Cp * rho * geopotential_depth)
    # print(dt_air.to_base_units().to(units.K / units.hour))
    # print(tt)
    print(solar_downwelling[-1], thermal_upwelling[-1],
            solar_downwelling[-1] - thermal_upwelling[-1])

    return dt_ground, dt_air, thermal_upwelling[-1]


def basic_grey_transmittances(t_lw, t_sw, geom):
    # 2.35
    e_n = 1 - t_lw ** (geom.dsig)
    e_n_sw = 1 - t_sw ** (geom.dsig)
    # print(geom.dsig)
    # exit()

    lw_transmittance = 1 - e_n
    sw_transmittance = 1 - e_n_sw

    return lw_transmittance, sw_transmittance


def basic_3_gas_absorbance(p, tp, tt, rho, q, geom):
    dp = p * geom.dsig
    geopotential_depth = (dp / (rho * G)).to_base_units()
    path_length = geopotential_depth

    oc = ozone_at(tp)
    sw_gasses = [
        # (oc, ozone_weight),
        # (q, h2o_weight), # page 51 in ad 2 says that this is just longwave
    ]
    sw_absorbance = compute_absorbance(sw_gasses, rho, path_length)

    lw_gasses = [
        (q, h2o_weight),
        # TODO CO2 distribution
        (co2_mmr, co2_weight),
    ]
    lw_absorbance = compute_absorbance(lw_gasses, rho, path_length)

    return lw_absorbance, sw_absorbance


def basic_grey_radiation(p, tp, tt, g, t_lw, t_sw, albedo, utc, geom):
    """
    implements the basic grey atmosphere from AD 2.7
    needs dsig to be constant
    """
    # 2.35
    e_n = 1 - t_lw ** (1 / geom.layers)
    e_n_sw = 1 - t_sw ** (1 / geom.layers)
    # print(e_n)
    # print("tt", tt)
    lw_transmittance, sw_transmittance = basic_grey_transmittances(t_lw, t_sw, geom)

    # 1) radiation emitted by each layer that reaches the surface
    # equation 2.25
    # lw_transmittance = np.full(tt.shape, 1 - e_n)
    # sw_transmittance = np.full(tt.shape, 1 - e_n_sw)
    emission = (1 - lw_transmittance) * sb_constant * tt ** 4
    # print(emission)
    cum_sw_trans_from_top = np.cumprod(sw_transmittance[::-1], axis=0)[::-1]
    cum_lw_trans_from_top = np.cumprod(lw_transmittance[::-1], axis=0)[::-1]
    cum_lw_trans_from_bottom = np.cumprod(lw_transmittance, axis=0)
    # print(cum_lw_trans_from_bottom)
    # print(cum_lw_trans_from_top)

    clw_b_div = cum_lw_trans_from_bottom / lw_transmittance

    # divide by transmittance to eliminate the current level's transmittance
    # from the cumprod, becuase we only want the levels below
    B = np.sum(emission * clw_b_div, axis=0)
    # print("B", B)

    # 2) radiation received from the sun
    # equation 2.26
    sza = zenith_angle(geom.long, geom.lat, utc, geom)
    # Sc = 342 * units.watt * units.m ** -2
    Sc = solar_constant * sza
    S = (1 - albedo) * Sc * cum_sw_trans_from_top[0]

    # 3) radiation emitted by the earth's surface
    # equation 2.27
    e_g = 1
    U_s = e_g * sb_constant * g.gt ** 4
    # print("U_s", U_s)

    # print("gt:", g.gt)
    dt_ground = (B + S - U_s) / Cg / (.1 * units.m)
    # print("dt_ground", dt_ground)
    # or do we want to do the balance form of the equation to speed convergence?
    # 2.28
    # nah this seems to make the sim explode
    # U_s = B + S
    # gt = (U_s / (e_g * sb_constant)) ** (1/4)
    # print(U_s)
    # print(gt)

    # now for the atomosphere

    # long wave from above
    LWA_a = None
    # ok, can I make a matrix that gives the transmission from a layer to a layer?
    """
    something like
      > to
    \/ from
    [[i1j1, i1j2, i1j2],
     [i2j1, i2j2, i2j2],
     [i3j1, i3j2, i3j2]]

    Have cumlw as this and its inverse:
    [[[0.1       ]]
     [[0.17782794]]
     [[0.31622777]]
     [[0.56234133]]]

    cumlw is:
    up: [t1, t1t2, t1t2t3, t1t2t3t4]
    dn: [t1t2t3t4, t2t3t4, t3t4, t4]
    so we want up / dn_rolled?
    dn_rolled: [0, t1t2t3t4, t2t3t4, t3t4]
    nah
    up_div: [1, t1, t1t2, t1t2t3]
    up_div_rolled: [0, 1, t1, t1t2]
    up_rolled: [0, t1, t1t2, t1t2t3]
    up_rolled_div_T: [0, 1, t2, t2t3], [0, t1/t2, t1, t1t3]
    [1, 0, 1, t3]
    dn_rolled_dn: [t2t3t4, t3t4, t4, 0]
    up_div_up_T: [1, t2, t2t3, t2t3t4], [t1/t2, 1, t3, t1t3ta4]
    dn_div_up_T: [t2t3t4, t2t3t4/t1, t3t4/t1, t4/t1
    dn_div_up_div_T: [t1t2t3t4, t3t4/t1, t4/t1t2, 1/t1t2t3
    
    up: [t1, t1t2, t1t2t3, t1t2t3t4]
    dn: [t1t2t3t4, t2t3t4, t3t4, t4]
    up_div_up_T_div_t: [1/t1, 1, t2, t2t3], [
    upz: [0, t1, t1t2, t1t2t3, t1t2t3t4]
    upz_div_upT: [0, 1, t2, t2t3], [0, 1/t2, 1, t3]




    [[0, 1, t2, t2t3],
     [1, 0, 1, t3],
     [t2, 1, 0, 1],
     [t2t3, t3, 1, 0],
     
    """

    flux_shape = (geom.layers + 1, geom.height, geom.width)
    upwelling = np.zeros(flux_shape) * Sc.u
    downwelling = np.zeros(flux_shape) * Sc.u
    downwelling_sw = np.zeros(flux_shape) * Sc.u
    # LWA_a = np.zeros(tt.shape) * Sc.u
    # for f in range(1, geom.layers):
        # for to in range(f):
            # cum_trans = np.full(tt[0].shape, 1.0)
            # for i in range(to+1, f):
                # cum_trans *= lw_transmittance[i]
# 
            # LWA_a[to] += emission[f] * (1 - lw_transmittance[to]) * cum_trans
            # LWA_a[to] += emission[f] * (1 - lw_transmittance[to]) * cum_lw_trans_from_top[to] / lw_transmittance[to] / cum_lw_trans_from_top[to]
            # LWA_a[to] += emission[f] * (1 - lw_transmittance[to]) * cum_lw_trans_from_bottom[to] / lw_transmittance[to] / cum_lw_trans_from_bottom[to]


    downwelling_sw[-1] = Sc
    absorbed_dw = np.zeros(tt.shape) * Sc.u
    absorbed_sw = np.zeros(tt.shape) * Sc.u
    for i in reversed(range(geom.layers)):
        # longwave
        absorbed_dw[i] = downwelling[i+1] * (1 - lw_transmittance[i])
        downwelling[i] = downwelling[i+1] * lw_transmittance[i] + emission[i]

        # shortwave
        absorbed_sw[i] = downwelling_sw[i+1] * (1 - sw_transmittance[i])
        downwelling_sw[i] = downwelling_sw[i+1] * sw_transmittance[i]

    LWA_a = absorbed_dw

    # print("LWA_a", LWA_a)
    # print("LWA_b", LWA_b)
    # print("serial", absorbed)
    # print("B", B)
    # print("lowest", downwelling[0])
    # print("S", S)
    # print("dw S", downwelling_sw[0] * (1 - albedo))

    # longwave from below
    # LWA_b = np.zeros(tt.shape) * Sc.u
    # for f in range(geom.layers):
        # for to in range(f+1, geom.layers):
            # cum_trans = np.full(tt[0].shape, 1.0)
            # for i in range(f+1, to):
                # cum_trans *= lw_transmittance[i]
            # LWA_b[to] += emission[f] * (1 - lw_transmittance[to]) * cum_trans
            # LWA_b[to] += emission[f] * (1 - lw_transmittance[to]) * cum_lw_trans_from_top[to] / lw_transmittance[to] / cum_lw_trans_from_top[to]
            # LWA_b[to] += emission[f] * (1 - lw_transmittance[to]) * cum_lw_trans_from_bottom[to] / lw_transmittance[to] / cum_lw_trans_from_bottom[to]

    absorbed = np.zeros(tt.shape) * Sc.u
    for i in range(geom.layers):
        absorbed[i] = upwelling[i] * (1 - lw_transmittance[i])
        upwelling[i + 1] = upwelling[i] * lw_transmittance[i] + emission[i]

    LWA_b = absorbed


    upwelling[0] = U_s
    absorbed_uw = np.zeros(tt.shape) * Sc.u
    for i in range(geom.layers):
        # longwave
        absorbed_uw[i] = upwelling[i] * (1 - lw_transmittance[i])
        upwelling[i + 1] = upwelling[i] * lw_transmittance[i] + emission[i]

    # print("LWA_a", LWA_a)
    # print("LWA_b", LWA_b)
    # print("serial", absorbed)

    # absorbed terrestrial radiation
    # 2.30
    U_n = clw_b_div * U_s * (1 - lw_transmittance)
    # print("U_n plus", U_n + LWA_b)
    # print("both absorbed", absorbed_uw)

    # absorbed solar radiation
    # 2.31
    S_n = (1 - sw_transmittance) * cum_sw_trans_from_top / sw_transmittance * Sc
    # print("S_n plus", S_n)
    # print("absorbed_sw", absorbed_sw)
    # exit()

    # emitted longwave radiation
    # 2.32
    B_n = emission
    # print("2B_n", 2*B_n)

    # print(LWA_a.u, LWA_b.u, U_n.u, S_n.u, B_n.u)

    # final temperature change
    # 2.34
    dTdt = (U_n + S_n - 2 * B_n + LWA_a + LWA_b) * (G / (Cp * p * geom.dsig))

    # print(dTdt.to_base_units())
    # assert g.gt > 0 * units.K
    # assert g.gt < 600 * units.K
    if(np.isnan(dTdt).any()):
        print("nanananan")
        exit()
    # exit()
    return dTdt, dt_ground









        












