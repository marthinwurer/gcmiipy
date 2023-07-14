"""
Discretizing on a 2.5d staggered grid with no flux limiting - just FTCS and matsuno on a C grid

grid is:
   i h ip
j  P U P
h  V   V
jp P U P

k is the vertical component
"""
import math
import unittest
from collections import namedtuple, defaultdict

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

import constants
from constants import *
import low_pass
from coordinates_3d import *
import temperature
import geometry
from geometry import *
from grey_solar import grey_solar, grey_radiation, basic_grey_radiation, zenith_angle
from humidity import *
from dynamics import *


old_t = standard_temperature


def calc_energy(p, u, v, t, q, g, geom):
    """
    https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001549
    """
    u_at_center = imh(u)
    v_at_center = jmh(v)
    mag = np.sqrt(u_at_center ** 2 + v_at_center ** 2)

    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    rho = tp / (constants.Rd * tt)
    dp = p * geom.dsig
    geopotential_depth = (dp / (rho * G)).to_base_units()

    airmass = rho * geopotential_depth * geom.area

    total_depth = np.cumsum(geopotential_depth, 0)
    geopotential = total_depth * airmass * G
    geo = np.sum(geopotential).to(units.J)

    ke = mag ** 2 * .5 * airmass
    ke = np.sum(ke).to(units.J)

    ate = tt * Cp * airmass
    ate = np.sum(ate).to(units.J)
    return ke, ate, geo, ke+ate+geo


STATS = defaultdict(list)


def solar_timestep(t, p, g, dt, utc, geom):
    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    dt_air, dt_ground = basic_grey_radiation(p, tp, tt, g, 0.1, 0.9, 0.3, utc, geom)
    gt_n = g.gt + dt_ground * dt
    # gt_n = 275 * units.K
    tt_n = tt + dt_air * dt
    t_n = temperature.to_potential_temp(tt_n, tp)
    g_n = GroundVars(gt_n, g.gw, g.snow, g.ice)
    return t_n, g_n



def full_timestep(p, u, v, t, q, g, dt, utc, geom):
    # atmosphere timestep
    # print("ke:", calc_energy(p, u, v, t, q, g, geom))
    # exit()
    p, u, v, t, q = matsuno_timestep(p, u, v, t, q, dt, geom)
    print("utc:", utc.to(units.days))
    STATS["u_max"].append(np.max(u))
    STATS["u_min"].append(np.min(u))
    STATS["v_max"].append(np.max(v))
    STATS["v_min"].append(np.min(v))
    # print("p:", p.m)
    # print("u:", u.m)
    STATS["ke"].append(calc_energy(p, u, v, t, q, g, geom))
    # print("ke:", calc_energy(p, u, v, t, q, g, geom))
    # print("v:", v.m)
    return p, u, v, t, q, g

    # physics timestep
    t_n, g_n = solar_timestep(t, p, g, dt, utc, geom)
    # print("p:", p.m)
    # print("p:", np.max(p), np.min(p))
    # print("u:", np.max(u), np.min(u))
    # print("v:", np.max(v), np.min(v))
    # print("u:", u.m)
    # print("v:", v.m)
    return p, u, v, t_n, q, g_n
    tt_n = tt + dt_air * dt
    t_n = temperature.to_potential_temp(tt_n, tp)
    return p, u, v, t_n, q, g
    exit()
    dt_ground, dt_air, upwelling = grey_radiation(p, q, tt, 0.0, g, None, dt, geom)
    gt_n = g.gt + dt_ground * dt
    # gt_n = 275 * units.K
    tt_n = tt + dt_air * dt
    t_n = temperature.to_potential_temp(tt_n, tp)
    g_n = GroundVars(gt_n, g.gw, g.snow, g.ice)
    return p, u, v, t_n, q, g_n







height = 24
width = 36
layers = 9
# layers = 3




def plot_callback(q):
    quantity = q
    plt.clf()
    plt.imshow(quantity)
    # plt.title('n = %s' % (i,))
    # ax = plt.gca()
    # ax.format_coord = lambda x, y: f'{int(x + .5)} {int(y + .5)} {quantity[int(y + .5), int(x + .5)]}'
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated


PrognosticVars = namedtuple("PrognosticVars", ("p", "u", "v", "t", "q", "gt", "gw", "sd", "id"))
GroundVars = namedtuple("GroundVars", ("gt", "gw", "snow", "ice"))


def gen_initial_conditions(geom):
    full = (geom.layers, geom.height, geom.width)
    surface = (geom.height, geom.width)
    p = np.full(surface, 1) * 100000 * units.Pa - geom.ptop
    u = np.full(full, 1) * 1.0 * units.m / units.s
    v = np.full(full, 1) * .0 * units.m / units.s
    # tt = np.full(full, 1) * standard_temperature
    tt = np.full(full, 1) * 360 * units.K
    tp = p * geom.sig + geom.ptop
    t = temperature.to_potential_temp(tt, tp)
    q = np.full(full, 1) * 0.000003 * units.kg * units.kg ** -1
    q = unit_maximum(q, rh_to_mmr(manabe_rh(geom), tp, tt))

    # init ground
    # gt = np.full(surface, 1) * standard_temperature
    gt = np.full(surface, 1) * 360 * units.K
    gw = np.zeros(surface) * units.m
    snow = np.zeros(surface) * units.m
    ice = np.zeros(surface) * units.m
    
    g = GroundVars(gt, gw, snow, ice)

    return p, u, v, t, q, g


# class TestBasicDiscretizaion(unittest.TestCase):
#     def test_timestep_u_changes(self):
#         geom = gen_geometry(height, width, layers)
#         p, u, v, t, q, _ = gen_initial_conditions(geom)
#         dx = 100 * units.m
#         dt = 60 * 15 * units.s
#
#         # p[10, 10, 0] *= 1.01
#         # u[0, 3, 0] *= 200
#         # t[0, 3, 0] *= 1.0001
#         p[10, 10] *= 1.01
#         u[0, :, 12] *= 2
#         # u[:, 0, 3] *= 2
#         u[0, 0, 3] *= 2
#         # v[0, 18, 18] = 1 * units.m / units.s
#         # t[0, 3, 3] *= 1.1
#
#         # geom.heightmap[2, 3] = 2 * units.m
#         # u[3] *= 2
#         # ok, CFL for this is sqrt(2)/4
#
#         # t[2] += 1 * standard_temperature.units
#         # q[side_len//4:side_len//2] = 1
#         # q[2] = 1
#         # u[1] += .1 * u.units
#
#         # orig_u = u
#         # u = low_pass.avrx(u, geom)
#         # plt.imshow((orig_u - u).m)
#         # plt.ioff()
#         # plt.show()
#         # plt.plot(u[1])
#         # plt.plot(orig_u[1])
#         # plt.show()
#
#         plt.ion()
#         for i in tqdm(range(100000)):
#             p, u, v, t, q = matsuno_timestep(p, u, v, t, q, dt, geom)
#
#             # plot_callback(temperature.to_true_temp(t, p).m)
#             # plot_callback((t[0]).m[ :, :])
#             plot_callback(p.m[:, :])
#             if np.isnan(u).any() != False:
#                 break
#
#         plt.ioff()
#         plt.show()


def run_model(height, width, layers, dt, timesteps, callback):
    geom = gen_geometry(height, width, layers, sig_func=geometry.manabe_sig)
    p, u, v, t, q, g = gen_initial_conditions(geom)
    utc = 0 * units.hours
    # u *= 0
    v[0,0,0] = 0.1 * v.u
    u *= 0

    # p[0, 0] *= 1.01

    for i in tqdm(range(timesteps)):
        p, u, v, t, q, g = full_timestep(p, u, v, t, q, g, dt, utc, geom)
        utc += dt
        if callback:
            callback(p, u, v, t, q)

    return p, u, v, t, q, g, geom




def four_band_lw():
    # defaults from MITgcm/aim:
    # these are layer absorptivities per dp = 10^5 Pa
    # the water vapor terms are expressed for dq = 1 g/kg
    ABLWIN = 0.7
    ABLCO2 = 4.0
    ABLWV1 = 0.7
    ABLWV2 = 50.0


def test_absorbtion_units():
    dp = 10 * units.hPa



def main():
    # test_humidity_calcs()
    # exit()


    # run_model(height, width, layers, 60 * 15 * units.s, 1000, None)
    # p, u, v, t, q, g, geom = run_model(1, 1, 18, 60 * 15 * units.s, 3, None)
    p, u, v, t, q, g, geom = run_model(8, 8, 3, 60 * 30 * units.s, 1200 * 12, None)
    print("ground temp:", g.gt)
    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    print("atmosphere temps:", tt)
    print("pressures:", tp)


if __name__ == "__main__":
    main()

