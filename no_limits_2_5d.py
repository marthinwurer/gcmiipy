"""
Discretizing on a 1d staggered grid with no flux limiting - just FTCS and matsuno on a C grid

grid is:
   i h ip
j  P U P
h  V   V
jp P U P
"""
import math
import unittest

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

import constants
from coordinates import *
from constants import *
import temperature


def calc_pu(p, u):
    pu = u * iph(p)
    return pu


def calc_pv(p, v):
    pv = v * jph(p)
    return pv


def un_pu(pu, p):
    u = pu / iph(p)
    return u


def un_pv(pv, p):
    v = pv / jph(p)
    return v


def advec_p(pu, pv, geom):

    dp = (pu - imj(pu)) / geom.dx_j + (pv - ijm(pv)) / geom.dy
    return dp


def advec_m(p, u, v, geom):
    """
       i h ip
    j  P U P
    h  V   V
    jp P U P
    """
    vph = iph(v)
    p_mid = iph(jph(p))

    puum = imh(u) ** 2 * p
    puup = ipj(puum)
    # puvm is at j-h, i+h
    puvm = jmh(u) * ijm(vph) * ijm(p_mid)
    puvp = ipj(puvm)

    dut = (puum - puup) / geom.dx_j + (puvm - puvp) / geom.dy

    pvvm = jmh(v) ** 2 * p
    pvvp = ijp(pvvm)
    # pvum is at i-h, j+h
    pvum = imj(p_mid) * imh(v) * imj(jph(u))
    pvup = ipj(pvum)

    dvt = (pvvm - pvvp) / geom.dy + (pvum - pvup) / geom.dx_h

    return (dut, dvt)


def pgf(p, t, geom):
    ppih = iph(p)
    ttu = temperature.to_true_temp(iph(t), ppih)
    rhou = ppih / (constants.Rd * ttu)

    pgfu = ppih / rhou * gradi(p, geom.dx_j)

    ppjh = jph(p)
    ttv = temperature.to_true_temp(jph(t), ppjh)
    rhov = ppjh / (constants.Rd * ttv)

    pgfv = ppjh / rhov * gradj(p, geom.dy)

    return pgfu, pgfv


def advec_t(pu, pv, t, geom):

    tpu = pu * iph(t)
    tpv = pv * jph(t)

    dt = (tpu - imj(tpu)) / geom.dx_j + (tpv - ijm(tpv)) / geom.dy

    return dt




def half_timestep(p, u, v, t, q, sp, su, sv, st, sq, dt, geom):
    pu = calc_pu(p, u)
    spu = calc_pu(sp, su)
    pv = calc_pv(p, v)
    spv = calc_pv(sp, sv)

    p_n = p - advec_p(spu, spv, geom) * dt

    dut, dvt = advec_m(sp, su, sv, geom)
    pgu, pgv = pgf(sp, st, geom)

    pu_n = pu - (dut + pgu) * dt
    pv_n = pv - (dvt + pgv) * dt

    u_n = un_pu(pu_n, p_n)
    v_n = un_pv(pv_n, p_n)

    t_n = t - (advec_t(spu, spv, st, geom) / p_n) * dt

    v_n[-1, :] *= 0
    # u_n[:, -1] *= 0

    return (p_n, u_n, v_n, t_n, q)


def matsuno_timestep(p, u, v, t, q, dt, geom):
    sp, su, sv, st, sq = half_timestep(p, u, v, t, q, p, u, v, t, q, dt, geom)
    return half_timestep(p, u, v, t, q, sp, su, sv, st, sq, dt, geom)






height = 24
width = 36
layers = 9


class Geom:
    def __init__(self):
        self.sigma_edges = []
        self.dsig = []
        # self.dx = 0
        self.dy = 0
        self.lat = []
        self.dx_j = 0
        self.dx_h = 0


def gen_geometry(height, width, layers):
    """"""
    """
      DATA SIGE/    1.,.948665,.866530,.728953,.554415,.390144,
     *  .251540,.143737,.061602,28*0./                        
    """
    sige = [1., .948665, .866530, .728953, .554415, .390144, .251540, .143737, .061602, 0.]

    """
C**** CALCULATE DSIG AND DSIGO                                           816.   
      DO 700 L=1,LM                                                      817.   
  700 DSIG(L)=SIGE(L)-SIGE(L+1)                                          818.   
      DO 710 L=1,LM-1                                                    819.   
  710 DSIGO(L)=SIG(L)-SIG(L+1)                                           820.   
    """

    # DATA PLB4 in R83ZAmacDBL.f seems to have the base layer values for pressure
    # I assume that sigma is calculated off of that?

    plb4_4 = [
             1013.2500, 1000.0000, 950.0000, 900.0000, 850.0000, 800.0000,
              750.0000, 700.0000, 650.0000, 600.0000, 550.0000, 500.0000,
              450.0000, 400.0000, 350.0000, 300.0000, 250.0000, 200.0000,
              150.0000, 100.0000, 50.0000, 20.0000, 10.0000, 5.0000,
                2.0000, 1.0000, 0.5000, 0.2000, 0.1000, 0.0500,
                0.0200, 0.0100, 0.0050, 0.0020, 0.0010, 1.E-05,
                0.
    ]

    mysig = []
    for i in range(layers+1):
        mysig.append(1 - i/(layers))

    circumference = 2 * radius * math.pi
    lat_j = np.zeros((height,))
    lat_h = np.zeros((height,))
    dlat = 180 / height
    dlong = 360 / width

    sin_j = np.zeros((height,))
    cos_j = np.zeros((height,))
    sin_h = np.zeros((height,))
    cos_h = np.zeros((height,))
    for i in range(height):
        lat_j[i] = 90 - (i+0.5) * dlat
        lat_h[i] = 90 - (i+1) * dlat

    cos_j = np.cos(lat_j * np.pi / 180)
    sin_j = np.sin(lat_j * np.pi / 180)
    cos_h = np.cos(lat_h * np.pi / 180)
    sin_h = np.sin(lat_h * np.pi / 180)
    dx_j = cos_j * circumference / width
    dx_h = cos_h * circumference / width

    print(lat_j)
    print(lat_h)
    print(cos_j)
    print(dx_j)
    print(dx_h)






    # plt.plot(plb4_4)
    # plt.plot(sige)
    # plt.plot(mysig)
    # plt.show()

    geom = Geom()
    geom.sigma_edges = mysig
    # TODO spherical geometry
    # geom.dx = circumference / width
    geom.dx_j = np.reshape(dx_j, (height, 1))
    geom.dx_h = np.reshape(dx_h, (height, 1))
    geom.dy = circumference / 2 / height

    return geom



def plot_callback(q):
    quantity = q
    plt.clf()
    plt.imshow(quantity)
    # plt.title('n = %s' % (i,))
    # ax = plt.gca()
    # ax.format_coord = lambda x, y: f'{int(x + .5)} {int(y + .5)} {quantity[int(y + .5), int(x + .5)]}'
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated


class TestBasicDiscretizaion(unittest.TestCase):
    def test_timestep_u_changes(self):
        geom = gen_geometry(height, width, layers)
        # p = np.full((height, width, layers), 1) * standard_pressure
        # u = np.full((height, width, layers), 1) * 1.0 * units.m / units.s
        # v = np.full((height, width, layers), 1) * .0 * units.m / units.s
        # q = np.full((height, width, layers), 1) * 0.1 * units.dimensionless
        # t = np.full((height, width, layers), 1) * temperature.to_potential_temp(standard_temperature, p)
        p = np.full((height, width), 1) * standard_pressure
        u = np.full((height, width), 1) * 1.0 * units.m / units.s
        v = np.full((height, width), 1) * .0 * units.m / units.s
        q = np.full((height, width), 1) * 0.1 * units.dimensionless
        t = np.full((height, width), 1) * temperature.to_potential_temp(standard_temperature, p)
        dx = 100 * units.m
        dt = 60 * units.s

        # p[10, 10, 0] *= 1.01
        # u[0, 3, 0] *= 200
        # t[3, 3, 0] *= 1.1
        p[10, 10] *= 1.01
        u[1, 3] *= 2
        t[3, 3] *= 1.1
        # u[3] *= 2
        # ok, CFL for this is sqrt(2)/4

        # t[2] += 1 * standard_temperature.units
        # q[side_len//4:side_len//2] = 1
        # q[2] = 1
        # u[1] += .1 * u.units

        plt.ion()
        for i in tqdm(range(100000)):
            p, u, v, t, q = matsuno_timestep(p, u, v, t, q, dt, geom)

            # plot_callback(temperature.to_true_temp(t, p).m)
            plot_callback(p.m[:, :])
            # plot_callback(p.m[:, :, 0])
            if np.isnan(u).any() != False:
                break

        plt.ioff()
        plt.show()
