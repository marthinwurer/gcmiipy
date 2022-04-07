"""
Discretizing on a 1d staggered grid with no flux limiting - just FTCS and matsuno on a C grid

grid is:
   i h ip
j  P U P
h  V   V
jp P U P
"""

import unittest
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


def advec_p(pu, pv, dx):

    dp = (pu - imj(pu)) / dx + (pv - ijm(pv)) / dx
    return dp


def advec_m(p, u, v, dx):
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

    dut = (puum - puup) / dx + (puvm - puvp) / dx

    pvvm = jmh(v) ** 2 * p
    pvvp = ijp(pvvm)
    # pvum is at i-h, j+h
    pvum = imj(p_mid) * imh(v) * imj(jph(u))
    pvup = ipj(pvum)

    dvt = (pvvm - pvvp) / dx + (pvum - pvup) / dx

    return (dut, dvt)


def pgf(p, t, dx):
    ppih = iph(p)
    ttu = temperature.to_true_temp(iph(t), ppih)
    rhou = ppih / (constants.Rd * ttu)

    pgfu = ppih / rhou * gradi(p, dx)

    ppjh = jph(p)
    ttv = temperature.to_true_temp(jph(t), ppjh)
    rhov = ppjh / (constants.Rd * ttv)

    pgfv = ppjh / rhov * gradj(p, dx)

    return pgfu, pgfv


def advec_t(pu, pv, t, dx):

    tpu = pu * iph(t)
    tpv = pv * jph(t)

    dt = (tpu - imj(tpu)) / dx + (tpv - ijm(tpv)) / dx

    return dt




def half_timestep(p, u, v, t, q, sp, su, sv, st, sq, dt, dx):
    pu = calc_pu(p, u)
    spu = calc_pu(sp, su)
    pv = calc_pv(p, v)
    spv = calc_pv(sp, sv)

    p_n = p - advec_p(spu, spv, dx) * dt

    dut, dvt = advec_m(sp, su, sv, dx)
    pgu, pgv = pgf(sp, st, dx)

    pu_n = pu - (dut + pgu) * dt
    pv_n = pv - (dvt + pgv) * dt

    u_n = un_pu(pu_n, p_n)
    v_n = un_pv(pv_n, p_n)

    t_n = t - (advec_t(spu, spv, st, dx) / p_n) * dt

    v_n[-1, :] *= 0
    u_n[:, -1] *= 0

    return (p_n, u_n, v_n, t_n, q)


def matsuno_timestep(p, u, v, t, q, dt, dx):
    sp, su, sv, st, sq = half_timestep(p, u, v, t, q, p, u, v, t, q, dt, dx)
    return half_timestep(p, u, v, t, q, sp, su, sv, st, sq, dt, dx)






height = 24
width = 36


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
        p = np.full((height, width), 1) * standard_pressure
        u = np.full((height, width), 1) * 1.0 * units.m / units.s
        v = np.full((height, width), 1) * .0 * units.m / units.s
        q = np.full((height, width), 1) * 0.1 * units.dimensionless
        t = np.full((height, width), 1) * temperature.to_potential_temp(standard_temperature, p)
        dx = 100 * units.m
        dt = .1 * units.s

        p[10, 10] *= 1.01
        u[0, 3] *= 200
        t[3, 3] *= 1.1
        # u[3] *= 2
        # ok, CFL for this is sqrt(2)/4

        # t[2] += 1 * standard_temperature.units
        # q[side_len//4:side_len//2] = 1
        # q[2] = 1
        # u[1] += .1 * u.units

        plt.ion()
        for i in tqdm(range(100000)):
            p, u, v, t, q = matsuno_timestep(p, u, v, t, q, dt, dx)

            # plot_callback(temperature.to_true_temp(t, p).m)
            plot_callback(p.m)
            if np.isnan(u).any() != False:
                break

        plt.ioff()
        plt.show()
