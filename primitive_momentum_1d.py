"""
Discretizing the full primitive equations in momentum form on a 1d staggered grid

grid is:
   i h ip
j  P U P

U is effectively at iph


Goal: get GCMII formulas working on a 1d grid.
Key variables:
U: velocity
P: Pressure
T: Potential Temperature
Q: Moisture content

Derived for calculations:
PU: momentum
PT: mixed temperature
PQ: mixed moisture

"""
import unittest
from tqdm import tqdm

from coordinates_1d import *
from constants import *


def limit_flux(q, u, dx):
    # determine the stable value for q at h, the interface.
    dir = u.m < 0

    dif = ip(q) - q
    q_h = np.copy(q.m) * q.units
    q_h[dir] = ip(q)[dir]
    return q_h * u


def advect_q(q_i, pu_h, dx):
    return div(limit_flux(q_i, pu_h, dx), dx)


def advect_u(u_h, pu_h, dx):
    #
    pu_ip = iph(pu_h)
    return div(limit_flux(u_h, pu_ip, dx), dx)

    return divu(pu_h * u_h, dx)


def half_timestep(p, u, t, q, sp, su, st, sq, dt, dx):
    # derive everything
    p_h = iph(p)
    sp_h = iph(sp)
    pu_h = p_h * u
    spu_h = sp_h * su
    pt_i = p * t
    pq_i = p * q

    p_n = p - dt * div(spu_h, dx)  # + sigma term

    rho_h = iph(sp / (Rd * potential_temp_to_temp(sp, st)))

    pu_n = pu_h - dt * (
            advect_u(su, spu_h, dx)
            # divu(spu_h * su, dx)  # advection term
            +
            (sp_h / rho_h) * gradh(sp, dx)  # pgf term
    )  # - divu(
    pt_n = pt_i - dt * advect_q(st, spu_h, dx)
    pq_n = pq_i - dt * advect_q(sq, spu_h, dx)
    u_n = pu_n / p_n
    t_n = pt_n / p_n
    q_n = pq_n / p_n

    return p_n, u_n, t_n, q_n


def matsuno_timestep(p, u, t, q, dt, dx):
    sp, su, st, sq = half_timestep(p, u, t, q, p, u, t, q, dt, dx)
    return half_timestep(p, u, t, q, sp, su, st, sq, dt, dx)


import matplotlib.pyplot as plt

side_len = 8


def plot_callback(q):
    quantity = q
    plt.clf()
    plt.plot(quantity)
    # plt.title('n = %s' % (i,))
    # ax = plt.gca()
    # ax.format_coord = lambda x, y: f'{int(x + .5)} {int(y + .5)} {quantity[int(y + .5), int(x + .5)]}'
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated


class TestPrimMomentum(unittest.TestCase):
    def test_timestep(self):
        p = np.full(side_len, 1) * standard_pressure
        u = np.full(side_len, 1) * 1.0 * units.m / units.s
        q = np.zeros(side_len) * units.dimensionless
        t = np.full(side_len, 1) * standard_temperature
        dx = 100 * units.m
        dt = 1 * units.s

        t[2] += 1 * standard_temperature.units
        q[2] = 1
        u[1] += .1 * u.units

        plt.ion()
        for i in tqdm(range(1000)):
            p, u, t, q = matsuno_timestep(p, u, t, q, dt, dx)

            plot_callback(u.m)
            if np.isnan(u).any() != False:
                break

        plt.ioff()
        plt.show()


