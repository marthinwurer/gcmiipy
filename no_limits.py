"""
Discretizing on a 1d staggered grid with no flux limiting - just FTCS and matsuno on a C grid

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
import matplotlib.pyplot as plt

import constants
import low_pass
from coordinates_1d import *
from constants import *
import temperature

"""

dpuq / dt = -del.puq

(puq_n - puq) / dt = 

I need divergence at the edge of the tile
that's the flux across the tile edge
That's PU


"""




def advec_q(u, q, dx):
    # usual central differencing (A scheme) uses the centers of the neighboring grids
    # C scheme uses the grid edges as the points for fluxes
    # that means we need to calculate the fluxes at the edges for Q
    q_ph = iph(q)
    q_mh = imh(q)
    u_m = im(u)
    # we want (ph - mh) / dx * u

    q_flux = ((q_ph * u) - (q_mh * u_m)) / dx

    return q_flux


def calc_pu(u, p):
    # PU is pressure times velocity, effectively momentum, at plus half
    return u * iph(p)

def un_pu(pu, p):
    return pu / iph(p)


def advec_p(pu, dx):

    return div(pu, dx)


def advec_pu(p, pu, u, dx):

    puu = pu * u
    puum = imh(u) ** 2 * p
    puup = iph(u) ** 2 * iph(p)

    # puumpu = imh(puu)
    # puuppu = iph(puu)
    #
    # respu = (puuppu - puumpu) / dx


    res = (puup - puum) / dx
    return res


def advec_t(pu, t, dx):

    return div(pu * iph(t), dx)




def pgf(p, t, dx):
    # find the change in velocity due to the pressure gradient
    # this occurs at half
    pph = iph(p)
    tph = iph(t)

    # tt == True Temperature, from potential temperature
    tt = temperature.to_true_temp(tph, pph)
    rho = pph / (constants.Rd * tt)

    pgf = pph / rho * gradh(p, dx)

    return pgf


def half_timestep(p, u, t, q, sp, su, st, sq, dt, dx):
    # start with advection of a tracer (q)
    # p_h = iph(p)
    # sp_h = iph(sp)
    # pq = p * q

    # div_pq = div(pq, dx)

    pu = calc_pu(u, p)
    spu = calc_pu(su, sp)

    q_n = q - advec_q(su, sq, dx) * dt

    # f_p = advec_p(spu, dx) * dt
    # print("max advec_p", max(f_p))
    p_n = p - advec_p(spu, dx) * dt

    # f_pu = un_pu(advec_pu(sp, spu, su, dx) * dt, p_n)
    # f_pu_old = un_pu(advec_pu(sp, spu, su, dx) * dt, p)
    # f_pgf = un_pu(pgf(sp, st, dx), p_n) * dt
    # print("max pgf", max(f_pgf))
    pu_n = pu - (advec_pu(sp, spu, su, dx) + pgf(sp, st, dx)) * dt

    u_n = un_pu(pu_n, p_n)
    # f_unold = un_pu(pu_n, p)

    t_n = t - (advec_t(spu, st, dx) / p_n) * dt





    return (p_n, u_n, t_n, q_n)


def matsuno_timestep(p, u, t, q, dt, dx):
    sp, su, st, sq = half_timestep(p, u, t, q, p, u, t, q, dt, dx)
    return half_timestep(p, u, t, q, sp, su, st, sq, dt, dx)



side_len = 128


def plot_callback(q):
    quantity = q
    plt.clf()
    plt.plot(quantity)
    # plt.title('n = %s' % (i,))
    # ax = plt.gca()
    # ax.format_coord = lambda x, y: f'{int(x + .5)} {int(y + .5)} {quantity[int(y + .5), int(x + .5)]}'
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated


class TestBasicDiscretizaion(unittest.TestCase):
    def test_advec_q(self):
        p = np.full(side_len, 1) * standard_pressure
        u = np.full(side_len, 1) * 1.0 * units.m / units.s
        q = np.full(side_len, 1) * 0.1 * units.dimensionless
        t = np.full(side_len, 1) * temperature.to_potential_temp(standard_temperature, p)
        dx = 100 * units.m
        dt = 1 * units.s

        # t[2] += 1 * standard_temperature.units
        q[side_len//4:side_len//2] = 1
        # u[1] += .1 * u.units

        advec_q(u, q, dx)

    def test_timestep(self):
        p = np.full(side_len, 1) * standard_pressure
        u = np.full(side_len, 1) * 1.0 * units.m / units.s
        q = np.zeros(side_len) * units.dimensionless
        t = np.full(side_len, 1) * temperature.to_potential_temp(standard_temperature, p)
        dx = 100 * units.m
        dt = 1 * units.s

        # t[2] += 1 * standard_temperature.units
        q[side_len//4:side_len//2] = 1
        # q[2] = 1
        # u[1] += .1 * u.units

        plt.ion()
        for i in tqdm(range(100000)):
            p, u, t, q = matsuno_timestep(p, u, t, q, dt, dx)

            plot_callback(t.m)
            if np.isnan(u).any() != False:
                break

        plt.ioff()
        plt.show()

    def test_timestep_u_changes(self):
        p = np.full(side_len, 1) * standard_pressure
        u = np.full(side_len, 1) * 1.0 * units.m / units.s
        q = np.full(side_len, 1) * 0.1 * units.dimensionless
        t = np.full(side_len, 1) * temperature.to_potential_temp(standard_temperature, p)
        dx = 70000 * units.m
        dt = 60. * 15 * units.s

        # ok, 0.17s at 100m and 2m/s is the edge of stability.
        # 0.1s and 334m/s is also the edge of stability.
        # What's the actual stability formula here?
        # well, with the standard cfl we get 0.0034 for 0.17, 100m, and 2m/s
        # now, with 1.7s and 100m and 2m/s, we're barely stable
        # with it bumped up to 4 m/s, it's still stable. interesting.
        # bumped up to 250m/s, still stable. this is wild.
        # 275 is barely stable, 280 barely counts, 290,
        # 295 is the last one that stabilizes
        # 10km, 10.7s, 296m/s
        # 17s fails, 16s barely stabilizes
        # 100km, 15min, 296 works.

        # p = p / 2.
        p[3] *= 1.00001
        # u[3] *= 1.001
        # ok, CFL for this is sqrt(2)/4
        # maybe CFL for pure advection, for gravity waves it's probably
        # sqrt(gH)*dt/dx
        # where gH is g * geopotential / g
        # geopotential is p / rho
        # the limiting factor will be the gravity waves CFL term.
        # in full it's actually c = u +- sqrt(gH)
        # from https://www2.atmos.umd.edu/~ekalnay/syllabi/AOSC614/NWP-CH03-2-4.pdf
        # We may assume u < 100 m/s and sqrt(gH) < 300 m/s, so a safe
        # maximum value for c is 400 m/s
        # GCMII solves this by low pass filtering the gravity waves when
        # the dx is too small for the timestep
        # ok, so the if we want to find dt, then we do CFL * dx / c

        # t[2] += 1 * standard_temperature.units
        q[side_len//4:side_len//2] = 1
        # q[2] = 1
        # u[1] += .1 * u.units

        plt.ion()
        for i in tqdm(range(100000)):
            # test filtering with an averaging filter
            print(max(u))
            # we want to filter out waves faster than the vertial difference
            # frequency, which at 15m and 1/24 of circumfrence is 833km over 15min,
            # which has
            # u = low_pass.butter_lowpass_filter(u,  )
            # u = (u + ip(u) + im(u)) / 3
            # u = (u + ip(u) + im(u)) / 3
            # u = (u + ip(u) + im(u)) / 3
            # u = (u + ip(u) + im(u)) / 3
            print(max(u))

            p, u, t, q = matsuno_timestep(p, u, t, q, dt, dx)
            print(max(p))

            plot_callback(p.m)
            if np.isnan(u).any() != False:
                break

        plt.ioff()
        plt.show()