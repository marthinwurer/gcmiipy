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

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

import constants
import low_pass
from coordinates_3d import *
from constants import *
import temperature
from geometry import *
from grey_solar import grey_solar


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


def aflux(pu, pv, geom):
    """
    compute pit (dp), conv, and sigma_dot
    """
    """
    So, the convergence seems to be u * surface pressure * dsig, which is interesting.
    I was thinking they'd be each layer * sigma and then integrated over that layer?
    So, what would that actually be, then? 
    the actual pressure is sig * p + ptop
    then that gets fed into rho for the pgf term
    well, isn't that that? dsig is the amount of pressure on that level, 
    so yeah, that would work.
    so then sigma dot is the diffence between what
    sigma would be for that level vs what it should be
    ok yeah, dsig * p is the amount of pressure for that level, and then that gets advected
    according to the velocity for that level
    That means that my layers are actually a bit thinner,
    so I might be able to take longer time steps.
    so, we have p which is at sig = 1
    then conv which is the change in pressure for each layer
    then pit which is the sum.
    """
    conv = ((pu - imj(pu)) / geom.dx_j + (pv - ijm(pv)) / geom.dy) * geom.dsig
    pit = np.sum(conv, 0)

    # sd1 = conv - geom.dsig * pit
    """
    Ok, so we get the change in pressure at each layer, right? That's conv
    we also get the change in surface pressure, that's pit
    Ah, so we do have to sum up the changes in the lower levels to get our 
    actual change at sigma.
    I think I need to work out a few layers out by hand
    
    
    """
    sd = np.cumsum(conv[::-1], 0)[::-1] - pit * geom.sigb
    # apply boundary condition to sd
    sd[0] = 0 * sd.u

    return (pit, sd)


def advec_sig(sd, q, geom):
    flux = kmh(q) * sd
    dq = (flux - kp(flux)) / geom.dsig
    return -dq




# def advec_p(pu, pv, geom):
    #
    # dp = (pu - imj(pu)) / geom.dx_j + (pv - ijm(pv)) / geom.dy
    # return dp


# def advec_m(p, u, v, geom):
#     """
#        i h ip
#     j  P U P
#     h  V   V
#     jp P U P
#     """
#     vph = iph(v)
#     p_mid = iph(jph(p))
#
#     puum = imh(u) ** 2 * p
#     puum = low_pass.avrx(puum, geom)
#     puup = ipj(puum)
#     # puvm is at j-h, i+h
#     puvm = jmh(u) * ijm(vph) * ijm(p_mid)
#     puvp = ipj(puvm)
#
#     dut = (puum - puup) / geom.dx_j + (puvm - puvp) / geom.dy
#
#     pvvm = jmh(v) ** 2 * p
#     pvvp = ijp(pvvm)
#     # pvum is at i-h, j+h
#     pvum = imj(p_mid) * imh(v) * imj(jph(u))
#     # TODO I might have to use pu for this function instead to deal with this stuff
#     pvum = low_pass.avrx(pvum, geom)
#     pvup = ipj(pvum)
#
#     dvt = (pvvm - pvvp) / geom.dy + (pvum - pvup) / geom.dx_h
#
#     return (dut, dvt)


def advec_m_pu(p, u, v, pu, pv, geom):
    # puum = imh(u) ** 2 * p
    puum = imh(u) * imh(pu)
    puup = ipj(puum)

    # # puvm is at j-h, i+h
    # puvp should be at jph, iph
    # pv is at i, j+h
    # puvm = jmh(u) * ijm(vph) * ijm(p_mid)
    # puvp = ipj(puvm)
    puvp = iph(pv) * jph(u)
    puvm = ijm(puvp)


    # v, pv: i, jph
    # pu: iph, j
    # pvvm i, j
    # pvvp i, jp1
    # pvum imh, jph
    # pvup iph, jph
    pvvm = jmh(v) * jmh(pv)
    pvvp = ijp(pvvm)
    pvup = iph(v) * jph(pu)
    pvum = imj(pvup)

    dut = (puum - puup) / geom.dx_j + (puvm - puvp) / geom.dy
    dvt = (pvvm - pvvp) / geom.dy + (pvum - pvup) / geom.dx_h

    return (dut, dvt)


def compute_geopotential(p, t, geom):
    # True pressure, true temperature, and density
    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    rho = tp / (constants.Rd * tt)

    sp = geom.sig * p

    spa = sp / rho
    s1 = spa * geom.dsig
    # this is the first phi in the original pgf
    pkdn = ((geom.sig * p + geom.ptop) / constants.P0) ** constants.kappa
    pkup = kp(pkdn)
    # for l+1
    stp = constants.Cp * kph(t) * (pkdn - pkup)
    s2 = geom.sigt * stp
    stp_n = km(stp)
    stp_n[0] = np.sum(s1 - s2, 0) + geom.heightmap * constants.G

    # this is the geopotential of each layer
    phi_theirs = np.cumsum(stp_n, 0)

    return phi_theirs



def pgf(p, t, geom):

    # True pressure, true temperature, and density
    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    rho = tp / (constants.Rd * tt)

    sp = geom.sig * p

    phi_theirs = compute_geopotential(p, t, geom)

    # phi at U and V
    phiu = iph(p) * gradi(phi_theirs, geom.dx_j)
    phiv = jph(p) * gradj(phi_theirs, geom.dy)

    ppih = iph(sp)
    rhou = iph(rho)
    # this part is basically SPA
    pgfu = ppih / rhou * gradi(p, geom.dx_j)

    ppjh = jph(sp)
    rhov = jph(rho)
    pgfv = ppjh / rhov * gradj(p, geom.dy)

    return pgfu, pgfv, phiu, phiv


def advec_t(pu, pv, t, geom):

    tpu = pu * iph(t)
    tpv = pv * jph(t)

    dt = (tpu - imj(tpu)) / geom.dx_j + (tpv - ijm(tpv)) / geom.dy

    return dt


def AzG(angle, top, bottom, urg, asg):
    """
    based off of formula 25 in manabe 64
    """

    return asg * (np.arccos(angle) * urg(top))


def solar(p, q, t, c, utc, geom):
    # do basic solar radiation, starting from the top of the atmosphere.
    # need more solar stuff

    # split between the cloudy and clear parts.

    # zenith_angle = 

    starting_flux = solar_constant * cos(zenith_angle)

    # pressure_depth = 

    # ok, so I need to be able to get optical depth from something.
    # I think we need to somehow get from whatever we have to cm/km
    # do we just have to mulitply that geopotential difference by the density?

    # do optical path from formula 44 of 
    # https://webspace.science.uu.nl/~delde102/AtmosphericDynamics%5b2020a%5dCh2.pdf

    # eventually take 4 band model from section 2.17 of that paper.



    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    rho = tp / (constants.Rd * tt)

    dp = p * geom.dsig

    # page 51
    sig_a = 0.125 * units.m ** 2 / units.kg

    optical_path = sig_a * rho_a / (G * rho) * dp






def half_timestep(p, u, v, t, q, sp, su, sv, st, sq, dt, geom):
    # u = low_pass.avrx(u, geom)
    # su = low_pass.avrx(su, geom)
    pu = calc_pu(p, u)
    spu_orig = calc_pu(sp, su)
    # pu = low_pass.avrx(pu_orig, geom)
    spu = low_pass.arakawa_1977(spu_orig, geom)
    pv = calc_pv(p, v)
    spv = calc_pv(sp, sv)

    pit, sd = aflux(spu, spv, geom)
    p_n = p - pit * dt

    # dut, dvt = advec_m(sp, su, sv, geom)
    dut, dvt = advec_m_pu(sp, su, sv, spu, spv, geom)
    pgu, pgv, phiu, phiv = pgf(sp, st, geom)
    dus = advec_sig(iph(sd), su, geom)
    dvs = advec_sig(jph(sd), sv, geom)

    pgfu = low_pass.arakawa_1977(pgu + phiu, geom)
    assert pgu.shape == pgfu.shape
    # phiu = low_pass.arakawa_1977(phiu, geom)

    pu_n = pu - (dut + dus + pgfu) * dt
    pv_n = pv - (dvt + dvs + phiv + pgv) * dt
    # pu_n = pu - (dut + pgu + dus) * dt
    # pv_n = pv - (dvt + pgv + dvs) * dt

    u_n = un_pu(pu_n, p_n)
    v_n = un_pv(pv_n, p_n)

    t_n = (t * p - (advec_t(spu, spv, st, geom) + advec_sig(sd, st, geom)) * dt) / p_n
    # t_n = t - (advec_t(spu, spv, st, geom) + advec_sig(sd, st, geom)) * dt / p_n

    v_n[:, -1, :] *= 0
    # u_n[:, -1] *= 0
    # print()
    # print(p[0, 3], p_n[0, 3], pgfu[0, 0, 3], dus[0, 0, 3])

    return (p_n, u_n, v_n, t_n, q)


def matsuno_timestep(p, u, v, t, q, dt, geom):
    sp, su, sv, st, sq = half_timestep(p, u, v, t, q, p, u, v, t, q, dt, geom)
    return half_timestep(p, u, v, t, q, sp, su, sv, st, sq, dt, geom)


def full_timestep(p, u, v, t, q, dt, geom):
    p, u, v, t, q = matsuno_timestep(p, u, v, t, q, dt, geom)
    grey_solar(p, q, t, 0.3, None, None, geom)
    return p, u, v, t, q 







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


class TestBasicDiscretizaion(unittest.TestCase):
    def test_timestep_u_changes(self):
        geom = gen_geometry(height, width, layers)
        p = np.full((height, width), 1) * standard_pressure - geom.ptop
        u = np.full((layers, height, width), 1) * 1.0 * units.m / units.s
        v = np.full((layers, height, width), 1) * .0 * units.m / units.s
        q = np.full((layers, height, width), 1) * 0.1 * units.dimensionless
        t = np.full((layers, height, width), 1) * temperature.to_potential_temp(standard_temperature, p)
        dx = 100 * units.m
        dt = 60 * 15 * units.s

        # p[10, 10, 0] *= 1.01
        # u[0, 3, 0] *= 200
        # t[0, 3, 0] *= 1.0001
        p[10, 10] *= 1.01
        u[0, :, 12] *= 2
        # u[:, 0, 3] *= 2
        u[0, 0, 3] *= 2
        # v[0, 18, 18] = 1 * units.m / units.s
        # t[0, 3, 3] *= 1.1

        # geom.heightmap[2, 3] = 2 * units.m
        # u[3] *= 2
        # ok, CFL for this is sqrt(2)/4

        # t[2] += 1 * standard_temperature.units
        # q[side_len//4:side_len//2] = 1
        # q[2] = 1
        # u[1] += .1 * u.units

        # orig_u = u
        # u = low_pass.avrx(u, geom)
        # plt.imshow((orig_u - u).m)
        # plt.ioff()
        # plt.show()
        # plt.plot(u[1])
        # plt.plot(orig_u[1])
        # plt.show()

        plt.ion()
        for i in tqdm(range(100000)):
            p, u, v, t, q = matsuno_timestep(p, u, v, t, q, dt, geom)

            # plot_callback(temperature.to_true_temp(t, p).m)
            # plot_callback((t[0]).m[ :, :])
            plot_callback(p.m[:, :])
            if np.isnan(u).any() != False:
                break

        plt.ioff()
        plt.show()


def run_model(height, width, layers, dt, timesteps, callback):
    geom = gen_geometry(height, width, layers)
    p = np.full((height, width), 1) * standard_pressure - geom.ptop
    u = np.full((layers, height, width), 1) * 1.0 * units.m / units.s
    v = np.full((layers, height, width), 1) * .0 * units.m / units.s
    q = np.full((layers, height, width), 1) * 0.1 * units.g * units.m ** -3
    tt = np.full((layers, height, width), 1) * standard_temperature
    t = temperature.to_potential_temp(tt, p * geom.sig + geom.ptop)
    # t = np.full((layers, height, width), 1) * temperature.to_potential_temp(standard_temperature, p)

    p[0, 0] *= 1.01

    for i in tqdm(range(timesteps)):
        p, u, v, t, q = full_timestep(p, u, v, t, q, dt, geom)
        if callback:
            callback(p, u, v, t, q)




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
    # run_model(height, width, layers, 60 * 15 * units.s, 1000, None)
    run_model(1, 1, 18, 60 * 15 * units.s, 1, None)


if __name__ == "__main__":
    main()

