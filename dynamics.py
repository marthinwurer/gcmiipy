import math

import numpy as np

import constants
from constants import *
import low_pass
from coordinates_3d import *
import temperature
import geometry
from geometry import *



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
    conv = ((pu - imj(pu)) / geom.dx_j + (pv - ijm(pv)) / geom.dy) * geom.dsig
    pit = np.sum(conv, 0)

    sd = np.cumsum(conv[::-1], 0)[::-1] - pit * geom.sigb
    # apply boundary condition to sd
    sd[0] = 0 * sd.u

    return (pit, sd)


def advec_sig(sd, q, geom):
    flux = kmh(q) * sd
    dq = (flux - kp(flux)) / geom.dsig
    return -dq


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

    # coriolis effect
    # need pu at pv and pv at pu
    # if geom.width != 1:
    if False:
        pu_at_pv = imh(jph(pu))
        pv_at_pu = iph(jmh(pv))

        w = 2 * math.pi / units.day

        cp_at_u = 2 * np.sin(geom.lat) * w
        cp_at_v = 2 * np.sin(jph(geom.lat)) * w

        coriolis_u = cp_at_u * -pv_at_pu
        coriolis_v = cp_at_v * pu_at_pv
    else:
        coriolis_u = 0 * units.kg * units.s ** -4
        coriolis_v = coriolis_u

    # print(w.to_base_units().u)
    # print(cp_at_u.to_base_units().u)
    # print(coriolis_u.to_base_units().u)



    dut = (puum - puup) / geom.dx_j + (puvm - puvp) / geom.dy + coriolis_u
    dvt = (pvvm - pvvp) / geom.dy + (pvum - pvup) / geom.dx_h + coriolis_v

    # print(dut.m)

    return (dut, dvt)


def compute_geopotential(p, t, geom):
    # True pressure, true temperature, and density
    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    rho = tp / (constants.Rd * tt)
    dp = p * geom.dsig
    geopotential_depth = (dp / (rho * G)).to_base_units()
    phi_mine = geom.heightmap + np.cumsum(geopotential_depth, 0) - (geopotential_depth / 2)
    phi_mine *= constants.G

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

    print(phi_theirs[0,0,0])
    print(phi_mine[0,0,0])
    
    assert phi_mine.to_base_units().u == phi_theirs.to_base_units().u

    return phi_theirs
    # return phi_mine



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
    assert pu.shape == pgfu.shape
    # phiu = low_pass.arakawa_1977(phiu, geom)

    pu_n = pu - (dut + dus + pgfu) * dt
    pv_n = pv - (dvt + dvs + phiv + pgv) * dt
    # pu_n = pu - (dut + pgu + dus) * dt
    # pv_n = pv - (dvt + pgv + dvs) * dt

    u_n = un_pu(pu_n, p_n)
    v_n = un_pv(pv_n, p_n)

    t_n = (t * p - (advec_t(spu, spv, st, geom) + advec_sig(sd, st, geom)) * dt) / p_n
    # t_n = t - (advec_t(spu, spv, st, geom) + advec_sig(sd, st, geom)) * dt / p_n

    # TODO advect water vapor as well
    # TODO might need to flux limit this
    q_n = (q * p - (advec_t(spu, spv, sq, geom) + advec_sig(sd, sq, geom)) * dt) / p_n

    # TODO do boundary conditions outside of the timestep
    v_n[:, -1, :] *= 0
    # u_n[:, -1] *= 0
    # print()
    # print(p[0, 3], p_n[0, 3], pgfu[0, 0, 3], dus[0, 0, 3])

    return (p_n, u_n, v_n, t_n, q_n)


def matsuno_timestep(p, u, v, t, q, dt, geom, boundary_conditions=None):
    sp, su, sv, st, sq = half_timestep(p, u, v, t, q, p, u, v, t, q, dt, geom)
    if boundary_conditions:
        sp, su, sv, st, sq = boundary_conditions(sp, su, sv, st, sq, dt, geom)
    op, ou, ov, ot, oq = half_timestep(p, u, v, t, q, sp, su, sv, st, sq, dt, geom)
    if boundary_conditions:
        op, ou, ov, ot, oq = boundary_conditions(op, ou, ov, ot, oq, dt, geom)
    return op, ou, ov, ot, oq

