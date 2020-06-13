
"""
The goal for this file is to implement Matsuno 1966 with a C grid.
The paper this is based off of:
https://www.jstage.jst.go.jp/article/jmsj1965/44/1/44_1_76/_pdf
"""
import numpy as np
import matplotlib.pyplot as plt
import tqdm as tqdm

from constants import units, unit_roll, G, get_total_variation
from coordinates import ijp, ijm, ipj, imj, imjp


def advection_of_velocity_u(u, v, dx):
    """
    grid is:
    i h ip
    P U P
    V   V
    P U P

    so we need to get u at half and v at half
    """
    """
    * = current U
    vu = V * U, both need to be interpolated
    uu = U * U, needs to be interpolated and squared
    im   i    ip
    U    U     U
    P  V vu V  P  j
    U uu *  uu U
    P  V vu V  P  jp
    U    U     U
    """
    # velocity terms
    u_ipj = (ipj(u) + u) / 2
    u_imj = (imj(u) + u) / 2
    # these names are relative to U, not to p
    v_ijm = (imj(v) + v) / 2
    v_ijp = (imjp(v) + ijp(v)) / 2

    # gradient terms
    du_ipj = (ipj(u) - u)
    du_imj = (u - imj(u))
    du_ijp = (ijp(u) - u)
    du_ijm = (u - ijm(u))

    finite = (u_ipj * du_ipj + u_imj * du_imj +
              v_ijp * du_ijp + v_ijm * du_ijm) / dx
    return finite


def advection_of_velocity_v(u, v, dx):
    """
    """
    """
    im   i    ip
    V    V    V
    P  U vv U P j
    V uv * uv V
    P  U vv U P jp
    V    V    V
    """
    # these names are relative to V, not to p
    # velocity terms
    v_ijp = (ijp(v) + v) / 2
    v_ijm = (ijm(v) + v) / 2
    u_ipj = (u + ijm(u)) / 2
    u_imj = (imj(u) + imjp(u)) / 2

    # gradient terms
    dv_ipj = (ipj(v) - v)
    dv_imj = (v - imj(v))
    dv_ijp = (ijp(v) - v)
    dv_ijm = (v - ijm(v))

    finite = (u_ipj * dv_ipj + u_imj * dv_imj +
              v_ijp * dv_ijp + v_ijm * dv_ijm) / dx
    return finite
    return NotImplemented
    # ipj is i+1/2,j ; imj is i-1/2,j
    u_ipj = (unit_roll(u, -1, 0) + u) / 2
    u_imj = (unit_roll(u, 1, 0) + u) / 2
    dv_ipj = unit_roll(v, -1, 0) - v
    dv_imj = v - unit_roll(v, 1, 0)
    v_ijp = (unit_roll(v, -1, 1) + v) / 2
    v_ijm = (unit_roll(u, 1, 1) + v) / 2
    dv_ijp = unit_roll(u, -1, 1) - v
    dv_ijm = v - unit_roll(v, 1, 1)

    finite = (u_ipj * dv_ipj + u_imj * dv_imj +
              v_ijp * dv_ijp + v_ijm * dv_ijm) / (2 * dx)
    return finite


def geopotential_gradient_u(p, dx):
    p_ipj = ipj(p)
    finite = (p_ipj - p) / dx * G
    return finite


def geopotential_gradient_v(p, dx):
    p_ijp = ijp(p)
    finite = (p_ijp - p) / dx * G
    return finite


def advection_of_geopotential(u, v, p, dx):
    u_imj = imj(u)
    v_ijm = ijm(v)
    up_imj = (imj(p) + p) / 2 * u_imj
    up_ipj = (ipj(p) + p) / 2 * u
    vp_ijm = (ijm(p) + p) / 2 * v_ijm
    vp_ijp = (ijp(p) + p) / 2 * v

    finite = (up_ipj - up_imj) / dx + (vp_ijp - vp_ijm) / dx
    return finite


def matsumo_scheme(u, v, p, dx, dt):
    u_star = u - dt * (advection_of_velocity_u(u, v, dx) +
                       geopotential_gradient_u(p, dx))
    # v_star = v
    v_star = v - dt * (advection_of_velocity_v(u, v, dx) +
                       geopotential_gradient_v(p, dx))
    p_star = p - dt * advection_of_geopotential(u, v, p, dx)

    geo_u_star = geopotential_gradient_u(p_star, dx)
    u_next = u - dt * (advection_of_velocity_u(u_star, v_star, dx) +
                       geo_u_star)
    # v_next = v
    v_next = v - dt * (advection_of_velocity_v(u_star, v_star, dx) +
                       geopotential_gradient_v(p_star, dx))
    pit_star = advection_of_geopotential(u_star, v_star, p_star, dx)
    p_next = p - dt * pit_star

    return u_next, v_next, p_next


def main():
    side_len = 64
    dx = 300 * units.km
    dt = 100 * units.s
    half = side_len // 2
    u = np.zeros((side_len, side_len)) * units.m / units.s
    v = np.zeros((side_len, side_len)) * units.m / units.s
    H = 8000 * units.m
    p = np.zeros((side_len, side_len), dtype=np.float) * units.m
    p[:] = H
    # p[half: half + 2, half:half+2] += 1 * units.m
    p[1, 2] += 1 * units.m
    u[half,half] += 1 * units.m / units.s

    plt.ion()
    plt.figure()
    plt.imshow(p)
    plt.title('Initial state')
    plt.show()

    initial_variation = get_total_variation(p)
    print("Initial Variation: %s" % (initial_variation,))

    for i in tqdm.tqdm(range(30000)):
        # print("iteration %s" % i)
        plt.clf()
        plt.imshow(p - H)
        plt.title('n = %s' % (i,))
        plt.show()
        plt.pause(0.001)  # pause a bit so that plots are updated

        u, v, p = matsumo_scheme(u, v, p, dx, dt)
        current_variation = get_total_variation(p)
        if current_variation.m > initial_variation.m + 0.1:
            print("iteration %s" % i)
            print("Variation too high: %s" % (current_variation,))
            # return False
            # break
    final_variation = get_total_variation(p)
    print("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

    plt.ioff()
    plt.show()


if __name__ == '__main__':
    main()


