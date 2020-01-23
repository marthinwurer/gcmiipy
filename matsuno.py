
"""
The goal for this file is to implement Matsuno 1966.
This paper has the dynamics of the primitive equations discretized.
https://www.jstage.jst.go.jp/article/jmsj1965/44/1/44_1_76/_pdf
"""
import numpy as np
import matplotlib.pyplot as plt

"""
U is x dimension velocity
V is y dimension velocity
arrays are [y,x]
           [i, j]

need separate functions for u and for v
"""

from constants import units, unit_roll, G, get_total_variation


def advection_of_velocity_u(u, v, dx):
    # ipj is i+1/2,j ; imj is i-1/2,j
    u_ipj = (unit_roll(u, -1, 0) + u) / 2
    u_imj = (unit_roll(u, 1, 0) + u) / 2
    du_ipj = unit_roll(u, -1, 0) - u
    du_imj = u - unit_roll(u, 1, 0)
    v_ijp = (unit_roll(v, -1, 1) + v) / 2
    v_ijm = (unit_roll(u, 1, 1) + v) / 2
    du_ijp = unit_roll(u, -1, 1) - u
    du_ijm = u - unit_roll(u, 1, 1)

    finite = (u_ipj * du_ipj + u_imj * du_imj +
              v_ijp * du_ijp + v_ijm * du_ijm) / (2 * dx)
    return finite


def geopotential_gradient_u(p, dx):
    # ipj = i+1,j; m is minus 1
    p_ipj = unit_roll(p, -1, 0)
    p_imj = unit_roll(p, 1, 0)

    finite = (p_ipj - p_imj) / (2 * dx) * G
    return finite


def advection_of_geopotential(u, v, p, dx):
    up = u * p
    vp = v * p
    up_ipj = unit_roll(up, -1, 0)
    up_imj = unit_roll(up, 1, 0)
    vp_ijp = unit_roll(vp, -1, 1)
    vp_ijm = unit_roll(vp, 1, 1)

    finite = (up_ipj - up_imj) / dx + (vp_ijp - vp_ijm) / dx
    return finite


def matsumo_scheme(u, v, p, dx, dt):
    u_star = u - dt * (advection_of_velocity_u(u, v, dx) +
                        geopotential_gradient_u(p, dx))
    v_star = v
    # v_star = v - dt * (advection_of_velocity_v(u, v, dx) +
    #                    geopotential_gradient_v(p, dx))
    p_star = p - dt * advection_of_geopotential(u, v, p, dx)

    u_next = u - dt * (advection_of_velocity_u(u_star, v_star, dx) +
                       geopotential_gradient_u(p_star, dx))
    v_next = v
    # v_next = v - dt * (advection_of_velocity_v(u_star, v_star, dx) +
    #                    geopotential_gradient_v(p_star, dx))
    p_next = p - dt * advection_of_geopotential(u_star, v_star, p_star, dx)

    return u_next, v_next, p_next


def main():
    side_len = 16
    dx = 300 * units.km
    dt = 900 * units.s
    half = side_len // 2
    u = np.zeros((side_len, side_len)) * units.m / units.s
    v = np.zeros((side_len, side_len)) * units.m / units.s
    H = 1000 * units.m
    p = np.zeros((side_len, side_len), dtype=np.float) * units.m
    p[:] = H
    p[half: half + 3, half:half+3] += 1 * units.m
    p[1, 2] += 1 * units.m

    plt.ion()
    plt.figure()
    plt.imshow(p)
    plt.title('Initial state')
    plt.show()

    initial_variation = get_total_variation(p)
    print("Initial Variation: %s" % (initial_variation,))

    for i in range(400):
        # print("iteration %s" % i)
        plt.clf()
        plt.imshow(p)
        plt.title('current_state...')
        plt.show()
        plt.pause(0.001)  # pause a bit so that plots are updated

        u, v, p = matsumo_scheme(u, v, p, dx, dt)
        current_variation = get_total_variation(p)
        if current_variation.m > initial_variation.m + 0.00001:
            print("iteration %s" % i)
            print("Variation too high: %s" % (current_variation,))
            # return False

    final_variation = get_total_variation(p)
    print("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

    plt.ioff()
    plt.show()



if __name__ == '__main__':
    main()


