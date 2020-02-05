"""
Using the matsumo scheme with temperature
"""
import matplotlib.pyplot as plt
import tqdm

from constants import *
from matsuno_c_grid import advection_of_velocity_u, advection_of_velocity_v, geopotential_gradient_u,\
    geopotential_gradient_v, advection_of_geopotential


def density_from(p, t):
    # compute real temperature from potential temperature
    pressure_ratio = (100000 * units.Pa / p)
    temp = t / (pressure_ratio ** (Rd / Cp))
    # rho = P/(RT)
    rho = p / (Rd * temp)
    return rho.to_base_units()


def geopotential_from(rho, p):
    geo = p / (G * rho)
    return geo.to_base_units()


def matsumo_scheme(u, v, p, t, dx, dt):
    density = density_from(p, t)
    geo = geopotential_from(density, p)
    u_star = u - dt * (advection_of_velocity_u(u, v, dx) +
                       geopotential_gradient_u(geo, dx))
    # v_star = v
    v_star = v - dt * (advection_of_velocity_v(u, v, dx) +
                       geopotential_gradient_v(geo, dx))
    p_star = p - dt * advection_of_geopotential(u, v, p, dx)
    t_star = t - dt * advection_of_geopotential(u, v, t, dx)

    density_star = density_from(p_star, t_star)
    geo_star = geopotential_from(density_star, p_star)
    geo_u_star = geopotential_gradient_u(geo_star, dx)
    u_next = u - dt * (advection_of_velocity_u(u_star, v_star, dx) +
                       geo_u_star)
    # v_next = v
    geo_v_star = geopotential_gradient_v(geo_star, dx)
    v_next = v - dt * (advection_of_velocity_v(u_star, v_star, dx) +
                       geo_v_star)
    pit_star = advection_of_geopotential(u_star, v_star, p_star, dx)
    p_next = p - dt * pit_star
    t_next = t - dt * advection_of_geopotential(u, v, t_star, dx)

    return u_next, v_next, p_next, t_next


def main():
    side_len = 64
    dx = 300 * units.km
    dt = 900 * units.s
    half = side_len // 2
    u = np.zeros((side_len, side_len)) * units.m / units.s
    v = np.zeros((side_len, side_len)) * units.m / units.s
    t = np.full((side_len, side_len), standard_temperature.m) * standard_temperature.u
    p = np.zeros((side_len, side_len), dtype=np.float) * standard_pressure.u
    H = standard_pressure
    p[:] = H
    # p[half: half + 2, half:half+2] += 1 * units.m
    p[1, 2] += 1 * p.u
    u[half, half] += 1 * units.m / units.s

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

        u, v, p, t = matsumo_scheme(u, v, p, t, dx, dt)
        current_variation = get_total_variation(p)
        if current_variation.m > initial_variation.m + 0.1:
            print("iteration %s" % i)
            print("Variation too high: %s" % (current_variation,))
            # return False
            # break
        if np.isnan(p).any():
            print("iteration %s" % i)
            print("NaN encountered")
            break
    final_variation = get_total_variation(p)
    print("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

    plt.ioff()
    plt.show()


if __name__ == '__main__':
    main()

