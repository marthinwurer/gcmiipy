import numpy as np
import matplotlib.pyplot as plt

from constants import units, unit_roll, get_total_variation, unit_maximum, unit_minimum, G


# world state


class WorldState(object):
    def __init__(self, spatial_change, V, p, rho, q):

        self.spatial_change = spatial_change

        self.V = V
        self.p = p
        self.rho = rho
        self.q = q

    def q_x(self, x):
        # have to use -x because we're rolling the selected index to the
        # current one, in the other direction
        return unit_roll(self.q, -x, 0)

    def dx(self):
        return self.spatial_change[-1]

    def v_x(self):
        return self.V[-1]

    def next_q(self, q_next):
        return WorldState(self.spatial_change, self.V, self.p, self.rho, q_next)


def get_initial_world(world_shape, spatial_change):
    V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
    p = np.zeros((*world_shape,)) * units.Pa
    rho = np.zeros((*world_shape,)) * units.kg * units.m ** -3
    q = np.zeros((*world_shape,)) * units.kg / units.kg

    # getting the dimensions
    half = world_shape[0] // 2
    quarter = half // 2

    # initial conditions
    q[quarter:half] = 1.0
    V[0][:] = 1.0 * units.m / units.s
    world = WorldState(spatial_change, V, p, rho, q)
    return world

# formulas


def advect_1d_fd(dt, world: WorldState):
    dx = world.dx()
    # have to use -1 because we're rolling the next index back to the current one
    q_p_1 = world.q_x(1)

    print(world.q)

    difference = (q_p_1 - world.q)
    print(difference)

    print(world.V[0])
    mult = difference * world.V[0]
    print(mult, "- Mult")

    step = dt / dx
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = world.q - finite
    print(new)
    return world.next_q(new)


def advect_1d_lf(dt, spatial_change, V, q, q_prev):
    dx = spatial_change[0]
    # have to use -1 because we're rolling the next index back to the current one
    q_p_1 = unit_roll(q, -1, 0)
    q_m_1 = unit_roll(q, 1, 0)

    difference = (q_p_1 - q_m_1)

    mult = difference * V[0]
    print(mult, "- Mult")

    step = (dt / dx) * 2
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = q_prev - finite
    print(new)
    return new


def upwind_spatial(spatial_change, V, q):
    # break out for method of lines: do spatial derivative, then solve for time.
    dx = spatial_change[0]

    q_p_1 = unit_roll(q, -1, 0)
    q_m_1 = unit_roll(q, 1, 0)

    zeroes = np.zeros(q.shape) * V.units  # need to cast zeros to unit type

    a_plus = unit_maximum(V[0], zeroes)
    a_minus = unit_minimum(V[0], zeroes)

    fd = (q_p_1 - q)
    bd = (q - q_m_1)

    mult = (fd * a_minus + bd * a_plus) * V[0] / dx
    print(mult, "- Mult")

    return mult


def upwind_with_spatial(dt, spatial_change, V, q):
    mult = upwind_spatial(spatial_change, V, q)
    finite = mult * dt

    print(finite)
    new = q - finite
    print(new)
    return new


def advect_1d_upwind(dt, spatial_change, V, q):
    dx = spatial_change[0]

    q_p_1 = unit_roll(q, -1, 0)
    q_m_1 = unit_roll(q, 1, 0)

    zeroes = np.zeros(q.shape) * units.meter / units.second

    a_plus = np.maximum(V[0], zeroes)  # need to cast zeros to unit tyoe
    a_minus = np.minimum(V[0], zeroes)

    fd = (q_p_1 - q)
    bd = (q - q_m_1)

    mult = (fd * a_minus + bd * a_plus) * V[0]
    print(mult, "- Mult")

    step = (dt / dx)
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = q - finite
    print(new)
    return new


def advect_1d_upwind_second(dt, spatial_change, V, q):
    dx = spatial_change[0]

    q_p_1 = unit_roll(q, -1, 0)
    q_p_2 = unit_roll(q, -2, 0)
    q_m_1 = unit_roll(q, 1, 0)
    q_m_2 = unit_roll(q, 2, 0)

    zeroes = np.zeros(q.shape)

    a_plus = np.maximum(V[0], zeroes)
    a_minus = np.minimum(V[0], zeroes)

    fd = (4 * q_p_1 - 3 * q - q_p_2)
    bd = (3 * q - 4 * q_m_1 + q_m_2)

    mult = (fd * a_minus + bd * a_plus) * V[0]
    print(mult, "- Mult")

    step = (dt / dx) / 2
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = q - finite
    print(new)
    return new


def advect_1d_upwind_third(dt, spatial_change, V, q):
    dx = spatial_change[0]

    q_p_1 = unit_roll(q, -1, 0)
    q_p_2 = unit_roll(q, -2, 0)
    q_m_1 = unit_roll(q, 1, 0)
    q_m_2 = unit_roll(q, 2, 0)

    zeroes = np.zeros(q.shape)

    a_plus = np.maximum(V[0], zeroes)
    a_minus = np.minimum(V[0], zeroes)

    bd = (2 * q_p_1 + 3 * q - 6 * q_m_1 + q_m_2)
    fd = (6 * q_p_1 - 3 * q - q_p_2 - 2 * q_m_1)

    mult = (fd * a_minus + bd * a_plus) * V[0]
    print(mult, "- Mult")

    step = (dt / dx) / 6
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = q - finite
    print(new)
    return new


def advection_1d_rk4(dt, spatial_change, V, q):
    dx = spatial_change[0]

    q_p_1 = unit_roll(q, -1, 0)
    q_m_1 = unit_roll(q, 1, 0)

    zeroes = np.zeros(q.shape)

    # TODO This gets the wrong units
    a_plus = np.maximum(V[0], zeroes)
    a_minus = np.minimum(V[0], zeroes)

    fd = (q_p_1 - q)
    bd = (q - q_m_1)

    mult = (fd * a_minus + bd * a_plus) * V[0]
    print(mult, "- Mult")

    step = (dt / dx)
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = q - finite
    print(new)
    return new


def central_spatial(spatial_change, V, q):
    # break out for method of lines: do spatial derivative, then solve for time.
    dx = spatial_change[0]

    q_p_1 = unit_roll(q, -1, 0)
    q_m_1 = unit_roll(q, 1, 0)

    cd = (q_p_1 - q_m_1)

    mult = cd * V[0] / dx
    print(mult, "- Mult")

    return mult


def forward_time(dt, spatial_change, V, q, spatial_func):
    mult = spatial_func(spatial_change, V, q)
    finite = mult * dt

    print(finite)
    new = q - finite
    print(new)
    return new


def ftcs_with_central(dt, spatial_change, V, q):
    return forward_time(dt, spatial_change, V, q, central_spatial)


def ft_with_upwind(dt, spatial_change, V, q):
    return forward_time(dt, spatial_change, V, q, upwind_spatial)


def lax_friedrichs(dt, spatial_change, V, q):

    # central spatial difference with averaged center
    dx = spatial_change[0]

    q_p_1 = unit_roll(q, -1, 0)
    q_m_1 = unit_roll(q, 1, 0)

    cd = (q_p_1 - q_m_1)

    mult = cd * V[0] / dx
    print(mult, "- Mult")

    q_avg = (q_p_1 + q_m_1) / 2
    finite = mult * dt

    print(finite)
    new = q_avg - finite
    print(new)
    return new


def run_1d_with_ft(initial_conditions, ft, steps=400, display_key="q", variation_key="q"):
    """
    Run and display a 2d system
    Args:
        initial_conditions: the initial conditions of the system. A dict.
        ft: the function that computes the values at the next time step
        display_key: The key of the dict to graph
        variation_key: the key of the dict to compute variation
    """
    current_conditions = initial_conditions
    plt.ion()
    plt.figure()
    plt.plot(current_conditions[display_key])
    plt.title('Initial state')
    plt.show()

    initial_variation = get_total_variation(current_conditions[variation_key])
    print("Initial Variation: %s" % (initial_variation,))

    for i in range(steps):
        # print("iteration %s" % i)
        plt.clf()
        plt.plot(current_conditions[display_key])
        plt.title('current_state...')
        plt.show()
        plt.pause(0.001)  # pause a bit so that plots are updated

        current_conditions = ft(**current_conditions)
        current_variation = get_total_variation(current_conditions[variation_key])
        # if initial_variation.m + 0.00001 < current_variation.m or \
        if initial_variation.m + 1000 < current_variation.m or \
                np.isnan(current_conditions[variation_key]).any():
            print("iteration %s" % i)
            print("Variation too high: %s" % (current_variation,))
            return False

    final_variation = get_total_variation(current_conditions[variation_key])
    print("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

    plt.ioff()
    plt.show()

    return True


def shallow_water_g_center_space(dt, spatial_change, h):
    dx = spatial_change[0]

    h_p_1 = unit_roll(h, -1, 0)
    h_m_1 = unit_roll(h, 1, 0)

    cd = (h_p_1 - h_m_1) / (2 * dx) * G * dt
    return cd


def shallow_water_h_center_space(dt, spatial_change, u, H):
    dx = spatial_change[-1]
    u = u[-1]

    u_p_1 = unit_roll(u, -1, -1)
    u_m_1 = unit_roll(u, 1, -1)

    cd = (u_p_1 - u_m_1) / (2 * dx) * H * dt
    return cd


def shallow_water_g_c_grid(dt, spatial_change, h):
    dx = spatial_change[0]

    h_p_1 = unit_roll(h, -1, 0)

    cd = (h_p_1 - h) / dx * G * dt
    return cd


def shallow_water_h_c_grid(dt, spatial_change, u, H):
    dx = spatial_change[-1]

    # each u space is u[j+1/2]
    u_p_1 = u[-1]

    # now u[j-1/2]
    u_m_1 = unit_roll(u[-1], 1, -1)

    cd = (u_p_1 - u_m_1) / dx * H * dt
    return cd




