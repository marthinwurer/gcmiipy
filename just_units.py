import numpy as np

from constants import units, unit_roll

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

    zeroes = np.zeros(q.shape)

    a_plus = np.maximum(V[0], zeroes)
    a_minus = np.minimum(V[0], zeroes)

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


def get_total_variation(q):
    q_p_1 = unit_roll(q, -1, 0)
    diff = q - q_p_1
    return np.sum(np.abs(diff.m)) * q.u

