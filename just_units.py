import numpy as np


# formula

def unit_roll(a, shift, axis=None):
    return np.roll(a, shift, axis=axis) * a.units


def advect_1d_fd(dt, spatial_change, V, q):
    dx = spatial_change[0]
    # have to use -1 because we're rolling the next index back to the current one
    q_p_1 = unit_roll(q, -1, 0)

    print(q)

    difference = (q_p_1 - q)
    print(difference)

    print(V[0])
    mult = difference * V[0]
    print(mult, "- Mult")

    step = dt / dx
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = q - finite
    print(new)
    return new


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

    zeroes = np.zeros(q.shape)

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


def advection_forward_differences(dt, spatial_change, V, q):
    dimensions = len(spatial_change)
    if dimensions == 1:
        advect_1d_fd(dt, spatial_change, V, q)
    # for dimension, delta in enumerate(spatial_change):


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
    return np.sum(np.abs(diff))

