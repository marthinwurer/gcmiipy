import numpy as np

from constants import units, unit_roll, unit_maximum, unit_minimum


def advect_nd_upwind(dt, spatial_change, V, q, dimensions=None):
    if not dimensions:
        # figure out dimensions
        pass

    # V needs to be of the shape dims,

    for axis in range(dimensions):
        dx = spatial_change[axis]

        dx = spatial_change[0]


        q_p_1 = unit_roll(q, -1, 0)
        q_m_1 = unit_roll(q, 1, 0)

        zeroes = np.zeros(q.shape) * units.meter / units.second

        a_plus = np.maximum(V[0], zeroes)  # need to cast zeros to unit tyoe
        a_minus = np.minimum(V[0], zeroes)

        fd = (q_p_1 - q)
        bd = (q - q_m_1)

        # mult = (fd * a_minus + bd * a_plus) * V[0]
        mult = (fd * a_minus + bd * a_plus) * V[0]
        print(mult, "- Mult")

        step = (dt / dx)
        print(step, "- Step")
        finite = mult * step

    print(finite)
    new = q - finite
    print(new)
    return new


def upwind_axis(dt, spatial_change, V, q, axis):
    """
    upwind along a single axis
    https://en.wikipedia.org/wiki/Upwind_scheme
    """
    dx = spatial_change[axis]
    zeroes = np.zeros(q.shape) * V.units  # need to cast zeros to unit type

    q_p_1 = unit_roll(q, -1, axis)
    q_m_1 = unit_roll(q, 1, axis)

    a_plus = unit_maximum(V[axis], zeroes)
    a_minus = unit_minimum(V[axis], zeroes)

    u_minus = (q - q_m_1)
    u_plus = (q_p_1 - q)

    mult = a_plus * u_minus + a_minus * u_plus
    step = (dt / dx)
    finite = mult * step
    print(finite.units)
    new = q - finite
    return new


def corner_transport_2d(dt, spatial_change, V, q):
    """
    Use dimensional splitting to implement the CTU advection scheme.
    https://ocw.mit.edu/courses/earth-atmospheric-and-planetary-sciences/12-950-atmospheric-and-oceanic-modeling-spring-2004/lecture-notes/lec10.pdf

    """
    # V[axis][
    q_star = q

    for axis in range(2):
        q_star = upwind_axis(dt, spatial_change, V, q_star, axis)

    return q_star

