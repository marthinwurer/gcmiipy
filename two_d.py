from functools import reduce
from operator import mul

import numpy as np
import matplotlib.pyplot as plt

from constants import units, unit_roll, unit_maximum, unit_minimum, unit_stack, R, Rd, kappa, P0, standard_temperature, \
    get_total_variation


def upwind_axis(dt, spatial_change, V, q, axis=0):
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
    new = q - finite
    return new


def upwind_axis_finite(dt, spatial_change, V, q, axis=0):
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
    return finite



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


def gradient(p, spatial_change, axis):
    p_p_1 = unit_roll(p, -1, axis)
    p_m_1 = unit_roll(p, 1, axis)
    return (p_p_1 - p_m_1) / (2 * spatial_change[axis])


def pressure_gradient(dt, spatial_change, p, t):
    """
    sigma * pi / rho Del pi
    """
    # gradient in x dimension
    x_grad = gradient(p, spatial_change, 0)
    y_grad = gradient(p, spatial_change, 1)

    grad = unit_stack([x_grad, y_grad])
    # SPA(I,J,LM)=SIG(LM)*SP*RGAS*T(I,J,LM)*PKDN/PDN
    # spa[y][x] = (p[y][x] + ptop)**kappa * t[y][x] * R_dry
    sig = 1
    true_t = t / (P0 / p) ** kappa
    # spa = p **
    rho = (p / (Rd * true_t)).to_base_units()

    # no pressure here because the original was for momentum
    # spa = p / rho

    output = grad / rho * dt
    return output


def fv_advect_axis_upwind(dt, spatial_change, V, p, axis=0):
    dx = spatial_change[axis]
    p_p_1 = unit_roll(p, -1, axis)

    zeroes = np.zeros(p.shape) * V.units  # need to cast zeros to unit type
    a_plus = unit_maximum(V[axis], zeroes)
    a_minus = unit_minimum(V[axis], zeroes)

    flux = (p * a_plus + p_p_1 * a_minus) * dt / dx

    f_m_1 = unit_roll(flux, 1, axis)

    p_next = p - flux + f_m_1
    return p_next


def fv_advect_axis_upwind_finite(dt, spatial_change, V, p, axis=0):
    dx = spatial_change[axis]
    p_p_1 = unit_roll(p, -1, axis)

    zeroes = np.zeros(p.shape) * V.units  # need to cast zeros to unit type
    a_plus = unit_maximum(V[axis], zeroes)
    a_minus = unit_minimum(V[axis], zeroes)

    flux = (p * a_plus + p_p_1 * a_minus) * dt / dx

    f_m_1 = unit_roll(flux, 1, axis)

    finite = f_m_1 - flux
    return finite


def fv_advect_axis_plain(dt, spatial_change, V, p, axis=0):
    dx = spatial_change[axis]
    volume = reduce(mul, [*spatial_change], 1)
    area = volume / dx
    p_p_1 = unit_roll(p, -1, axis)

    average_at_edge = (p + p_p_1) / 2

    # flux is the velocity times the quantity
    flux = V[axis] * average_at_edge * dt * area

    f_m_1 = unit_roll(flux, 1, axis)

    p_next = p - (flux - f_m_1) / volume
    return p_next


def fv_advect_axis_plain_finite(dt, spatial_change, V, p, axis=0):
    dx = spatial_change[axis]
    volume = reduce(mul, [*spatial_change], 1)
    area = volume / dx
    p_p_1 = unit_roll(p, -1, axis)

    average_at_edge = (p + p_p_1) / 2

    # flux is the velocity times the quantity
    flux = V[axis] * average_at_edge * dt * area

    f_m_1 = unit_roll(flux, 1, axis)

    finite = f_m_1 - flux
    return finite



"""
The discretization scheme will be a C grid.

    P,Q,T | V[0]
    ------+
    V[1]

You can split this into a grid of cells designated as A, B, C, and D.

A B A B A B
C D C D C D
A B A B A B
C D C D C D
A B A B A B
C D C D C D

P, Q, and T are stored at A cells. Velocities are stored at B and C for U and V respectively. 
Fluxes for U and V need to be calculated by averaging U and V values at the A and D cells.


    velocity is the outward flow at the south and east face.

    P,Q,T | V[0]
    ------+
    V[1]
"""


def finite_volume_advection(dt, spatial_change, V, p):
    """
    Uses dimensional splitting.
    """
    p_star = p

    for axis in range(2):
        p_star = fv_advect_axis_upwind(dt, spatial_change, V, p_star, axis)

    return p_star


def pgf_c_grid_axis(p, spatial_change, axis=0):
    # TODO rename this
    """
    Compute the pressure gradient along an axis
    """
    dx = spatial_change[axis]
    p_p_1 = unit_roll(p, -1, axis)

    pgf = (p_p_1 - p) / dx

    return pgf


def pgf_c_grid(dt, spatial_change, p, t):
    """
    sigma * pi / rho Del pi
    """
    # gradient in x dimension
    x_grad = pgf_c_grid_axis(p, spatial_change, 0)
    y_grad = pgf_c_grid_axis(p, spatial_change, 1)

    grad = unit_stack([x_grad, y_grad])
    # SPA(I,J,LM)=SIG(LM)*SP*RGAS*T(I,J,LM)*PKDN/PDN
    # spa[y][x] = (p[y][x] + ptop)**kappa * t[y][x] * R_dry
    sig = 1

    # convert from potential temperature to real temperature
    true_t = t / (P0 / p) ** kappa
    # spa = p **
    rho = (p / (Rd * true_t)).to_base_units()

    # no pressure here because the original was for momentum
    # spa = p / rho

    output = grad / rho * dt
    return output


def pgf_templess(dt, spatial_change, p):
    """
    Assumes dry air and standard temperature
    """
    # gradient in x dimension
    x_grad = pgf_c_grid_axis(p, spatial_change, 0)
    y_grad = pgf_c_grid_axis(p, spatial_change, 1)

    p_edge = pressure_at_edge(p)
    d_edge = p_edge / (Rd * standard_temperature)

    grad = unit_stack([x_grad, y_grad])
    output = grad * dt / d_edge
    return output


def pressure_at_edge(p):
    # axes in spatial change are (x, y)
    p_east = (unit_roll(p, -1, 0) + p) / 2
    p_south = (unit_roll(p, -1, 1) + p) / 2
    return unit_stack([p_east, p_south])


def pressure_at_edge_one_d(p):
    # axes in spatial change are (x, y)
    p_east = (unit_roll(p, -1, 0) + p) / 2
    return p_east


def advect_with_momentum(dt, spatial_change, V, p):
    p_edge = pressure_at_edge(p)
    momentum = V * p_edge

    # I compute flux at a different point.
    # I need to split out directional flux computation to it's own function
    # Then I can use that flux in all the other functions
    # or I just need to divide the pgf by the mass or something
    # PV, THE SOUTH-NORTH MASS FLUX
    # PV(I,J,L)=.25*DXV(J)*(V(I,J,L)+V(IM1,J,L))*(P(I,J)+P(I,J-1))
    # CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))*DSIG(L)
    # PIT(I,J)=PIT(I,J)+CONV(I,J,L)
    # PA(I,J)=P(I,J)+(DT1*PIT(I,J)/DXYP(J)+PTRUNC)
    p_next = finite_volume_advection(dt, spatial_change, momentum, p)

    return p_next


def pgf_one_d(dt, dx, p, axis=0):
    p_p_1 = unit_roll(p, -1, axis)
    pressure_gradient = (p_p_1 - p) / dx

    d_edge = pressure_at_edge_one_d(p) / (Rd * standard_temperature)

    pgf = pressure_gradient * dt / d_edge

    return pgf


def run_2d_with_ft(initial_conditions, ft, steps=400, display_key="q", variation_key="q"):
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
    plt.imshow(current_conditions[display_key])
    plt.title('Initial state')
    plt.show()

    initial_variation = get_total_variation(current_conditions[variation_key])
    print("Initial Variation: %s" % (initial_variation,))

    for i in range(steps):
        # print("iteration %s" % i)
        plt.clf()
        plt.imshow(current_conditions[display_key])
        plt.title('current_state...')
        plt.show()
        plt.pause(0.001)  # pause a bit so that plots are updated

        current_conditions = ft(**current_conditions)
        current_variation = get_total_variation(current_conditions[variation_key])
        if current_variation.m > initial_variation.m + 0.00001:
            print("iteration %s" % i)
            print("Variation too high: %s" % (current_variation,))
            # return False

    final_variation = get_total_variation(current_conditions[variation_key])
    print("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

    plt.ioff()
    plt.show()

    return True






