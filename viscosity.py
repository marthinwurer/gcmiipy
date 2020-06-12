"""
Starting based on an idea I got from the disqus conversation.
My system is closed, and there are numerical errors that can add up.
I need some form of entropy to drive the energy in the system down.
Viscosity is part of the Navier-Stokes equations and can diffuse energy.
I want to implement it so I don't need filtering.
I want the system to damp itself instead of using external forcings.
"""
from coordinates import ijp, ijm, ipj, imj, imjp


def finite_laplacian_2d(q, dx):
    """
    Use the five point stencil to compute the laplacian of q
    """
    top = ijp(q) + ijm(q) + ipj(q) + imj(q) - 4 * q
    output = top / (dx * dx)

    return output


def incompressible_viscosity_2d(u, mu, dx):
    lap = finite_laplacian_2d(u, dx)
    output = mu * lap
    return output.to_base_units()





