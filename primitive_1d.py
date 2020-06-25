"""
Discretizing the full primitive equations on a 1d staggered grid

grid is:
   i h ip
j  P U P

U is effectively at iph
"""
from coordinates_1d import *
from constants import *


def advect_v_u(u, dx):
    # get staggered velocity
    u_iph = iph(u)
    u_imh = imh(u)

    du_ip = (ip(u) - u)
    du_im = (u - im(u))

    finite = (u_iph * du_ip + u_imh * du_im) / dx
    return finite




def advect_forward_euler(rho, u, dx, dt):
    rho_next = rho - advect_rho(rho, u, dx) * dt
    u_next = u
    return rho_next, u_next




# need to try formula at end of https://climate.ucdavis.edu/AOSS605-NumericalMethodsLectures.pdf



# ok, flux form should work great. input - output.
# flux form of drho/dt
def advect_rho(rho, u, dx):
    rho_ph = iph(rho)

    urho = u * rho_ph

    finite = (urho - im(urho)) / dx
    return finite


def advect_u_scaled(ut, u, p, pa, dt, dx):
    # TODO fix this
    ut_scaled = iph(p) * ut

    u_at_h = imh(u)

    adv_val = u_at_h * u_at_h * p

    # todo this was backwards
    adv_diff = (adv_val - ip(adv_val)) / dx

    geo_val = p * p * G / 2

    geo_diff = (geo_val - ip(geo_val)) / dx

    output = ((ut * p) - (adv_diff + geo_diff) * dt) / pa

    return output


def advect_matsumo(rho, u, dt, dx):
    rho_star = rho - advect_rho(rho, u, dx) * dt
    u_star = u
    rho_next = rho - advect_rho(rho_star, u_star, dx) * dt
    u_next = u

    return rho_next, u_next


def shallow_water_matsuno(h, u, dt, dx):
    h_star = h - advect_rho(h, u, dx) * dt
    u_star = advect_u_scaled(u, u, h, h_star, dt, dx)
    u_star[-1] = 0 * u.u
    h_next = h - advect_rho(h_star, u_star, dx) * dt
    u_next = advect_u_scaled(u, u_star, h_star, h_next, dt, dx)
    u_next[-1] = 0 * u.u

    return h_next, u_next


def advect_maccormack(rho, u, dt, dx):
    rho_star = rho - advect_rho(rho, u, dx) * dt
    u_star = u
    rho_next = ((rho + rho_star) - advect_rho(rho_star, u_star, dx) * dt) / 2
    u_next = u

    return rho_next, u_next


def shallow_water_maccormack(h, u, dt, dx):
    # https://en.wikipedia.org/wiki/MacCormack_method
    pass


def advect_lax_friedrichs(rho, u, dt, dx):
    rho_next = (ip(rho) + im(rho)) / 2 - advect_rho(rho, u, dx) * dt
    u_next = u

    return rho_next, u_next


from flux_limiter import donor_cell_flux, donor_cell_advection

def advect_upwind(rho, u, dt, dx):
    flux = donor_cell_flux(rho, u)
    #todo finish this

















