import numpy as np
from constants import units
from two_d import corner_transport_2d

dx = 10 * units.m
dy = 10 * units.m
dt = 1 * units.s
side_length = 16
world_shape = (side_length, side_length * 2)
half = world_shape[0] // 2
quarter = half // 2
spatial_change = (dx, dy)
num_steps = 400


def get_initial_conditions():
    V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s

    p = np.zeros((*world_shape,)) * units.Pa

    rho = np.zeros((*world_shape,)) * units.kg * units.m ** -3

    q = np.zeros((*world_shape,)) * units.grams

    t = np.full((*world_shape,), 273.15) * units.Kelvin

    # initial conditions
    q[quarter:half, quarter:half] = 1.0 * units.grams
    # V[0][:] = -1.0 * units.m / units.s
    # V[1][:] = -1.0 * units.m / units.s
    V[0][half] = 1.0 * units.m / units.s

    return V, q, p, rho, t

def main():

    """
    old fortran code:
      CALL AFLUX (UT,VT,PA)  # calculates fluxes of air mass
      CALL ADVECM (PC,P,DTLF)  # advect pressure
      CALL ADVECV (PC,U,V,P,UT,VT,PA,DTLF)  # advect velocity
      CALL ADVECT (PC,T,P,TT,DTLF)  # advect temperature
      CALL ADVECQ (PC,Q,P,QT,DTLF)  # advect moisture
      CALL PGF (U,V,P,UT,VT,TT,PA,DTLF)  # calculate the pressure gradient force
    """

    V, q, p, rho, t = get_initial_conditions()

    p_next = corner_transport_2d(dt, spatial_change, V, p)
    Vx_next = corner_transport_2d(dt, spatial_change, V, V[0])
    Vy_next = corner_transport_2d(dt, spatial_change, V, V[1])
    t_next = corner_transport_2d(dt, spatial_change, V, t)
    q_next = corner_transport_2d(dt, spatial_change, V, q)


if __name__ == '__main__':
    main()
