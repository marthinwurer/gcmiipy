"""
Time to decode what GCMII does for advecting temperature:
"""
import math

import numpy as np

import constants

v_dim = 16
h_dim = 32

# DLAT=.5*TWOPI/(JM-1)
# d_lat is the latitude per grid cell. they had the -1 there for their polar cells, I think.
d_lat = math.pi / (v_dim)
print("d_lat", d_lat)

lat = np.zeros((v_dim))
lat[0] = -math.pi/2
lat[-1] = -lat[0]

# FJEQ is j at the equator.
fj_eq = (1 + v_dim) / 2
print("fj_eq", fj_eq)

for i in range(1, v_dim - 1):
    j = i + 1
    lat[i] = d_lat * (j - fj_eq)

print(lat)

deg_lat = lat / math.pi * 180
print(deg_lat)

# COSP and SINP are just cos and sin at those latitudes.
# DXP is more interesting. it's the radius times d_long times the cosine.
# I think it's the distance between cell centers.
# oh yeah it's dlong, of course it is
radius = constants.radius

# Is this just hungarian notation?
# P for at P, V for at V?
# damn it probably is.
# hell I was thinking of doing that too. Damn.

# Yeah, DXYP is the area of the cell at P.

# aight, the lines that I want to figure out:
# FD(I,J)=PA(I,J)*DXYP(J)
# TT(I,J,L)=TT(I,J,L)*FD(I,J)
# And at the end:
# FD(I,J)=PB(I,J)*DXYP(J)
# RFDP=1./FD(I,J)
# TT(I,J,L)=TT(I,J,L)*RFDP
#
# Looks like they do the same thing with water content.

# PA is pressure, I'm pretty sure. : 'kilogram / meter / second ** 2'
# DXYP is area. : meter ** 2
# FD ends up being: 'kilogram * meter / second ** 2'
# That's almost joules, which are 'kilogram * meter ** 2 / second ** 2'

# oh, then they multiply it by the mass flux, not just by the velocity
# that gets calculated in:
# SUBROUTINE AFLUX (U,V,P)
# which modifies -> SPA, PU, PV, CONV, PIT, SD,

"""
Time to break down GCMII's main time step code:
      INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)
      CALL AFLUX (U,V,P)
      CALL ADVECM (P,PB,DTFS)
      CALL ADVECV (P,UX,VX,PB,U,V,P,DTFS)
      CALL ADVECT (P,TX,PB,T,DTFS)
      CALL PGF (UX,VX,PB,U,V,T,P,DTFS)

U: horizontal velocity
V: vertical velocity
P: pressure
PB: next pressure. Starts as a copy of P.
DTFS: Change in time for the forward step
UX: Copy of U
VX: A global of shape V

CALCULATES THE AIR MASS FLUXES
AFLUX(U. V, P) -> SPA, PU, PV, CONV, PIT, SD
U: horizontal velocity
V: vertical velocity
P: pressure
SPA: Used to calc PU. Gets smoothed.
PU: East-west mass flux
PV: North-south mass flux
CONV: Horizontal mass convergence
PIT: Pressure tendency
SD: Sigma dot

CALCULATES UPDATED COLUMN PRESSURES
ADVECM (P,PA,DT1): -> PA
PA: New surface pressure

ADVECTS MOMENTUM (INCLUDING THE CORIOLIS FORCE)
ADVECV (PA,UT,VT,PB,U,V,P,DT1): -> FD, UT, DUT
Uses PU and PV from AFLUX
FD: Scaling parameter. PA(I,J)*DXYP(J)
UT: Scaled UT, momentum component
DUT: change in UT
VT: Not scaled? 

ADVECTS POTENTIAL TEMPERATURE
ADVECT (PA,TT,PB,T,DT1) ->
Uses PU and PV from AFLUX
They are used to advect the scaled potential temperature.
FD: Scaling parameter. Scales with PA then unscales with PB


"""





