"""
Holds all the computation relevant to coordinate systems.

"""
from constants import unit_roll

"""
U is x dimension velocity, with the i component
V is y dimension velocity, with the j component
arrays are [y, x]
           [j, i]

h = half
m = minus one
p = plus one

Coordinates:
Pij  Uij  Pijp
Vij       Vijp
Pipj Uipj Pipjp

grid is:
   i h ip
j  P U P
h  V   V
jp P U P
"""

def ipj(q):
    """
    i+1,j
    """
    return unit_roll(q, -1, 1)


def imj(q):
    return unit_roll(q, 1, 1)


def ijp(q):
    return unit_roll(q, -1, 0)


def ijm(q):
    return unit_roll(q, 1, 0)


def imjp(q):
    return imj(ijp(q))


def iph(q):
    return (q + ipj(q)) / 2


def imh(q):
    return (q + imj(q)) / 2


def jph(q):
    return (q + ijp(q)) / 2


def jmh(q):
    return (q + ijm(q)) / 2


def gradi(q_i, dx):
    """
    Gradient at h of q_i
    """
    return (ipj(q_i) - q_i) / dx


def gradj(q_j, dy):
    """
    Gradient at h of q_i
    """
    return (ijp(q_j) - q_j) / dy
