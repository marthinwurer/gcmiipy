"""
Holds all the computation relevant to a 1-d coordinate system

"""
from constants import unit_roll

"""
U is x dimension velocity, with the i component
arrays are [x]
           [i]

h = half
m = minus one
p = plus one

Coordinates:
Pij  Uij  Pijp

grid is:
   i h ip
j  P U P
"""


def ip(q):
    return unit_roll(q, -1, 0)


def im(q):
    return unit_roll(q, 1, 0)


def iph(q):
    return (q + ip(q)) / 2


def imh(q):
    return (q + im(q)) / 2


def div(q_h, dx):
    return (q_h - im(q_h)) / dx


def divu(q_h, dx):
    return (ip(q_h) - im(q_h)) / (2 * dx)


def gradh(q_i, dx):
    """
    Gradient at h of q_i
    """
    return (ip(q_i) - q_i) / dx
