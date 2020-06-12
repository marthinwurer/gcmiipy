from matsumo_temp import *

side_len = 5
u, v, p, t = gen_initial_conditions(side_len)
dx = 300 * units.km
dt = 900 * units.s
half = side_len // 2
H = standard_pressure
# p[half: half + 2, half:half+2] += 1 * units.m
# p[1, 2] += 1 * p.u
u[half, half] += 1 * units.m / units.s
