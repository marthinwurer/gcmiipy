
import numpy as np
from constants import units


# just do conservation of momentum


dx = 1000 * units.km

dy = 1000 * units.km


world_shape = (16,)

spatial_change = (dx,)

V = np.zeros((*world_shape, len(world_shape))) * units.m / units.s

p = np.zeros((*world_shape,)) * units.Pa

rho = np.zeros((*world_shape,)) * units.kg * units.m ** -3


