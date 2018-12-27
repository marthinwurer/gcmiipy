from itertools import count

import numpy as np
from constants import units

import matplotlib
import matplotlib.pyplot as plt

# just do conservation of momentum

# model constants

dx = 10000 * units.m

dy = 10000 * units.m

dt = 900 * units.s

world_shape = (160,)

spatial_change = (dx,)


# states

V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s

p = np.zeros((*world_shape,)) * units.Pa

rho = np.zeros((*world_shape,)) * units.kg * units.m ** -3

q = np.zeros((*world_shape,)) * units.kg / units.kg


# formula

def unit_roll(a, shift, axis=None):
    return np.roll(a, shift, axis=axis) * a.units


def advect_1d_fd(dt, spatial_change, V, q):
    dx = spatial_change[0]
    # have to use -1 because we're rolling the next index back to the current one
    q_p_1 = unit_roll(q, -1, 0)

    print(q)


    difference = (q_p_1 - q)
    print(difference)

    print(V[0])
    mult = difference * V[0]
    print(mult, "- Mult")

    step = dt / dx
    print(step, "- Step")
    finite = mult * step

    print(finite)
    new = q - finite
    print(new)
    return new







def advection_forward_differences(dt, spatial_change, V, q):
    dimensions = len(spatial_change)
    if dimensions == 1:
        advect_1d_fd(dt, spatial_change, V, q)
    # for dimension, delta in enumerate(spatial_change):


# initial conditions


q[3:7] = 1.0
V[0][:] = 1.0 * units.m / units.s
# V[0][6:7] = 2.0 * units.m / units.s


plt.ion()
plt.figure()
plt.plot(q)
plt.title('Initial state')
plt.show()



# forward step the model

qs = [q]

q_prev = q

done = False
for i in range(30):
    print("iteration %s" % i)
    plt.clf()
    plt.plot(q_prev)
    plt.title('current_state...')
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated

    q_next = advect_1d_fd(dt, spatial_change, V, q_prev)
    qs.append(q_next)
    q_prev = q_next

plt.ioff()
plt.show()

