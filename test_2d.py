import logging
import unittest

import numpy as np
import matplotlib.pyplot as plt

from constants import units
from just_units import get_total_variation
from two_d import *

logger = logging.getLogger(__name__)

def setUpModule():
    FORMAT = '%(asctime)-15s | %(filename)s:%(lineno)s | %(message)s'
    logging.basicConfig(format=FORMAT, level=logging.INFO)

dx = 10 * units.m
dy = 10 * units.m
dt = 1 * units.s
side_length = 160
world_shape = (side_length, side_length)
half = world_shape[0] // 2
quarter = half // 2
spatial_change = (dx, dy)
num_steps = 400


def get_initial_conditions():
    V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s

    p = np.zeros((*world_shape,)) * units.Pa

    rho = np.zeros((*world_shape,)) * units.kg * units.m ** -3

    q = np.zeros((*world_shape,)) * units.grams

    # initial conditions
    q[quarter:half, quarter:half] = 1.0 * units.grams
    V[0][:] = -1.0 * units.m / units.s
    V[1][:] = -1.0 * units.m / units.s

    return V, q


class MyTestCase(unittest.TestCase):
    def test_ctu(self):
        V, q = get_initial_conditions()
        plt.ion()
        plt.figure()
        plt.imshow(q)
        plt.title('Initial state')
        plt.show()
        q_prev = q

        initial_variation = get_total_variation(q)
        logger.info("Initial Variation: %s" % (initial_variation,))

        for i in range(num_steps):
            print("iteration %s" % i)
            plt.clf()
            plt.imshow(q_prev)
            plt.title('current_state...')
            plt.show()
            plt.pause(0.001)  # pause a bit so that plots are updated

            q_next = corner_transport_2d(dt, spatial_change, V, q_prev)
            q_prev = q_next
            current_variation = get_total_variation(q_prev)
            if current_variation.m > initial_variation.m + 0.00001:
                print("Variation too high: %s" % (current_variation,))
                self.assertEqual(True, False)

        else:
            q_next = q

        final_variation = get_total_variation(q_next)
        logger.info("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

        plt.ioff()
        plt.show()


if __name__ == '__main__':
    unittest.main()
