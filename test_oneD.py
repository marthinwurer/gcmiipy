import logging
import numpy as np
import matplotlib.pyplot as plt
import unittest

from just_units import advect_1d_fd, advect_1d_lf, advect_1d_upwind, advect_1d_upwind_second, advect_1d_upwind_third, \
    upwind_with_spatial, ftcs_with_central, ft_with_upwind, lax_friedrichs, get_total_variation, WorldState, \
    get_initial_world
from constants import units

logger = logging.getLogger(__name__)

def setUpModule():
    FORMAT = '%(asctime)-15s | %(filename)s:%(lineno)s | %(message)s'
    logging.basicConfig(format=FORMAT, level=logging.DEBUG)


dx = 10 * units.m

dy = 10 * units.m

dt = 1 * units.s

world_shape = (160,)

half = world_shape[0] // 2

quarter = half // 2

spatial_change = (dx,)

num_steps = 400


class TestOneD(unittest.TestCase):

    def get_initial_conditions(self):
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s

        p = np.zeros((*world_shape,)) * units.Pa

        rho = np.zeros((*world_shape,)) * units.kg * units.m ** -3

        q = np.zeros((*world_shape,)) * units.kg / units.kg

        # initial conditions
        q[quarter:half] = 1.0
        V[0][:] = 1.0 * units.m / units.s

        return V, q

    def basic_run(self, func):
        V, q = self.get_initial_conditions()
        plt.ion()
        plt.figure()
        plt.plot(q)
        plt.title('Initial state')
        plt.show()
        q_prev = q

        initial_variation = get_total_variation(q)
        logger.info("Initial Variation: %s" % (initial_variation,))

        for i in range(num_steps):
            print("iteration %s" % i)
            plt.clf()
            plt.plot(q_prev)
            plt.title('current_state...')
            plt.show()
            plt.pause(0.001)  # pause a bit so that plots are updated

            q_next = func(dt, spatial_change, V, q_prev)
            q_prev = q_next
        else:
            q_next = q

        final_variation = get_total_variation(q_next)
        logger.info("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

        plt.ioff()
        plt.show()

    def world_run(self, func, world):
        plt.ion()
        plt.figure()
        plt.plot(world.q)
        plt.title('Initial state')
        plt.show()

        initial_variation = get_total_variation(world.q)
        logger.info("Initial Variation: %s" % (initial_variation,))

        for i in range(num_steps):
            print("iteration %s" % i)
            plt.clf()
            plt.plot(world.q)
            plt.title('current_state...')
            plt.show()
            plt.pause(0.001)  # pause a bit so that plots are updated

            world = func(dt, world)

        final_variation = get_total_variation(world.q)
        logger.info("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

        plt.ioff()
        plt.show()

    def test_forward(self):
        self.world_run(advect_1d_fd, get_initial_world(world_shape, spatial_change))

    def test_upwind(self):
        self.basic_run(advect_1d_upwind)

    def test_upwind_second(self):
        self.basic_run(advect_1d_upwind_second)

    def test_upwind_third(self):
        self.basic_run(advect_1d_upwind_third)

    def test_upwind_spatial(self):
        self.basic_run(upwind_with_spatial)

    def test_ftcs_central(self):
        self.basic_run(ftcs_with_central)

    def test_ftcs_upwind(self):
        self.basic_run(ft_with_upwind)

    def test_lax_friedrichs(self):
        self.basic_run(lax_friedrichs)

    def test_leapfrog(self):
        V, q = self.get_initial_conditions()
        plt.ion()
        plt.figure()
        q_prev = q
        q_next = advect_1d_fd(dt, spatial_change, V, q_prev)
        q_current = q_next

        for i in range(num_steps):
            print("iteration %s" % i)
            plt.clf()
            plt.plot(q_prev)
            plt.title('current_state...')
            plt.show()
            plt.pause(0.001)  # pause a bit so that plots are updated

            q_next = advect_1d_lf(dt, spatial_change, V, q_current, q_prev)
            q_prev = q_current
            q_current = q_next

        plt.ioff()
        plt.show()

