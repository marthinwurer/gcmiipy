import logging
import numpy as np
import matplotlib.pyplot as plt
import unittest

from just_units import *
from constants import units, get_total_variation, standard_pressure
from two_d import upwind_axis, fv_advect_axis_plain, fv_advect_axis_upwind, upwind_axis_finite, pgf_one_d, \
    fv_advect_axis_upwind_finite

logger = logging.getLogger(__name__)

def setUpModule():
    FORMAT = '%(asctime)-15s | %(filename)s:%(lineno)s | %(message)s'
    logging.basicConfig(format=FORMAT, level=logging.INFO)


dx = 10 * units.m

dy = 10 * units.m

dt = 1 * units.s

world_shape = (161,)

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
        V[0][:] = -1.0 * units.m / units.s
        # q[:] = 1.0
        # q[half] = .5
        # V[0][quarter:half] = 1.0 * units.m / units.s

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

    def test_upwind_axis(self):
        self.basic_run(upwind_axis)

    def test_fv_upwind(self):
        self.basic_run(fv_advect_axis_upwind)

    def test_fv_plain(self):
        self.basic_run(fv_advect_axis_plain)

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

    def test_upwind_axis_run(self):
        # uses a backward euler or Matsuno scheme
        V, q = self.get_initial_conditions()
        # V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        # q = np.full((*world_shape,), standard_pressure.m, dtype=np.float) * units.Pa
        # t = np.full((*world_shape,), 273.15) * units.kelvin
        # p[half, half] += 0.0001 * units.Pa

        state = {"V": V, "q": q}

        def ft(V, q):
            q_next = upwind_axis(dt, spatial_change, V, q)
            state = {"V": V, "q": q_next}
            return state

        stable = run_1d_with_ft(state, ft)
        self.assertTrue(stable)

    def test_upwind_axis_backwards_euler(self):
        # uses a backward euler or Matsuno scheme
        V, q = self.get_initial_conditions()
        # V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        # q = np.full((*world_shape,), standard_pressure.m, dtype=np.float) * units.Pa
        # t = np.full((*world_shape,), 273.15) * units.kelvin
        # p[half, half] += 0.0001 * units.Pa

        state = {"V": V, "q": q}

        def ft(V, q):
            q_star = q - upwind_axis_finite(dt, spatial_change, V, q)
            q_next = q - upwind_axis_finite(dt, spatial_change, V, q_star)
            state = {"V": V, "q": q_next}
            return state

        stable = run_1d_with_ft(state, ft)
        self.assertTrue(stable)

    def test_pgf_matsuno(self):
        # uses a backward euler or Matsuno scheme
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        q = np.full((*world_shape,), standard_pressure.m, dtype=np.float) * units.Pa
        q[half] += 0.001 * units.Pa

        state = {"V": V, "q": q}

        def ft(V, q):
            q_star = q - fv_advect_axis_upwind_finite(dt, spatial_change, V, q)
            v_star = V + pgf_one_d(dt, spatial_change[0], q)
            q_next = q - fv_advect_axis_upwind_finite(dt, spatial_change, v_star, q_star)
            v_next = V + pgf_one_d(dt, spatial_change[0], q_star)
            # q_next = q_star
            # v_next = v_star
            state = {"V": v_next, "q": q_next}
            return state

        stable = run_1d_with_ft(state, ft)
        self.assertTrue(stable)

    def test_forward_backward(self):
        # uses the scheme from https://www.youtube.com/watch?v=1wWHqltukXo
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        H = 1 * units.m
        q = np.zeros((*world_shape,), dtype=np.float) * units.m
        q[:] = H
        q[half: half+3] += 0.1 * units.m

        state = {"V": V, "q": q}

        def ft(V, q):
            v_next = V - shallow_water_g_center_space(dt, spatial_change, q)
            q_next = q - shallow_water_h_center_space(dt, spatial_change, v_next, H)
            state = {"V": v_next, "q": q_next}
            return state

        stable = run_1d_with_ft(state, ft)
        self.assertTrue(stable)

    def test_forward_backward_advection(self):
        # uses the scheme from https://www.youtube.com/watch?v=1wWHqltukXo
        # the only issue with this scheme is that sqrt(g*H)dt/dx needs to be < 2,
        # which may be unfeasable for the grid size that I want to use.
        # minimum grid size = sqrt(g*h)dt = sqrt(10*1000)*900s = 90000m
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        H = 1000 * units.m
        q = np.zeros((*world_shape,), dtype=np.float) * units.m
        q[:] = H
        q[half: half+3] += 0.1 * units.m

        state = {"V": V, "q": q}
        spatial_change = (100000 * units.m,)
        dt = 900 * units.s

        def ft(V, q):
            v_next = V - shallow_water_g_center_space(dt, spatial_change, q)
            q_next = q - shallow_water_h_center_space(dt, spatial_change, v_next, q)
            state = {"V": v_next, "q": q_next}
            return state

        stable = run_1d_with_ft(state, ft)
        self.assertTrue(stable)

    def test_forward_backward_c_grid(self):
        # uses the scheme from https://www.youtube.com/watch?v=1wWHqltukXo
        # the only issue with this scheme is that sqrt(g*H)dt/dx needs to be < 2,
        # which may be unfeasable for the grid size that I want to use.
        # minimum grid size = sqrt(g*h)dt = sqrt(10*1000)*900s = 90000m
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        H = 1000 * units.m
        q = np.zeros((*world_shape,), dtype=np.float) * units.m
        q[:] = H
        q[half: half+3] += 1 * units.m

        state = {"V": V, "q": q}
        spatial_change = (100000 * units.m,)
        dt = 900 * units.s

        def ft(V, q):
            v_next = V - shallow_water_g_c_grid(dt, spatial_change, q)
            q_next = q - shallow_water_h_c_grid(dt, spatial_change, v_next, q)
            state = {"V": v_next, "q": q_next}
            return state

        stable = run_1d_with_ft(state, ft, steps=10000)
        self.assertTrue(stable)

