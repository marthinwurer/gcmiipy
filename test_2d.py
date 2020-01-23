import logging
import unittest

import numpy as np
import matplotlib.pyplot as plt

from constants import units, standard_pressure
from two_d import *

logger = logging.getLogger(__name__)

def setUpModule():
    FORMAT = '%(asctime)-15s | %(filename)s:%(lineno)s | %(message)s'
    logging.basicConfig(format=FORMAT, level=logging.INFO)

dx = 10 * units.m
dy = 10 * units.m
dt = 1 * units.s
# side_length = 160
# side_length = 16
side_length = 4
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
    V[0][:] = 2.0 * units.m / units.s
    V[1][:] = -2.0 * units.m / units.s

    return V, q


class Test2d(unittest.TestCase):
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

    def test_advect_v(self):
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        p = np.full((*world_shape,), standard_pressure.m, dtype=np.float) * units.Pa
        t = np.full((*world_shape,), 273.15) * units.kelvin

        # initial conditions
        # V[0][quarter:half, quarter:half] = 1.0 * units.m / units.s
        # V[1][half] = -1.0 * units.m / units.s

        p[half, half] += 0.0001 * units.Pa

        plt.ion()
        plt.figure()
        plt.imshow(V[0])
        plt.title('Initial state')
        plt.show()
        p_prev = p
        v_prev = V

        initial_variation = get_total_variation(p)
        logger.info("Initial Variation: %s" % (initial_variation,))

        for i in range(num_steps):

            p_next = finite_volume_advection(dt, spatial_change, v_prev, p_prev)
            # v_x_next = fv_advect_axis_upwind(dt, spatial_change, v_prev, v_prev[0])
            # v_y_next = fv_advect_axis_upwind(dt, spatial_change, v_prev, v_prev[1])
            pgf = pgf_c_grid(dt, spatial_change, p_prev, t)
            # v_next = unit_stack([v_x_next, v_y_next]) - pgf
            v_next = v_prev - pgf
            # t_next = finite_volume_advection(dt, spatial_change, v_prev, t)
            # v_next = -pgf
            p_prev = p_next
            v_prev = v_next
            # t = t_next
            current_variation = get_total_variation(p_prev)
            if current_variation.m > initial_variation.m + 0.00001:
                print("Variation too high: %s" % (current_variation,))
                # self.assertEqual(True, False)

            print("iteration %s" % i)
            plt.clf()
            plt.imshow(pgf[0])
            # plt.imshow(p_next)
            plt.title('current_state...')
            plt.show()
            plt.pause(0.001)  # pause a bit so that plots are updated

        else:
            p_next = p

        final_variation = get_total_variation(p_next)
        logger.info("Initial Variation: %s Final Variation: %s" % (initial_variation, final_variation))

        plt.ioff()
        plt.show()

    def test_fv(self):
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

            q_next = finite_volume_advection(dt, spatial_change, V, q_prev)
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

    def test_pressure_at_edges(self):
        p = np.full((3, 3), 1, dtype=np.float) * units.Pa
        p[1, 1] = 0 * units.Pa
        p_edge = pressure_at_edge(p)
        self.assertEqual(p.u, p_edge.u)
        self.assertEqual(0.5 * units.Pa, p_edge[0][0][1])


    def test_advection_of_momentum(self):
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        p = np.full((*world_shape,), standard_pressure.m, dtype=np.float) * units.Pa
        # t = np.full((*world_shape,), 273.15) * units.kelvin

        p_next = advect_with_momentum(dt, spatial_change, V, p)
        p_next = p_next.to_base_units()

        self.assertEqual(p.u, p_next.u)

    def test_pgf_one_d(self):
        U = np.zeros((1, 3)) * units.m / units.s
        p = np.full((3,), 1, dtype=np.float) * units.Pa
        p[1, 1] = 0.1 * units.Pa
        p_edge = pressure_at_edge_one_d(p)

        V = U * p_edge

        pgf = pgf_one_d(dt, dx, p).to_base_units()

        v_next = V - pgf * p_edge

        # We just want this to work without exceptions due to units
        self.assertIsNotNone(v_next)

    def test_pgf_units(self):
        U = np.zeros((2, 3, 3)) * units.m / units.s
        p = np.full((3, 3), 1, dtype=np.float) * units.Pa
        p[1, 1] = 0.1 * units.Pa
        p_edge = pressure_at_edge(p)

        V = U * p_edge

        pgf = pgf_templess(dt, spatial_change, p)

        v_next = V - pgf * p_edge

        # We just want this to work without exceptions due to units
        self.assertIsNotNone(v_next)

    def test_pgf_temp_units(self):
        U = np.zeros((2, 3, 3)) * units.m / units.s
        p = np.full((3, 3), 1, dtype=np.float) * units.Pa
        t = np.full((3, 3), standard_temperature.m) * units.K
        p[1, 1] = 0.1 * units.Pa
        p_edge = pressure_at_edge(p)

        V = U * p_edge

        pgf = pgf_c_grid(dt, spatial_change, p, t)

        v_next = V - pgf * p_edge

        # We just want this to work without exceptions due to units
        self.assertIsNotNone(v_next)

    def test_run_func(self):
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        q = np.zeros((*world_shape,), dtype=np.float) * units.grams
        q[quarter:half, quarter:half] = 1.0 * units.grams
        V[0][:] = 2.0 * units.m / units.s
        V[1][:] = -2.0 * units.m / units.s

        state = {"V": V, "q": q}

        def ft(V, q):
            q_next = corner_transport_2d(dt, spatial_change, V, q)
            state = {"V": V, "q": q_next}
            return state

        stable = run_2d_with_ft(state, ft)
        self.assertTrue(stable)

    def test_iterative_pgf(self):
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        p = np.full((*world_shape,), standard_pressure.m, dtype=np.float) * units.Pa
        t = np.full((*world_shape,), 273.15) * units.kelvin
        p[half, half] += 0.0001 * units.Pa

        state = {"V": V, "p": p, "t": t}

        def ft(V, p, t):
            p_next = finite_volume_advection(dt, spatial_change, V, p)
            pgf = pgf_c_grid(dt, spatial_change, p, t)
            v_next = V - pgf
            state = {"V": v_next, "p": p_next, "t": t}
            return state

        stable = run_2d_with_ft(state, ft, display_key="p", variation_key="p")
        self.assertTrue(stable)

    def test_backward_pgf(self):
        # uses a backward euler or Matsuno scheme
        V = np.zeros((len(world_shape), *world_shape)) * units.m / units.s
        p = np.full((*world_shape,), standard_pressure.m, dtype=np.float) * units.Pa
        t = np.full((*world_shape,), 273.15) * units.kelvin
        p[half, half] += 0.0001 * units.Pa

        state = {"V": V, "p": p, "t": t}

        def ft(V, p, t):
            p_star = finite_volume_advection(dt, spatial_change, V, p)
            pgf_star = pgf_c_grid(dt, spatial_change, p, t)
            v_star = V - pgf_star
            p_next = finite_volume_advection(dt, spatial_change, v_star, p_star)
            pgf_next = pgf_c_grid(dt, spatial_change, p_star, t)
            v_next = V - pgf_next
            state = {"V": v_next, "p": p_next, "t": t}
            return state

        stable = run_2d_with_ft(state, ft, display_key="p", variation_key="p")
        self.assertTrue(stable)




if __name__ == '__main__':
    unittest.main()
