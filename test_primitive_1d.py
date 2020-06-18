import unittest

import matplotlib.pyplot as plt
from tqdm import tqdm

from constants import *
from matsuno_c_grid import courant_number
from primitive_1d import *

side_len = 8

def plot_callback(q):
    quantity = q
    plt.clf()
    plt.plot(quantity)
    # plt.title('n = %s' % (i,))
    # ax = plt.gca()
    # ax.format_coord = lambda x, y: f'{int(x + .5)} {int(y + .5)} {quantity[int(y + .5), int(x + .5)]}'
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated


class TestPrimOneD(unittest.TestCase):
    def test_advect_rho(self):
        # stable, but not TVD
        rho = np.full((side_len,), rd.m) * rd.u
        u = np.full((side_len,), 10) * units.m / units.s
        dx = 100 * units.m
        dt = 1 * units.s

        rho[25:50] *= 2

        plt.ion()
        for i in range(1000):
            rho_next, u_next = advect_matsumo(rho, u, dx, dt)

            rho = rho_next
            plot_callback(rho)

        plt.ioff()
        plt.show()

    def test_advect_rho_fe(self):
        # stable, but not TVD
        rho = np.full((side_len,), rd.m) * rd.u
        height = np.full((side_len,), 1000.0) * units.m
        u = np.full((side_len,), 10) * units.m / units.s
        dx = 100000 * units.m
        dt = 300 * units.s

        height[25:50] *= 1.5

        plt.ion()
        initial_variation = get_total_variation(rho)

        for i in tqdm(range(100000)):
            height_next, u_next = advect_matsumo(height, u, dx, dt)

            height = height_next
            plot_callback(height)
            print(courant_number(height, u, dx, dt))
            current_variation = get_total_variation(rho)
            if current_variation.m > initial_variation.m + 0.1:
            #     print("iteration %s" % i)
                print("Variation too high: %s" % (current_variation,))
            if np.isnan(rho).any():
                # print("iteration %s" % i)
                print("NaN encountered")
                break

        plt.ioff()
        plt.show()

    def test_aflux_units(self):
        # just test to make sure that units are good.
        u = np.full((side_len,), 10) * units.m / units.s
        p = np.full((side_len,), standard_pressure.m) * standard_pressure.u
        dx = 100000 * units.m
        dy = 100000 * units.m
        dt = 300 * units.s

        pu, conv, pit = aflux(u, p, dy)

        pa = advecm(p, pit, dt, dx * dy)

        self.assertEqual(p.u, pa.u)

    def test_advecv_units(self):
        # just test to make sure that units are good.
        u = np.full((side_len,), 10) * units.m / units.s
        p = np.full((side_len,), standard_pressure.m) * standard_pressure.u
        dx = 100000 * units.m
        dy = 100000 * units.m
        dt = 300 * units.s

        pu, conv, pit = aflux(u, p, dy)

        pa = advecm(p, pit, dt, dx * dy)

        u_next = advecv(u, pu, p, pa, dt, dx)

        self.assertEqual(u.u, u_next.u)

    def test_pgf_units(self):
        # just test to make sure that units are good.
        # side_len = 5
        u = np.full((side_len,), 10) * units.m / units.s
        p = np.full((side_len,), standard_pressure.m) * standard_pressure.u
        t = np.full((side_len,), standard_temperature.m) * standard_temperature.u
        dx = 100000 * units.m
        dy = 100000 * units.m
        dt = 300 * units.s

        p[2] *= 1.1

        pu, conv, pit = aflux(u, p, dy)

        pa = advecm(p, pit, dt, dx * dy)

        u_next = advecv(u, pu, p, pa, dt, dx)

        spa, theta, phi, geo, pg, u_next = pgf(u, u_next, p, pa, t, dt, dx)
        print("spa", spa[0].to_base_units())
        print("theta", theta[0].to_base_units())
        print("phi", phi[0].to_base_units())

        self.assertEqual(u.u, u_next.u)

    def test_advect_units(self):
        u = np.full((side_len,), 10) * units.m / units.s
        p = np.full((side_len,), standard_pressure.m) * standard_pressure.u
        t = np.full((side_len,), standard_temperature.m) * standard_temperature.u
        dx = 100000 * units.m
        dy = 100000 * units.m
        dt = 300 * units.s

        p[2] *= 1.1

        pu, conv, pit = aflux(u, p, dy)

        pa = advecm(p, pit, dt, dx * dy)

        t_next = advect(pu, p, t, pa, t, dt, dx)

        self.assertEqual(t.u, t_next.u)

    def test_advecq_units(self):

        u = np.full((side_len,), 10) * units.m / units.s
        p = np.full((side_len,), standard_pressure.m) * standard_pressure.u
        t = np.full((side_len,), standard_temperature.m) * standard_temperature.u
        q = np.full((side_len,), 1) * units.gram / units.kg  # specific humidity is mostly dimensionless
        dx = 100000 * units.m
        dy = 100000 * units.m
        dt = 300 * units.s

        q[2] *= 1.1

        pu, conv, pit = aflux(u, p, dy)

        pa = advecm(p, pit, dt, dx * dy)

        q_next = advecq(pu, p, q, pa, q, dt, dx)

        self.assertEqual(q.u, q_next.u)

    def test_advecq(self):
        # initial conditions
        side_len = 100
        u = np.full((side_len,), 10) * units.m / units.s
        p = np.full((side_len,), standard_pressure.m) * standard_pressure.u
        t = np.full((side_len,), standard_temperature.m) * standard_temperature.u
        q = np.full((side_len,), 0) * units.gram / units.kg  # specific humidity is mostly dimensionless
        dx = 100000 * units.m
        dy = 100000 * units.m
        dt = 300 * units.s

        q[25:50] = 1.1 * units.gram / units.kg
        plt.ion()
        for i in range(10000):
            pu, conv, pit = aflux(u, p, dy)

            pa = advecm(p, pit, dt, dx * dy)

            q_star = advecq(pu, p, q, pa, q, dt, dx)

            q_next = advecq(pu, p, q, pa, q_star, dt, dx)

            q = q_next
            plot_callback(q)

        plt.ioff()
        plt.show()





if __name__ == '__main__':
    unittest.main()
