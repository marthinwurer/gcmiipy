import unittest

import matplotlib.pyplot as plt
import numpy as np

from coordinates_1d import *
from constants import units


def van_leer(r):
    return (r + np.abs(r)) / (1 + np.abs(r))


def calc_r(q):
    units = q.u
    q = q.m
    a = q - im(q)
    b = ip(q) - q
    c = np.divide(a, b, out=np.zeros_like(a), where=(b != 0))
    return c * units


def donor_cell_flux(q, u):
    q_edge = np.where(u > 0, q, ip(q)) * q.u
    flux = q_edge * u

    return flux

def donor_cell_advection(q, u, dx, dt):
    flux = donor_cell_flux(q, u)
    return q + (im(flux) - flux) * dt / dx


def plot_callback(q):
    quantity = q
    plt.clf()
    plt.plot(quantity)
    # plt.title('n = %s' % (i,))
    # ax = plt.gca()
    # ax.format_coord = lambda x, y: f'{int(x + .5)} {int(y + .5)} {quantity[int(y + .5), int(x + .5)]}'
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated

class TestFluxLimiter(unittest.TestCase):
    def test_van_leer(self):
        self.assertEqual(1, van_leer(1))
        self.assertEqual(0, van_leer(0))

    def test_calc_r(self):
        side_len = 16
        u = np.full((side_len,), 10) * units.m / units.s
        q = np.full((side_len,), 0) * units.gram / units.kg  # specific humidity is mostly dimensionless
        q[4:8] = 1.1 * units.gram / units.kg

        r = calc_r(q)

        plt.plot(r)
        plt.show()

    def test_van_leer_with_r(self):
        side_len = 16
        u = np.full((side_len,), 10.0) * units.m / units.s
        q = np.full((side_len,), 0.0) * units.gram / units.kg  # specific humidity is mostly dimensionless
        q[4:8] = 1 * units.gram / units.kg

        r = calc_r(q)

        print(r)
        plt.plot(r)
        plt.show()

    def test_donor_cell_flux_positive(self):
        side_len = 16
        u = np.full((side_len,), 10.0) * units.m / units.s
        q = np.full((side_len,), 0.0) * units.gram / units.kg  # specific humidity is mostly dimensionless
        q[4:8] = 1 * q.u
        dx = 100 * units.m
        dt = 1 * units.s

        plt.ion()
        for i in range(100):
            q_next = donor_cell_advection(q, u, dx, dt)

            q = q_next
            print(q)
            plot_callback(q)

        plt.ioff()

    def test_donor_cell_flux_negative(self):
        side_len = 160
        u = np.full((side_len,), -10.0) * units.m / units.s
        q = np.full((side_len,), 0.0) * units.gram / units.kg  # specific humidity is mostly dimensionless
        q[40:80] = 1 * units.gram / units.kg
        dx = 100 * units.m
        dt = 1 * units.s

        plt.ion()
        for i in range(100):
            q_next = donor_cell_advection(q, u, dx, dt)

            q = q_next
            print(q)
            plot_callback(q)

        plt.ioff()




if __name__ == '__main__':
    unittest.main()
