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

    def run_shallow(self, count, func, h, u, dx, dt):
        plt.ion()
        for i in tqdm(range(count)):
            # h, u = advect_matsumo(h, u, dt, dx)
            h, u = func(h, u, dt, dx)
            plot_callback(h)
            if np.isnan(h).any():
                # print("iteration %s" % i)
                print("NaN encountered")
                break

            c = courant_number(h, u, dx, dt)
            if c > .35:
                print("courant too high:", c)

        c = courant_number(h, u, dx, dt)
        print("final courant:", c)

        plt.ioff()
        plt.show()
        return h, u

    def test_shallow_1d(self):
        side_len = 100
        u = np.full((side_len,), 0.0) * units.m / units.s
        h = np.full((side_len,), 10.0) * units.m
        dx = 300000 * units.m
        dt = 900 * units.s

        u[-1] = 0 * u.u
        h[:50] = 20 * h.u

        c = courant_number(h, u, dx, dt)
        print("initial courant:", c)
        self.assertLess(c.m, 0.35)

        plt.ion()
        for i in tqdm(range(10000)):
            # h, u = advect_matsumo(h, u, dt, dx)
            h, u = shallow_water_matsuno(h, u, dt, dx)
            plot_callback(h)
            if np.isnan(h).any():
                # print("iteration %s" % i)
                print("NaN encountered")
                break

            c = courant_number(h, u, dx, dt)
            if c > .35:
                print("courant too high:", c)

        c = courant_number(h, u, dx, dt)
        print("final courant:", c)

        plt.ioff()
        plt.show()


    def test_advect_maccormack(self):
        side_len = 100
        u = np.full((side_len,), 10.0) * units.m / units.s
        h = np.full((side_len,), 10.0) * units.m
        dx = 300000 * units.m
        dt = 900 * units.s

        # u[-1] = 0 * u.u
        h[25:50] = 20 * h.u

        c = courant_number(h, u, dx, dt)
        print("initial courant:", c)
        h, u = self.run_shallow(1000, advect_maccormack, h, u, dx, dt)



"""
Some good tests from a real CFD guy:

You want to show that the theoretical convergence rate matches the computed rate.
So, you need an exact solution of the Navier-Stokes equations.
These are not easy to come by but there are some for different geometries.
For 1D there are a lot.
For 2D incompressible, I've used variants of the Taylor-Green vortex problem.
https://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex
When you have your exact solution, you compute the "error", which is the difference between the exact and computational result.
Then you vary your grid spacing while keeping the time step constant. You need to pick a time step that will be stable for all grids tested.
Plot the error as a function of the grid spacing.
For a second order scheme you should see that the error asymptotically approaches O(dx^2).
If you don't see that it's possible that you have a bug or need to use a smaller grid.
This is called "verification testing".
There's also "validation testing" and "uncertainty quantification".
And you can apply ideas like unit tests, integration tests, etc., but those are orthogonal to everything else I just mentioned.

This one is famous for compressible flow but I can't vouch for how useful it is: https://en.wikipedia.org/wiki/Sod_shock_tube
http://www.cfdbooks.com/
Volume 1 listed here has a ton of exact solutions that you can use for code testing.


Stuff I found:
https://www.sciencedirect.com/science/article/pii/S1674237018300255
Has initial conditions for a dam break

https://www.reading.ac.uk/web/files/maths/02-99.pdf
Has a few problems with initial conditions and analytical solutions

https://pdfs.semanticscholar.org/46ef/e430659baa9dcd3e138e6bf836dca91fa544.pdf
discusses some

https://www.cambridge.org/core/journals/acta-numerica/article/finitevolume-schemes-for-shallowwater-equations/AE2AC80D1E6E9F6BC0E68496A1C3EC52/core-reader
interesting, discusses the effects of water edges.
q is discharge = h * v (defined after 2.4)
talks about dealing with shallow water edges. draining time steps looks cool.

https://sci-hub.tw/https://ascelibrary.org/doi/full/10.1061/%28ASCE%29HY.1943-7900.0001683
analytical solution for dam break




"""




if __name__ == '__main__':
    unittest.main()
