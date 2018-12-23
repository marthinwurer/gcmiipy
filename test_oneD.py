import logging
import numpy as np
import unittest

from oneD import VectorField, ScalarField, Equation, MaterialDerivative, Divergence, Gradient, Differential
from constants import ureg


def setUpModule():
    FORMAT = '%(asctime)-15s | %(filename)s:%(lineno)s | %(message)s'
    logging.basicConfig(format=FORMAT, level=logging.DEBUG)


class TestOneD(unittest.TestCase):

    def test_gradient(self):
        pass

    def test_finite_volume(self):
        """
        Goals: 1d version of primitive equations from hansen et al
        with units
        with finite volumes
        start with forward euler
        start with the basic math stuff

        for finite volume method:
            find average values in cell for conserved variables (integrate over cell volume)
            integrate over cell surfaces
                for each surface:

        divergence integrated over volume is dot with normal of surface

        so, find average value for the cell for conserved variable

        interpolate to find fluxes at the edges (or use upwinding scheme to pick which fluxes to use)

        use those to solve for next step

        for each edge:
            interpolate value of fluxes
            dot product them with normal



        """

        world_shape = (16,)

        V = VectorField(np.zeros(world_shape) * ureg.meters / ureg.second)

        p = ScalarField(np.zeros(world_shape) * ureg.pascals)

        rho = ScalarField(np.zeros(world_shape) * ureg.kg / ureg.meters ** 3)

        dt = Differential(1 * ureg.seconds)

        # conservation of momentum

        conservation_of_momentum = Equation(
            MaterialDerivative(V, dt)
            ,
            1 / rho * Gradient(p)
        )

