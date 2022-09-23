import numpy as np

from constants import *
import temperature

pressures = np.asarray([0.3734, 3.9564, 66.939, 110.91, 868.02, 5474.9, 22632, 108900]) * units.Pa
temperatures = np.asarray([-86.28, -58.5, -2.5, -2.5, -44.5, -56.5, -56.5, 19.0]) * units.celsius


def temp_at(p):
    return np.interp(p, pressures, temperatures).to_base_units()



if __name__ == "__main__":
    import no_limits_2_5d
    height = 1
    width = 1
    layers = 18
    geom = no_limits_2_5d.gen_geometry(height, width, layers)
    p = np.full((height, width), 1) * standard_pressure - geom.ptop
    u = np.full((layers, height, width), 1) * 1.0 * units.m / units.s
    v = np.full((layers, height, width), 1) * .0 * units.m / units.s
    q = np.full((layers, height, width), 1) * 0.1 * units.dimensionless

    tp = p * geom.sig + geom.ptop
    t = temperature.to_potential_temp(temp_at(tp), p)
    t = temperature.to_potential_temp(temp_at(tp), tp)
    print("Pressure:")
    print(tp)
    print("temperature:")
    print(t)

    # t = np.full((layers, height, width), 1) * temperature.to_potential_temp(standard_temperature, p)


    geopotential = no_limits_2_5d.compute_geopotential(p, t, geom).to_base_units()

    print("geopotential:")
    print(geopotential / G)








