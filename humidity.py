from constants import *


def manabe_rh(geom):
    # relative humidity distribution from Manabe 1967
    rh = 0.77 * (geom.sig - 0.02) / (1 - 0.02)
    return rh


def saturation_vapor_pressure(tt):
    # use the buck equation to estimate the svp
    # https://en.wikipedia.org/wiki/Vapour_pressure_of_water#Accuracy_of_different_formulations
    t = tt.to(units.celsius).m
    return 0.61121 * units.kPa * np.exp((18.678 - t / 234.5) * (t / (257.14 + t)))


def w_s_at(tp, tt):
    e_s = saturation_vapor_pressure(tt)
    w_s = (Rd / Rv) * e_s / (tp - e_s)
    return w_s


def vmr_from_mmr(mmr, mmg, mma):
    return mma / mmg * mmr


def rh_to_mmr(rh, tp, tt):
    # using formulats from https://earthscience.stackexchange.com/a/5077
    e_s = saturation_vapor_pressure(tt)

    e = rh * e_s
    # This doesn't work because I'm not using Clausius-Clapeyron, I'm using buck.
    w = e * Rd / (Rv * (tp - e))
    w_s = w_s_at(tp, tt)
    # w = rh * w_s
    q = w / (w + 1)
    return q


def mmr_to_rh(mmr, tp, tt):
    e_s = saturation_vapor_pressure(tt)
    # https://www.e-education.psu.edu/meteo300/node/519
    # w is the water vapor mixing ratio: density of water over density of dry air without water
    w = mmr / (1 - mmr)
    # w = e * Rd / (Rv * (tp - e))
    """
    w = e * Rd / (Rv * (tp - e))
    w = e * eps / (p - e)
    w * (p - e) = e * eps
    w * p - w * e = e * eps
    w * p = e * eps + w * e
    w * p = e * (eps + w)
    w * p / (eps + w) = e
    """ 
    e = w * tp / (Rd / Rv + w)
    # e = w * (p - 
    # w_s = w_s_at(tp, tt)
    # rh = w / w_s
    rh = e / e_s
    return rh


def test_humidity_calcs():
    # rh = 0.5
    # t = units.Quantity(50, units.celsius).to_base_units()
    # p = standard_pressure
    # rho = (p / (constants.Rd * t)).to_base_units()
    # mmr = rh_to_mmr(rh, p, t).to_base_units().m
    # vmr = vmr_from_mmr(mmr, M_water, Md).to_base_units().m
    # print(mmr, vmr, rho, rho * vmr, rho * mmr)
    # exit()

    for i in tqdm(range(101)):
        t = units.Quantity(i, units.celsius).to_base_units()
        for j in range(100):
            p = (j + 1) * 10 * units.hPa
            for k in range(10):
                rh = (k + 1) / 10
                mmr = rh_to_mmr(rh, p, t)
                rh_back = mmr_to_rh(mmr, p, t)
                within = np.abs(rh - rh_back) < 1e-6
                if not within:
                    print(rh, rh_back, mmr, t, p, saturation_vapor_pressure(t), w_s_at(p, t).to_base_units(), rh/w_s_at(p, t))
                assert within

