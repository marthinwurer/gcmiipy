import matplotlib.pyplot as plt

from no_limits_2_5d import *


def run_model(height, width, layers, dt, timesteps, callback):
    geom = gen_geometry(height, width, layers, sig_func=geometry.manabe_sig)
    p, u, v, t, q, g = gen_initial_conditions(geom)
    utc = 0 * units.hours
    # u *= 0
    v[0,0,0] = 0.1 * v.u
    u *= 0
    geom.heightmap.m[0, 8] = 1000

    # p[0, 0] *= 1.01

    for i in tqdm(range(timesteps)):
        p, u, v, t, q, g = full_timestep(p, u, v, t, q, g, dt, utc, geom)
        utc += dt
        if callback:
            callback(p, u, v, t, q)

    return p, u, v, t, q, g, geom


def plot_callback(p, u, v, t, q):
    quantity = u[0]
    if len(STATS["ke"]) % 100 != 99:
        return
    plt.clf()
    # plt.imshow(quantity)
    plt.plot([i[3].m for i in STATS["ke"]])
    # plt.title('n = %s' % (i,))
    # ax = plt.gca()
    # ax.format_coord = lambda x, y: f'{int(x + .5)} {int(y + .5)} {quantity[int(y + .5), int(x + .5)]}'
    plt.show()
    plt.pause(0.001)  # pause a bit so that plots are updated


def main():
    # test_humidity_calcs()
    # exit()


    # run_model(height, width, layers, 60 * 15 * units.s, 1000, None)
    # p, u, v, t, q, g, geom = run_model(1, 1, 18, 60 * 15 * units.s, 3, None)
    plt.ion()
    try:
        p, u, v, t, q, g, geom = run_model(1, 16, 17, 60 * 30 * units.s, 1200 * 12, plot_callback)
    except:
        print("exception")
        raise
    plt.ioff()
    print("ground temp:", g.gt)
    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    print("atmosphere temps:", tt)
    print("pressures:", tp)


if __name__ == "__main__":
    main()
