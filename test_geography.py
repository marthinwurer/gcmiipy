import matplotlib.pyplot as plt

from no_limits_2_5d import *

def plot_callback(p, u, v, t, q):
    quantity = p
    plt.clf()
    plt.imshow(quantity)
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
    p, u, v, t, q, g, geom = run_model(8, 8, 3, 60 * 30 * units.s, 1200 * 12, None)
    plt.ioff()
    print("ground temp:", g.gt)
    tp = p * geom.sig + geom.ptop
    tt = temperature.to_true_temp(t, tp)
    print("atmosphere temps:", tt)
    print("pressures:", tp)


if __name__ == "__main__":
    main()
