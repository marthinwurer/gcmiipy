import math

import numpy as np

from coordinates_3d import kp
import constants
from constants import radius, units, unit_roll, G, Md, R, Rd

class Geom:
    def __init__(self, height, width, layers):
        self.height = height
        self.width = width
        self.layers = layers
        self.sige = []
        self.dsig = []
        self.sigb = []
        self.sigt = []
        self.sig = []
        # self.dx = 0
        self.dy = 0
        self.lat = []
        self.long = []
        self.dx_j = 0
        self.dx_h = 0
        self.ptop = 0
        self.heightmap = []


# functions for determining sigma 
def manabe_sig(s):
    return s ** 2 * (3 - 2 * s)


def equal_sig(s):
    return s


def gen_geometry(height, width, layers, sig_func=equal_sig,
                 north_edge=90, south_edge=-90, west_edge=-180, east_edge=180):
    """"""
    """
      DATA SIGE/    1.,.948665,.866530,.728953,.554415,.390144,
     *  .251540,.143737,.061602,28*0./                        
    """
    sige = np.asarray([1., .948665, .866530, .728953, .554415, .390144, .251540, .143737, .061602, 0.])

    """
C**** CALCULATE DSIG AND DSIGO                                           816.   
      DO 700 L=1,LM                                                      817.   
  700 DSIG(L)=SIGE(L)-SIGE(L+1)                                          818.   
      DO 710 L=1,LM-1                                                    819.   
  710 DSIGO(L)=SIG(L)-SIG(L+1)                                           820.   
    """

    # DATA PLB4 in R83ZAmacDBL.f seems to have the base layer values for pressure
    # I assume that sigma is calculated off of that?

    plb4_4 = [
        1013.2500, 1000.0000, 950.0000, 900.0000, 850.0000, 800.0000,
        750.0000, 700.0000, 650.0000, 600.0000, 550.0000, 500.0000,
        450.0000, 400.0000, 350.0000, 300.0000, 250.0000, 200.0000,
        150.0000, 100.0000, 50.0000, 20.0000, 10.0000, 5.0000,
        2.0000, 1.0000, 0.5000, 0.2000, 0.1000, 0.0500,
        0.0200, 0.0100, 0.0050, 0.0020, 0.0010, 1.E-05,
        0.
    ]



    geom = Geom(height, width, layers)
    mysig = []

    for i in range(layers+1):
        mysig.append(sig_func(1 - i/(layers)))

    def rs(arr):
        return np.reshape(arr, (arr.shape[0], 1, 1))

    geom.sige = rs(np.asarray(mysig))
    geom.sigt = rs(np.asarray(mysig[1:]))
    geom.sigb = rs(np.asarray(mysig[:-1]))

    geom.dsig = geom.sigb - geom.sigt
    geom.sig = (geom.sigb + geom.sigt) / 2
    geom.dsigv = kp(geom.sig) - geom.sig
    # TODO I might need another dsig for inbetween layers for vertical advection

    circumference = 2 * radius * math.pi
    lat_j = np.zeros((height,))
    lat_h = np.zeros((height,))
    dlat = (north_edge - south_edge) / height
    dlong = (east_edge - west_edge) / width

    sin_j = np.zeros((height,))
    cos_j = np.zeros((height,))
    sin_h = np.zeros((height,))
    cos_h = np.zeros((height,))
    for i in range(height):
        lat_j[i] = north_edge - (i+0.5) * dlat
        lat_h[i] = north_edge - (i+1) * dlat

    long_k = np.zeros((width,))
    for i in range(width):
        long_k[i] = west_edge + (i+0.5) * dlong


    geom.lat = lat_j.reshape((height, -1)) * units.degrees
    geom.long = long_k * units.degrees

    cos_j = np.cos(lat_j * np.pi / 180)
    sin_j = np.sin(lat_j * np.pi / 180)
    cos_h = np.cos(lat_h * np.pi / 180)
    sin_h = np.sin(lat_h * np.pi / 180)
    dx_j = cos_j * circumference / width
    dx_h = cos_h * circumference / width

    print(lat_j)
    print(lat_h)
    print(cos_j)
    print(cos_h)
    print(dx_j)
    print(dx_h)






    # plt.plot(plb4_4)
    # plt.plot(sige)
    # plt.plot(mysig)
    # plt.show()

    # TODO spherical geometry
    # geom.dx = circumference / width
    geom.dx_j = np.reshape(dx_j, (1, height, 1))
    geom.dx_h = np.reshape(dx_h, (1, height, 1))
    geom.dy = circumference / 2 / height

    # aproximate area with trapezoids
    top = (unit_roll(dx_h, 1, axis=0) + dx_h) * geom.dy * 0.5
    geom.area = top
    print(top)
    # exit()

    # geom.ptop = 10 * units.hPa
    geom.ptop = 0 * units.hPa

    geom.heightmap = np.zeros((height, width)) * units.m

    return geom


def gen_square_geometry(height, width, layers, dx, dy, sig_func=equal_sig):
    geom = Geom(height, width, layers)
    geom.ptop = 0 * units.hPa

    mysig = []

    for i in range(layers+1):
        mysig.append(sig_func(1 - i/(layers)))

    def rs(arr):
        return np.reshape(arr, (arr.shape[0], 1, 1))

    geom.sige = rs(np.asarray(mysig))
    geom.sigt = rs(np.asarray(mysig[1:]))
    geom.sigb = rs(np.asarray(mysig[:-1]))

    geom.dsig = geom.sigb - geom.sigt
    geom.sig = (geom.sigb + geom.sigt) / 2
    geom.dsigv = kp(geom.sig) - geom.sig

    geom.lat = 0 * units.degrees
    geom.long = 0 * units.degrees

    geom.dx_j = np.full((1, height, 1), dx.m) * units.m
    geom.dx_h = np.full((1, height, 1), dx.m) * units.m
    geom.dy = dy

    geom.heightmap = np.zeros((height, width)) * units.m
    return geom


def pressure_from_heightmap(height, sea_level_pressure, sea_level_temp):
    """
    https://en.wikipedia.org/wiki/Barometric_formula
    """
    top = -G * Md * height
    bottom = R * sea_level_temp
    div = top / bottom
    print(top.to(units.J / units.mol))
    print(bottom)
    print(div.to_base_units())
    print(Rd * sea_level_temp)

    # second attempt:
    # assuming constant temperature
    # left = np.sqrt(height * G * Rd * sea_level_temp + sea_level_pressure ** 2)
    # third attempt
    # suddenly quadratic
    # a = 1
    # b = sea_level_pressure
    # c = height * G / (Rd * sea_level_temp)
    # print(a)
    # print(b.to_base_units())
    # print(c.to_base_units())
    # val = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

    # fourth attempt
    # p1 = sea_level_pressure * Rd * sea_level_temp / (height * G)

    # fifth attempt, ft sympy
    p1 = -Rd * sea_level_temp * sea_level_pressure / (height * G - Rd * sea_level_temp)
    m = Rd * sea_level_temp
    print("m", m)

    # sixth attempt, with more sympy:
    # [p2*(2*Rd*tt + g*h1 - g*h2)/(2*Rd*tt - g*h1 + g*h2)]
    p2 = sea_level_pressure
    tt = sea_level_temp
    h1 = height
    h2 = 0 * units.m
    g = G
    p1 = p2*(2*Rd*tt + g*h1 - g*h2)/(2*Rd*tt - g*h1 + g*h2)


    wiki_val = sea_level_pressure * np.exp((-G * Md * height) / (R * sea_level_temp))
    print(wiki_val)
    print("p1", p1.to(units.Pa))
    # exit()

    return wiki_val
    # return p1.to(units.Pa)

