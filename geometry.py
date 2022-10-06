import math

import numpy as np

from coordinates_3d import kp
from constants import radius, units

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


def gen_geometry(height, width, layers):
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

    def manabe_sig(s):
        return s ** 2 * (3 - 2 * s)

    for i in range(layers+1):
        mysig.append(manabe_sig(1 - i/(layers)))

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
    dlat = 180 / height
    dlong = 360 / width

    sin_j = np.zeros((height,))
    cos_j = np.zeros((height,))
    sin_h = np.zeros((height,))
    cos_h = np.zeros((height,))
    for i in range(height):
        lat_j[i] = 90 - (i+0.5) * dlat
        lat_h[i] = 90 - (i+1) * dlat

    long_k = np.zeros((width,))
    for i in range(width):
        long_k[i] = -180 + (i+0.5) * dlong


    geom.lat = lat_j
    geom.long = long_k

    cos_j = np.cos(lat_j * np.pi / 180)
    sin_j = np.sin(lat_j * np.pi / 180)
    cos_h = np.cos(lat_h * np.pi / 180)
    sin_h = np.sin(lat_h * np.pi / 180)
    dx_j = cos_j * circumference / width
    dx_h = cos_h * circumference / width

    print(lat_j)
    print(lat_h)
    print(cos_j)
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

    # geom.ptop = 10 * units.hPa
    geom.ptop = 0 * units.hPa

    geom.heightmap = np.zeros((height, width)) * units.m

    return geom

