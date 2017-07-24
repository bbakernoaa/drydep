__author__ = 'barry'


# constants to be used throughout the file
from numpy import pi

kb = 1.38e-23
mua = 1.89E-5
nua = 1.57E-5
g = 9.81
kappa = 0.4
cpa = 1000.
rhoa = 1.2

amfp = 0.067e-6


def calc_monin_obukhov(ustar, t, h):
    # monin obhukov length
    L = -ustar ** 3 / kappa * t / g * rhoa * cpa / h
    return L


def cfac(d):
    # cunningham correction factor #taken from Zhang 2001
    from numpy import exp

    return 1 + (2 * amfp / d) * (1.257 + 0.4 * exp(-0.55 * d / amfp))


def brownian_diffusion(d, t):
    cu = cfac(d)
    db = (cu * kb * t) / (3. * pi * mua * d)
    return db


def schmidt_number(d, t):
    """

    @rtype : object
    """
    db = brownian_diffusion(d, t)
    sc = nua / db
    return sc


def vg_calc(rho, d, cc):
    return (rho * 9.81 * (d ** 2) * cc) / (18 * mua)


def kinematic_viscosity(temp):
    return 145.8 * 1.e-8 * (temp ** 1.5) / (temp + 110.4)


def calc_reynolds(z0, ustar, anu):
    return z0 * ustar / anu


def dynamic_viscosity(amu):
    return amu / rhoa


def relaxation_time(d, rhop):
    cu = cfac(d)
    trel = rhop * (d ** 2.) * cu / (18. * mua)
    return trel


def ws(d, rhop):
    trel = relaxation_time(d, rhop)
    ws = g * trel
    return ws