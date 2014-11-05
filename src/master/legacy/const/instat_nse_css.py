# -*- coding: utf-8 -*-
from dolfin import __version__
from dolfin import grad,Identity


DOLFIN_VERSION      = __version__


"""output for velocity and pressure """
def INSTAT_NSE_CSS_VELOCITY_PVD_FILE(RE,refalg,reflevel):
    return "instat_nse_css/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/velo.pvd"%(refalg,reflevel,RE)

def INSTAT_NSE_CSS_PRESSURE_PVD_FILE(RE,refalg,reflevel):
    return "instat_nse_css/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/pressure.pvd"%(refalg,reflevel,RE)

def INSTAT_NSE_CSS_VELOCITY_TIMESERIES_FILE(RE,refalg,reflevel):
    return "instat_nse_css/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/velo.hdf5"%(refalg,reflevel,RE)

def INSTAT_NSE_CSS_PRESSURE_TIMESERIES_FILE(RE,refalg,reflevel):
    return "instat_nse_css/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/pressure.hdf5"%(refalg,reflevel,RE)


def epsilon(u):
    "Return symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, nu):
    "Return stress tensor."
    return 2*nu*epsilon(u) - p*Identity(u.cell().d)
