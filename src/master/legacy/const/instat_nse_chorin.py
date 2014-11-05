# -*- coding: utf-8 -*-
from dolfin import __version__

DOLFIN_VERSION      = __version__


"""output for velocity and pressure """
def INSTAT_NSE_CHORIN_VELOCITY_PVD_FILE(RE,refalg,reflevel):
    return "instat_nse_chorin/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/velo.pvd"%(refalg,reflevel,RE)

def INSTAT_NSE_CHORIN_PRESSURE_PVD_FILE(RE,refalg,reflevel):
    return "instat_nse_chorin/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/pressure.pvd"%(refalg,reflevel,RE)

def INSTAT_NSE_CHORIN_VELOCITY_TIMESERIES_FILE(RE,refalg,reflevel):
    return "instat_nse_chorin/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/velo.hdf5"%(refalg,reflevel,RE)

def INSTAT_NSE_CHORIN_PRESSURE_TIMESERIES_FILE(RE,refalg,reflevel):
    return "instat_nse_chorin/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/pressure.hdf5"%(refalg,reflevel,RE)


