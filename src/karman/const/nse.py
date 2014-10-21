# -*- coding: utf-8 -*-
from dolfin import __version__

DOLFIN_VERSION      = __version__


"""output for velocity and pressure """
def NSE_VELOCITY_PVD_FILE(RE,refalg,reflevel):
    return "nse/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/velo.pvd"%(refalg,reflevel,RE)

def NSE_PRESSURE_PVD_FILE(RE,refalg,reflevel):
    return "nse/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/pressure.pvd"%(refalg,reflevel,RE)



