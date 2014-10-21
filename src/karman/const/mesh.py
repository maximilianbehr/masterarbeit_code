# -*- coding: utf-8 -*-
from dolfin import __version__
from dolfin.mesh import refinement


"""constants for rectangular domain"""
rect         = {"x0":0,"x1":4,"y0":0,"y1":1}

"""constants for circualar domain"""
circ         = {"x0":1.5,"y0":0.5,"r":0.15,"fragments":15}

"""the resolution of the macro karman mesh"""
initial_resolution                      = 5

"""Indices for Boundary parts"""
GAMMALEFT                               = 1
GAMMALOWER                              = 2
GAMMARIGHT                              = 3
GAMMAUPPER                              = 4
GAMMABALL                               = 5
GAMMABALLCTRLUPPER                      = 6
GAMMABALLCTRLLOWER                      = 7
GAMMAINNER                              = 0

"""dolfin version"""
DOLFIN_VERSION                          = __version__

"""macro karman mesh files and boundary functions"""
KARMAN_MACRO_MESH_FILE                  = "karman/"+DOLFIN_VERSION+"/macro/macro.xml.gz"
KARMAN_MACRO_MESH_FILE_XDMF             = "karman/"+DOLFIN_VERSION+"/macro/macro.xdmf"
KARMAN_MACRO_BOUNDARY_FILE              = "karman/"+DOLFIN_VERSION+"/macro/boundary.xml.gz"
KARMAN_MACRO_BOUNDARY_FILE_XDMF         = "karman/"+DOLFIN_VERSION+"/macro/boundary.xdmf"
KARMAN_MACRO_MESH_PVD_FILE              = "karman/"+DOLFIN_VERSION+"/macro/macro.pvd"
KARMAN_MACRO_BOUNDARY_PVD_FILE          = "karman/"+DOLFIN_VERSION+"/macro/boundary.pvd"

"""output for refined mesh files"""
def KARMAN_REFN_MESH_FILE(refinementalg, num):
    return "karman/"+DOLFIN_VERSION+"/%s/ref_%d/mesh.xml.gz"%(refinementalg,num)

def KARMAN_REFN_BOUNDARY_FILE(refinementalg, num):
    return "karman/"+DOLFIN_VERSION+"/%s/ref_%d/boundary.xml.gz"%(refinementalg,num)

def KARMAN_REFN_MESH_FILE_XDMF(refinementalg, num):
    return "karman/"+DOLFIN_VERSION+"/%s/ref_%d/mesh.xdmf"%(refinementalg,num)

def KARMAN_REFN_BOUNDARY_FILE_XDMF(refinementalg, num):
    return "karman/"+DOLFIN_VERSION+"/%s/ref_%d/boundary.xdmf"%(refinementalg,num)

def KARMAN_REFN_MESH_PVD_FILE(refinementalg,num):
    return  "karman/"+DOLFIN_VERSION+"/%s/ref_%d/mesh.pvd"%(refinementalg,num)

def KARMAN_REFN_BOUNDARY_PVD_FILE(refinementalg, num):
    return "karman/"+DOLFIN_VERSION+"/%s/ref_%d/boundary.pvd"%(refinementalg,num)

