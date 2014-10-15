# -*- coding: utf-8 -*-
import karman
from dolfin import __version__
from dolfin.mesh import refinement


"""the resolution of the macro karman mesh"""
initial_resolution                      = 5


"""dolfin version"""
DOLFIN_VERSION                          = __version__


"""macro karman mesh files and boundary functions"""
KARMAN_MACRO_MESH_FILE                  = "../karman_mesh/"+DOLFIN_VERSION+"/macro/macro.xml.gz"
KARMAN_MACRO_BOUNDARY_FILE              = "../karman_mesh/"+DOLFIN_VERSION+"/macro/boundary.xml.gz"
KARMAN_MACRO_MESH_PVD_FILE              = "../karman_mesh/"+DOLFIN_VERSION+"/macro/macro.pvd"
KARMAN_MACRO_BOUNDARY_PVD_FILE          = "../karman_mesh/"+DOLFIN_VERSION+"/macro/boundary.pvd"


"""output for refined mesh files"""
def KARMAN_REFN_MESH_FILE(refinementalg, num):
    return "../karman_mesh/"+DOLFIN_VERSION+"/%s/ref_%d/mesh.xml.gz"%(refinementalg,num)

def KARMAN_REFN_BOUNDARY_FILE(refinementalg, num):
    return "../karman_mesh/"+DOLFIN_VERSION+"/%s/ref_%d/boundary.xml.gz"%(refinementalg,num)

def KARMAN_REFN_MESH_PVD_FILE(refinementalg,num):
    return  "../karman_mesh/"+DOLFIN_VERSION+"/%s/ref_%d/mesh.pvd"%(refinementalg,num)

def KARMAN_REFN_BOUNDARY_PVD_FILE(refinementalg, num):
    return "../karman_mesh/"+DOLFIN_VERSION+"/%s/ref_%d/boundary.pvd"%(refinementalg,num)
