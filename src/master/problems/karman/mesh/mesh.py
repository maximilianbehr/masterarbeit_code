# -*- coding: utf-8 -*-
from dolfin import *
from dolfin import __version__
import os
import socket

"""constants for rectangular domain"""
rect         = {"x0":0,"x1":4,"y0":0,"y1":1}

"""constants for circualar domain"""
circ         = {"x0":1.5,"y0":0.5,"r":0.15,"fragments":15}

"""the resolution of the macro master mesh"""
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

"""macro master mesh files and boundary functions"""
if socket.gethostname() == "pc747":
    OUTPUTLOCATION 				= os.path.abspath("/scratch/behr/masters/src/master/data/karman")
elif socket.gethostname() == "pc800":
    pass
elif socket.gethostname() == "pc785":
    pass
elif socket.gethostname() == "pc633":
    pass
elif socket.gethostname() == "jack":
    OUTPUTLOCATION              = os.path.abspath("/Users/daniels/Documents/LiClipseWorkspace/master/src/master/data/karman")

KARMAN_MACRO_MESH_FILE                  = os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,"macro","macro.xml.gz")
KARMAN_MACRO_MESH_FILE_XDMF             = os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,"macro","macro.xdmf")
KARMAN_MACRO_BOUNDARY_FILE              = os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,"macro","boundary.xml.gz")
KARMAN_MACRO_BOUNDARY_FILE_XDMF         = os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,"macro","boundary.xdmf")
KARMAN_MACRO_MESH_PVD_FILE              = os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,"macro","macro.pvd")
KARMAN_MACRO_BOUNDARY_PVD_FILE          = os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,"macro","boundary.pvd")

"""output for refined mesh files"""
def KARMAN_REFN_MESH_FILE(refinementalg, num):
    return os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,refinementalg,"ref_%s"%num,"mesh.xml.gz")

def KARMAN_REFN_BOUNDARY_FILE(refinementalg, num):
    return os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,refinementalg,"ref_%s"%num,"boundary.xml.gz")

def KARMAN_REFN_MESH_FILE_XDMF(refinementalg, num):
    return os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,refinementalg,"ref_%s"%num,"mesh.xdmf")

def KARMAN_REFN_BOUNDARY_FILE_XDMF(refinementalg, num):
    return os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,refinementalg,"ref_%s"%num,"boundary.xdmf")

def KARMAN_REFN_MESH_PVD_FILE(refinementalg,num):
    return os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,refinementalg,"ref_%s"%num,"mesh.pvd")

def KARMAN_REFN_BOUNDARY_PVD_FILE(refinementalg, num):
    return os.path.join(OUTPUTLOCATION,DOLFIN_VERSION,refinementalg,"ref_%s"%num,"boundary.pvd")

class GammaLeft(SubDomain):
    """Left Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],rect["x0"])

class GammaLower(SubDomain):
    """Lower Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],rect["y0"])

class GammaRight(SubDomain):
    """Right Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],rect["x1"])

class GammaUpper(SubDomain):
    """Upper Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],rect["y1"])

class GammaBall(SubDomain):
    """Ball"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]
    tresh   = 1.5   #slightly larger use on boundary argument to ensure we get only the boundary of the Ball

    def inside(self,x,on_boundary):
        r = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        return on_boundary and r < self.tresh*self.radius

class GammaBallCtrlUpper(SubDomain):
    """Ball Upper"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]

    uppery  = centery + 6.0/8.0*radius
    lowery  = centery + 1.0/8.0*radius

    def inside(self,x,on_boundary):
        return GammaBall().inside(x, on_boundary) and between(x[1],(self.lowery,self.uppery)) and self.centerx<x[0]

class GammaBallCtrlLower(SubDomain):
    """Ball Lower"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]

    uppery  = centery-1.0/8.0*radius
    lowery  = centery-6.0/8.0*radius

    def inside(self,x,on_boundary):
        return GammaBall().inside(x, on_boundary) and between(x[1],(self.lowery,self.uppery)) and self.centerx<x[0]

class BallProjection(SubDomain):
    """Ball Projection"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]
    tresh   = 1.1   #slightly larger

    def inside(self,x,on_boundary):
        r       = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        return r < self.tresh*self.radius #dont use on_boundary argument to ensure that a slightly larger Ball is transformed

    def snap(self,x):
        r       = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        if r < self.tresh*self.radius:
            x[0]    = self.centerx + (self.radius/r)*(x[0]-self.centerx)
            x[1]    = self.centery + (self.radius/r)*(x[1]-self.centery)

def macrokarman():

    # use Constants related to the geometry for building CSG
    Rectangledomain     = Rectangle(rect["x0"],rect["y0"],rect["x1"],rect["y1"])
    Balldomain          = Circle(circ["x0"],circ["y0"],circ["r"],circ["fragments"])
    domain              = Rectangledomain - Balldomain

    # Mark boundary parts
    mesh                = Mesh(domain, initial_resolution)
    boundaryfunction    = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    boundaryfunction.set_all(GAMMAINNER)
    GammaLeft().mark(boundaryfunction,GAMMALEFT)
    GammaLower().mark(boundaryfunction,GAMMALOWER)
    GammaRight().mark(boundaryfunction,GAMMARIGHT)
    GammaUpper().mark(boundaryfunction,GAMMAUPPER)
    GammaBall().mark(boundaryfunction,GAMMABALL)
    GammaBallCtrlLower().mark(boundaryfunction,GAMMABALLCTRLUPPER)
    GammaBallCtrlUpper().mark(boundaryfunction,GAMMABALLCTRLLOWER)

    return mesh,boundaryfunction

def refinekarman(mesh):

    #refine given master
    refmesh = refine(mesh)

    #hole
    hole    = BallProjection()
    refmesh.snap_boundary(hole)

    #mark boundary parts
    boundaryfunction = MeshFunction("size_t", refmesh, refmesh.topology().dim()-1)
    boundaryfunction.set_all(GAMMAINNER)
    GammaLeft().mark(boundaryfunction,GAMMALEFT)
    GammaLower().mark(boundaryfunction,GAMMALOWER)
    GammaRight().mark(boundaryfunction,GAMMARIGHT)
    GammaUpper().mark(boundaryfunction,GAMMAUPPER)
    GammaBall().mark(boundaryfunction,GAMMABALL)
    GammaBallCtrlLower().mark(boundaryfunction,GAMMABALLCTRLUPPER)
    GammaBallCtrlUpper().mark(boundaryfunction,GAMMABALLCTRLLOWER)

    return refmesh,boundaryfunction


if __name__=="__main__":

    #for refalg in ["bisection", "iterative_bisection", "recursive_bisection", "regular_cut"]:
    for refalg in ["bisection"]:

        parameters["refinement_algorithm"]=refalg

        mesh,boundaryfunction = macrokarman()

        File(KARMAN_MACRO_MESH_FILE,"compressed")<<mesh
        File(KARMAN_MACRO_BOUNDARY_FILE,"compressed")<<boundaryfunction

        File(KARMAN_MACRO_MESH_FILE_XDMF,"compressed")<<mesh
        File(KARMAN_MACRO_BOUNDARY_FILE_XDMF,"compressed")<<boundaryfunction

        File(KARMAN_MACRO_MESH_PVD_FILE,"compressed")<<mesh
        File(KARMAN_MACRO_BOUNDARY_PVD_FILE,"compressed")<<boundaryfunction

        begin("Refinement with %s"%parameters["refinement_algorithm"])
        for ref in range(1,5):
            info("Level %d"%ref)
            mesh, boundaryfunction = refinekarman(mesh)

            File(KARMAN_REFN_MESH_FILE(parameters["refinement_algorithm"],ref),"compressed")<<mesh
            File(KARMAN_REFN_BOUNDARY_FILE(parameters["refinement_algorithm"],ref),"compressed")<<boundaryfunction

            File(KARMAN_REFN_MESH_FILE_XDMF(parameters["refinement_algorithm"],ref),"compressed")<<mesh
            File(KARMAN_REFN_BOUNDARY_FILE_XDMF(parameters["refinement_algorithm"],ref),"compressed")<<boundaryfunction

            File(KARMAN_REFN_MESH_PVD_FILE(parameters["refinement_algorithm"],ref),"compressed")<<mesh
            File(KARMAN_REFN_BOUNDARY_PVD_FILE(parameters["refinement_algorithm"],ref),"compressed")<<boundaryfunction

        end()



